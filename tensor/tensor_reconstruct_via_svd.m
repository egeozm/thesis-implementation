function recon = tensor_reconstruct_via_svd(cpd_result, x, signal_type, components, params, fs)
%TENSOR_RECONSTRUCT_VIA_SVD Reconstruct tensor CPD via per-stride Hankel SVD.
%
%   recon = tensor_reconstruct_via_svd(cpd_result, x, signal_type, components, params, fs)
%
%   Reassembles a full Hankel approximation for each stride from the CPD
%   block tensor, performs an SVD on each reconstructed Hankel matrix, groups
%   those eigentriples using the requested grouping mode, then maps the
%   per-stride grouped reconstructions back to the original sample grid.

    if nargin < 6 || isempty(fs)
        fs = 1;
    end

    params = tensor_merge_params(params);
    params.fs = fs;
    signal_type = string(lower(string(signal_type)));
    x = x(:);
    recon = init_recon();

    if ~isfield(cpd_result, 'success') || ~cpd_result.success
        recon.failure_reason = "CPD result is not successful; cannot reconstruct components.";
        return;
    end
    if ~isfield(cpd_result, 'build_metadata') || ~isstruct(cpd_result.build_metadata) || ...
            ~isfield(cpd_result.build_metadata, 'stride_info') || ...
            ~isfield(cpd_result.build_metadata, 'block_info')
        recon.failure_reason = "Missing tensor build metadata in CPD result.";
        return;
    end

    build_metadata = cpd_result.build_metadata;
    stride_info = build_metadata.stride_info;
    block_info = build_metadata.block_info;
    if isempty(stride_info) || isempty(block_info)
        recon.failure_reason = "Tensor build metadata does not contain strides/blocks.";
        return;
    end

    [true_components, component_names] = truth_components(signal_type, components);
    num_components = size(true_components, 2);
    signal_length = numel(x);
    est_components = zeros(signal_length, num_components);
    component_weights = zeros(signal_length, num_components);

    stride_results = repmat(init_stride_result(), numel(stride_info), 1);
    num_successful_strides = 0;
    for i = 1:numel(stride_info)
        stride_result = reconstruct_stride(cpd_result, x, stride_info(i), block_info, ...
            signal_type, true_components, component_names, params);
        stride_results(i) = stride_result;

        if ~stride_result.success
            continue;
        end

        num_successful_strides = num_successful_strides + 1;
        sample_idx = stride_info(i).sample_indices(:);
        for g = 1:num_components
            if isempty(stride_result.est_components_decimated)
                continue;
            end
            component_vals = stride_result.est_components_decimated(:, g);
            if ~any(abs(component_vals) > 0)
                continue;
            end
            est_components(sample_idx, g) = est_components(sample_idx, g) + ...
                stride_result.combine_weight * component_vals;
            component_weights(sample_idx, g) = component_weights(sample_idx, g) + ...
                stride_result.combine_weight;
        end
    end

    if num_successful_strides == 0
        recon.failure_reason = "No stride produced a successful SVD-based reconstruction.";
        return;
    end

    for g = 1:num_components
        valid = component_weights(:, g) > 0;
        if any(valid)
            est_components(valid, g) = est_components(valid, g) ./ component_weights(valid, g);
        end
    end

    if any(sum(component_weights, 2) == 0)
        recon.failure_reason = "SVD-based reconstruction left uncovered samples.";
        return;
    end

    residual_full = x - sum(est_components, 2);

    recon.success = true;
    recon.failure_reason = "";
    recon.support_indices = (1:signal_length).';
    recon.support_length = signal_length;
    recon.reference_stride = 1;
    recon.reference_slice_index = 1;
    recon.support_counts = sum(component_weights, 2);
    recon.rank = cpd_result.rank;
    recon.est_components_support = est_components;
    recon.est_components_full = est_components;
    recon.residual_support = residual_full;
    recon.residual_full = residual_full;
    recon.slice_components = {};
    recon.component_names = component_names;
    recon.grouping = struct( ...
        'mode', string(params.grouping_mode), ...
        'per_stride', stride_results, ...
        'num_successful_strides', num_successful_strides);
    recon.per_stride = stride_results;
end

function stride_result = reconstruct_stride(cpd_result, x, stride_info, block_info, signal_type, true_components, component_names, params)
    stride_result = init_stride_result();
    stride_result.stride = stride_info.stride;
    stride_result.decimated_length = stride_info.decimated_length;
    stride_result.raw_rows = stride_info.raw_rows;
    stride_result.raw_cols = stride_info.raw_cols;
    stride_result.sample_indices = stride_info.sample_indices(:);

    [H_hat, H_count] = assemble_hankel_for_stride(cpd_result, stride_info, block_info);
    if any(H_count(:) == 0)
        stride_result.failure_reason = "Reassembled stride Hankel has uncovered entries.";
        return;
    end
    x_dec = x(stride_info.sample_indices);
    H_true = embed_hankel(x_dec, stride_info.raw_rows);
    stride_result.hankel_reconstructed = H_hat;
    stride_result.hankel_counts = H_count;
    stride_result.matrix_rel_error = norm(H_true(:) - H_hat(:)) / max(norm(H_true(:)), eps);
    stride_result.matrix_fit = 1 - stride_result.matrix_rel_error;
    stride_result.signal_rel_error = norm(x_dec - diagonal_average(H_hat, stride_info.decimated_length)) / max(norm(x_dec), eps);
    stride_result.combine_weight = max(stride_result.matrix_fit, 0.05);

    features = ssa_features_from_matrix(H_hat, stride_info.decimated_length, params.fs / stride_info.stride, params.group_r_max);
    stride_result.singular_values = features.singular_values;
    stride_result.features = features;

    [component_groups, grouping_info, success, failure_reason] = select_groups_for_stride( ...
        features, signal_type, true_components(stride_info.sample_indices, :), component_names, params);
    if ~success
        stride_result.failure_reason = failure_reason;
        return;
    end

    est_components_decimated = ssa_reconstruct_groups(features, component_groups);
    stride_result.success = true;
    stride_result.failure_reason = "";
    stride_result.component_groups = component_groups;
    stride_result.grouping = grouping_info;
    stride_result.est_components_decimated = est_components_decimated;
end

function [H_hat, H_count] = assemble_hankel_for_stride(cpd_result, stride_info, block_info)
    A = cpd_result.factors{1};
    B = cpd_result.factors{2};
    C = cpd_result.factors{3};
    lambda = cpd_result.lambda(:);

    H_sum = zeros(stride_info.raw_rows, stride_info.raw_cols);
    H_count = zeros(stride_info.raw_rows, stride_info.raw_cols);
    stride_block_ids = find([block_info.stride] == stride_info.stride);

    for block_id = stride_block_ids
        coeff = lambda .* C(block_id, :).';
        H_block = A * diag(coeff) * B.';
        row_range = block_info(block_id).row_start:block_info(block_id).row_end;
        col_range = block_info(block_id).col_start:block_info(block_id).col_end;
        H_sum(row_range, col_range) = H_sum(row_range, col_range) + H_block;
        H_count(row_range, col_range) = H_count(row_range, col_range) + 1;
    end

    H_hat = H_sum;
    valid = H_count > 0;
    H_hat(valid) = H_sum(valid) ./ H_count(valid);
end

function features = ssa_features_from_matrix(H_hat, N_decimated, fs, r_max)
    [U, S, V] = svd(H_hat, 'econ');
    singular_values = diag(S);
    rank_H = numel(singular_values);
    r_use = min(r_max, rank_H);

    elementary = zeros(N_decimated, r_use);
    dominant_freqs = nan(1, r_use);
    peak_values = nan(1, r_use);
    for r = 1:r_use
        Hr = S(r, r) * (U(:, r) * V(:, r).');
        elementary(:, r) = diagonal_average(Hr, N_decimated);
        [dominant_freqs(r), ~, peak_values(r)] = estimate_dominant_frequency(elementary(:, r), fs);
    end

    features = struct();
    features.x = diagonal_average(H_hat, N_decimated);
    features.N = N_decimated;
    features.L = size(H_hat, 1);
    features.K = size(H_hat, 2);
    features.rank = rank_H;
    features.U = U;
    features.S = S;
    features.V = V;
    features.singular_values = singular_values;
    features.r_max = r_use;
    features.elementary = elementary;
    features.dominant_freqs = dominant_freqs;
    features.peak_values = peak_values;
    features.fs = fs;
end

function [component_groups, grouping_info, success, failure_reason] = select_groups_for_stride(features, signal_type, true_components_dec, component_names, params)
    mode_name = lower(string(params.grouping_mode));
    switch mode_name
        case "oracle"
            [component_groups, grouping_info] = truth_guided_groups(features, signal_type, true_components_dec);
            success = ~isempty(component_groups) && all(~cellfun(@isempty, component_groups));
            if success
                failure_reason = "";
            else
                failure_reason = "truth-guided grouping returned empty groups";
            end
        case "benchmark"
            [component_groups, grouping_info] = benchmark_groups(features, signal_type, params);
            success = ~isempty(component_groups) && all(~cellfun(@isempty, component_groups));
            if success
                failure_reason = "";
            else
                failure_reason = "benchmark grouping returned empty groups";
            end
        case "auto"
            grouping_auto = auto_group_ssa(features, auto_group_params(params));
            [component_groups, success, failure_reason] = groups_for_signal(signal_type, grouping_auto);
            grouping_info = grouping_auto;
            if success
                grouping_info.component_names = component_names;
            end
        otherwise
            success = false;
            failure_reason = sprintf('Unsupported SVD reconstruction grouping mode %s.', params.grouping_mode);
            component_groups = {};
            grouping_info = struct();
    end
end

function params_out = auto_group_params(params)
    params_out = struct();
    params_out.r_max = params.group_r_max;
    params_out.sv_ratio_thresh = params.group_sv_ratio_thresh;
    params_out.freq_tol = params.group_freq_tol;
    params_out.min_osc_freq = params.group_min_osc_freq;
end

function [component_groups, grouping] = truth_guided_groups(features, signal_type, true_components_dec)
    corr_to_truth = zeros(size(true_components_dec, 2), features.r_max);
    for r = 1:features.r_max
        for c = 1:size(true_components_dec, 2)
            corr_to_truth(c, r) = abs(pearson_raw(features.elementary(:, r), true_components_dec(:, c)));
        end
    end

    switch lower(string(signal_type))
        case "multiscale"
            used = false(1, features.r_max);
            slow_idx = select_oscillatory_pair(corr_to_truth(2, :), used);
            used(slow_idx) = true;
            fast_idx = select_oscillatory_pair(corr_to_truth(3, :), used);
            used(fast_idx) = true;
            trend_idx = select_trend_group(corr_to_truth(1, :), used);
            component_groups = {sort(trend_idx), sort(slow_idx), sort(fast_idx)};
        case "close"
            used = false(1, features.r_max);
            comp1_idx = select_oscillatory_pair(corr_to_truth(1, :), used);
            used(comp1_idx) = true;
            comp2_idx = select_oscillatory_pair(corr_to_truth(2, :), used);
            component_groups = {sort(comp1_idx), sort(comp2_idx)};
        otherwise
            error('tensor_reconstruct_via_svd:UnknownSignal', 'Unknown signal type %s.', signal_type);
    end

    grouping = struct();
    grouping.mode = "oracle";
    grouping.corr_to_truth = corr_to_truth;
    grouping.component_groups = component_groups;
end

function [component_groups, grouping] = benchmark_groups(features, signal_type, params)
    specs = benchmark_specs(signal_type, features.N, features.fs);
    num_groups = numel(specs);
    component_groups = cell(1, num_groups);
    scores = zeros(features.r_max, num_groups);

    for r = 1:features.r_max
        term = features.elementary(:, r);
        dom_freq = features.dominant_freqs(r);
        for g = 1:num_groups
            scores(r, g) = benchmark_score(term, dom_freq, specs(g), params);
        end
        [best_score, best_idx] = max(scores(r, :));
        if best_score >= params.benchmark_assign_min_score
            component_groups{best_idx} = [component_groups{best_idx}, r]; %#ok<AGROW>
        end
    end

    for g = 1:num_groups
        if isempty(component_groups{g})
            [~, fallback_idx] = max(scores(:, g));
            component_groups{g} = fallback_idx;
        end
        component_groups{g} = unique(component_groups{g}, 'stable');
    end

    grouping = struct();
    grouping.mode = "benchmark";
    grouping.scores = scores;
    grouping.specs = specs;
    grouping.component_groups = component_groups;
end

function specs = benchmark_specs(signal_type, N, fs)
    n = (1:N).';
    n_norm = n / N;
    switch lower(string(signal_type))
        case "multiscale"
            specs = [ ...
                make_spec("trend", nan, [ones(N, 1), n_norm, n_norm.^2], true), ...
                make_spec("slow", 0.01, sinusoid_basis(n, 0.01, fs), false), ...
                make_spec("fast", 0.08, sinusoid_basis(n, 0.08, fs), false) ...
                ];
        case "close"
            specs = [ ...
                make_spec("comp1", 0.05, sinusoid_basis(n, 0.05, fs), false), ...
                make_spec("comp2", 0.055, sinusoid_basis(n, 0.055, fs), false) ...
                ];
        otherwise
            error('tensor_reconstruct_via_svd:UnknownSignal', 'Unknown signal type %s.', signal_type);
    end
end

function spec = make_spec(name, target_freq, basis, is_trend)
    spec = struct();
    spec.name = char(name);
    spec.target_freq = target_freq;
    spec.basis = orthonormal_basis(basis);
    spec.is_trend = is_trend;
end

function basis = sinusoid_basis(n, freq, fs)
    omega = 2 * pi * freq / max(fs, eps);
    basis = [sin(omega * n), cos(omega * n)];
end

function score = benchmark_score(term, dominant_freq, spec, params)
    projected = project_onto_basis(term, spec.basis);
    proj_ratio = sum(projected.^2) / max(sum(term.^2), eps);
    if spec.is_trend
        if dominant_freq < params.benchmark_trend_freq_max
            freq_weight = 1;
        else
            freq_weight = 0.1;
        end
    else
        if isnan(dominant_freq) || dominant_freq <= 0
            freq_weight = 0.1;
        else
            freq_diff = abs(dominant_freq - spec.target_freq);
            freq_weight = 1 / (1 + freq_diff / max(params.benchmark_freq_tol, eps));
        end
    end
    score = proj_ratio * freq_weight;
end

function Q = orthonormal_basis(B)
    [Q, ~] = qr(B, 0);
end

function projected = project_onto_basis(x, basis)
    projected = basis * (basis' * x(:));
end

function [component_groups, success, failure_reason] = groups_for_signal(signal_type, grouping)
    pairs = grouping.oscillatory_pairs;
    switch lower(string(signal_type))
        case "multiscale"
            if numel(pairs) < 2
                component_groups = {};
                success = false;
                failure_reason = "fewer than two oscillatory pairs";
                return;
            end
            if isempty(grouping.trend_idx)
                component_groups = {};
                success = false;
                failure_reason = "empty trend group";
                return;
            end
            component_groups = {grouping.trend_idx, pairs{1}, pairs{2}};
            success = true;
            failure_reason = "";
        case "close"
            if numel(pairs) < 2
                component_groups = {};
                success = false;
                failure_reason = "fewer than two oscillatory pairs";
                return;
            end
            component_groups = {pairs{1}, pairs{2}};
            success = true;
            failure_reason = "";
        otherwise
            error('tensor_reconstruct_via_svd:UnknownSignal', 'Unknown signal type %s.', signal_type);
    end
end

function idx = select_oscillatory_pair(scores, used)
    available = find(~used);
    if isempty(available)
        idx = [];
        return;
    end
    [~, best_pos] = max(scores(available));
    seed = available(best_pos);
    remaining = available(available ~= seed);
    neighbors = remaining(abs(remaining - seed) == 1);
    if ~isempty(neighbors)
        [~, partner_pos] = max(scores(neighbors));
        partner = neighbors(partner_pos);
    elseif ~isempty(remaining)
        [~, partner_pos] = max(scores(remaining));
        partner = remaining(partner_pos);
    else
        partner = [];
    end
    idx = unique(sort([seed, partner]));
end

function idx = select_trend_group(scores, used)
    available = find(~used);
    if isempty(available)
        idx = [];
        return;
    end
    available_scores = scores(available);
    [sorted_scores, order] = sort(available_scores, 'descend');
    threshold = max(0.15 * sorted_scores(1), 0.05);
    keep_mask = sorted_scores >= threshold;
    keep = available(order(keep_mask));
    if isempty(keep)
        keep = available(order(1));
    end
    idx = sort(keep(1:min(4, numel(keep))));
end

function [true_components, component_names] = truth_components(signal_type, components)
    switch lower(string(signal_type))
        case "multiscale"
            true_components = [components.trend, components.slow, components.fast];
            component_names = {"trend", "slow", "fast"};
        case "close"
            true_components = [components.comp1, components.comp2];
            component_names = {"comp1", "comp2"};
        otherwise
            error('tensor_reconstruct_via_svd:UnknownSignal', 'Unknown signal type %s.', signal_type);
    end
end

function stride_result = init_stride_result()
    stride_result = struct( ...
        'success', false, ...
        'failure_reason', "", ...
        'stride', nan, ...
        'decimated_length', nan, ...
        'raw_rows', nan, ...
        'raw_cols', nan, ...
        'sample_indices', [], ...
        'hankel_reconstructed', [], ...
        'hankel_counts', [], ...
        'matrix_rel_error', nan, ...
        'matrix_fit', nan, ...
        'signal_rel_error', nan, ...
        'combine_weight', nan, ...
        'features', struct(), ...
        'singular_values', [], ...
        'component_groups', {{}}, ...
        'grouping', struct(), ...
        'est_components_decimated', []);
end

function recon = init_recon()
    recon = struct();
    recon.success = false;
    recon.failure_reason = "";
    recon.support_indices = [];
    recon.support_length = nan;
    recon.reference_stride = nan;
    recon.reference_slice_index = nan;
    recon.support_counts = [];
    recon.rank = nan;
    recon.est_components_support = [];
    recon.est_components_full = [];
    recon.residual_support = [];
    recon.residual_full = [];
    recon.slice_components = {};
    recon.component_names = {};
    recon.grouping = struct();
    recon.per_stride = struct([]);
end

function c = pearson_raw(a, b)
    a = a(:) - mean(a);
    b = b(:) - mean(b);
    den = norm(a) * norm(b);
    if den == 0
        c = 0;
    else
        c = (a.' * b) / den;
    end
end
