function outcome = evaluate_tensor_cpd(signal_type, x, components, params, fs)
%EVALUATE_TENSOR_CPD Evaluate grouped CPD reconstructions against synthetic truth.
%
%   outcome = evaluate_tensor_cpd(signal_type, x, components, params, fs)
%
%   Reconstructs CP rank-1 terms onto the full original grid, groups those
%   terms according to params.grouping_mode, appends a residual column, and
%   evaluates the matched estimates against the full-length synthetic truth.

    if nargin < 5 || isempty(fs)
        fs = 1;
    end

    signal_type = string(lower(string(signal_type)));
    params = tensor_merge_params(params);
    params.fs = fs;

    tensor_out = tensor_build_decimated_hankel(x, params);
    cpd_out = tensor_run_cpd(tensor_out, params);

    outcome = struct();
    outcome.signal_type = signal_type;
    outcome.tensor_output = tensor_out;
    outcome.cpd_result = cpd_out;
    outcome.component_names = {};
    outcome.est_component_names = {};
    outcome.success = false;
    outcome.failure_reason = "";
    outcome.reconstruction = struct();
    outcome.grouping = struct();
    outcome.grouping_mode = string(params.grouping_mode);
    outcome.reconstruction_mode = string(params.reconstruction_mode);
    outcome.support_indices = [];
    outcome.support_length = nan;
    outcome.est_components = [];
    outcome.metrics = empty_metrics(0);
    outcome.mean_nmse = nan;
    outcome.mean_abs_rho = nan;
    outcome.mean_peak_error = nan;
    outcome.fit = cpd_out.fit;
    outcome.rel_error = cpd_out.rel_error;
    outcome.score = inf;

    if ~cpd_out.success
        outcome.failure_reason = cpd_out.failure_reason;
        return;
    end

    [true_components, component_names] = truth_components(signal_type, components);
    use_cpd_reconstruction = strcmpi(params.reconstruction_mode, 'cpd');
    use_svd_reconstruction = strcmpi(params.reconstruction_mode, 'svd') && ...
        ~strcmpi(params.grouping_mode, 'none');

    if use_cpd_reconstruction
        recon = tensor_reconstruct_via_cpd(cpd_out, x, params);
        outcome.reconstruction = recon;
        outcome.support_indices = recon.support_indices;
        outcome.support_length = recon.support_length;
        if ~recon.success
            outcome.failure_reason = recon.failure_reason;
            return;
        end
        [est_components, grouping, est_component_names] = apply_grouping( ...
            recon.est_components_full, true_components, cpd_out.lambda, params, signal_type);
        [est_components, grouping] = refine_grouped_components_from_residual( ...
            x(:), est_components, grouping, signal_type, params);
        residual_full = x(:) - sum(est_components, 2);
    elseif use_svd_reconstruction
        recon = tensor_reconstruct_via_svd(cpd_out, x, signal_type, components, params, fs);
        outcome.reconstruction = recon;
        outcome.support_indices = recon.support_indices;
        outcome.support_length = recon.support_length;
        if ~recon.success
            outcome.failure_reason = recon.failure_reason;
            return;
        end
        est_components = recon.est_components_full;
        residual_full = recon.residual_full;
        grouping = recon.grouping;
        est_component_names = recon.component_names;
    else
        recon = tensor_reconstruct_components(cpd_out, x);
        outcome.reconstruction = recon;
        outcome.support_indices = recon.support_indices;
        outcome.support_length = recon.support_length;
        if ~recon.success
            outcome.failure_reason = recon.failure_reason;
            return;
        end
        [est_components, grouping, est_component_names] = apply_grouping( ...
            recon.est_components_full, true_components, cpd_out.lambda, params, signal_type);
        residual_full = x(:) - sum(est_components, 2);
    end

    est_support = [est_components, residual_full];
    [assignment, ~] = match_components(true_components, est_support);
    metrics = compute_metrics(true_components, est_support, assignment, fs);

    outcome.component_names = component_names;
    outcome.est_component_names = [est_component_names, {"residual"}];
    outcome.success = true;
    outcome.failure_reason = "";
    outcome.grouping = grouping;
    outcome.grouping_mode = string(grouping.mode);
    outcome.est_components = est_support;
    outcome.metrics = metrics;
    outcome.mean_nmse = mean(metrics.nmse, 'omitnan');
    outcome.mean_abs_rho = mean(abs(metrics.pearson), 'omitnan');
    outcome.mean_peak_error = mean(metrics.peak_freq_error, 'omitnan');
    outcome.score = outcome.mean_nmse;
end

function [est_components, grouping, est_component_names] = apply_grouping(rank1_components, true_components, lambda, params, signal_type)
    mode_name = lower(string(params.grouping_mode));
    switch mode_name
        case "none"
            est_components = rank1_components;
            grouping = build_none_grouping(lambda, size(rank1_components, 2));
            est_component_names = component_labels("term", size(est_components, 2));
        case "oracle"
            [est_components, grouping] = group_cp_terms_oracle(rank1_components, true_components, params);
            est_component_names = component_labels("oracle_group", size(est_components, 2));
        case "benchmark"
            [est_components, grouping, est_component_names] = ...
                group_cp_terms_benchmark(rank1_components, lambda, signal_type, params);
        case "auto"
            [est_components, grouping] = group_cp_terms_auto(rank1_components, lambda, params);
            est_component_names = component_labels("group", size(est_components, 2));
        otherwise
            error('evaluate_tensor_cpd:UnknownGroupingMode', ...
                'Unsupported grouping mode "%s".', params.grouping_mode);
    end
end

function grouping = build_none_grouping(lambda, num_terms)
    grouping = struct();
    grouping.mode = "none";
    grouping.lambda = lambda(:);
    grouping.num_grouped_components = num_terms;
    grouping.clusters = struct([]);
end

function labels = component_labels(prefix, count)
    labels = cell(1, count);
    for k = 1:count
        labels{k} = char(prefix + string(k));
    end
end

function [true_components, component_names] = truth_components(signal_type, components)
    switch signal_type
        case "multiscale"
            true_components = [components.trend, components.slow, components.fast];
            component_names = {"trend", "slow", "fast"};
        case "close"
            true_components = [components.comp1, components.comp2];
            component_names = {"comp1", "comp2"};
        otherwise
            error('evaluate_tensor_cpd:UnknownSignal', 'Unknown signal type %s', signal_type);
    end
end

function metrics = empty_metrics(num_components)
    metrics = struct();
    metrics.nmse = nan(num_components, 1);
    metrics.pearson = nan(num_components, 1);
    metrics.peak_freq_error = nan(num_components, 1);
    metrics.leakage = nan(max(num_components, 1), 1);
    metrics.dominant_freq_true = nan(num_components, 1);
    metrics.dominant_freq_est = nan(num_components, 1);
    metrics.true_idx = (1:num_components).';
    metrics.est_idx = nan(num_components, 1);
end

function [est_components, grouping] = refine_grouped_components_from_residual(x, est_components, grouping, signal_type, params)
    signal_type = lower(string(signal_type));
    if signal_type ~= "multiscale" || size(est_components, 2) < 3
        return;
    end

    % In the pure-CPD path the oscillatory components are often recovered well,
    % while the low-order trend remains in the residual. Reassign that smooth
    % residual energy to the dedicated trend component.
    residual = x(:) - sum(est_components, 2);
    trend_basis = trend_polynomial_basis(numel(x));
    trend_update = trend_basis * (trend_basis' * residual);
    est_components(:, 1) = est_components(:, 1) + trend_update;

    if isstruct(grouping)
        grouping.trend_refine_mode = "residual_polynomial_projection";
        grouping.trend_refine_energy = norm(trend_update);
    end
end

function Q = trend_polynomial_basis(N)
    n = (1:N).';
    n_norm = n / max(N, 1);
    [Q, ~] = qr([ones(N, 1), n_norm, n_norm.^2], 0);
end
