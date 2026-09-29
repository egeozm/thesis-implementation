function recon = tensor_reconstruct_via_cpd(cpd_result, x, params)
%TENSOR_RECONSTRUCT_VIA_CPD Reconstruct full-length CP terms without post-CPD SVD.
%
%   recon = tensor_reconstruct_via_cpd(cpd_result, x, params)
%
%   Each CP rank-1 term is first reassembled into a full Hankel matrix within
%   each stride by averaging overlapping CP-generated blocks at their native
%   row/column locations. That stride-level term Hankel is then diagonal-
%   averaged to the decimated time grid, and the resulting per-stride term
%   estimates are fused back onto the original grid. This keeps the inverse
%   path fully CPD-based while avoiding the old blockwise 1D overlap-add.

    if nargin < 3
        params = struct();
    end

    params = tensor_merge_params(params);
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

    U = cpd_result.factors;
    lambda = cpd_result.lambda(:);
    if numel(U) ~= 3
        recon.failure_reason = "Current reconstruction expects a third-order CPD.";
        return;
    end

    A = U{1};
    B = U{2};
    block_factor = U{3};
    rank_k = numel(lambda);
    signal_length = numel(x);
    est_components = zeros(signal_length, rank_k);
    component_weights = zeros(signal_length, rank_k);
    stride_results = repmat(init_stride_result(), numel(stride_info), 1);
    num_successful_strides = 0;

    for i = 1:numel(stride_info)
        stride_result = reconstruct_stride( ...
            i, stride_info(i), block_info, x, A, B, block_factor, lambda);
        stride_results(i) = stride_result;
        if ~stride_result.success
            continue;
        end

        num_successful_strides = num_successful_strides + 1;
        sample_idx = stride_result.sample_indices(:);
        for r = 1:rank_k
            term_vals = stride_result.term_components_decimated(:, r);
            est_components(sample_idx, r) = est_components(sample_idx, r) + ...
                stride_result.combine_weight * term_vals;
            component_weights(sample_idx, r) = component_weights(sample_idx, r) + ...
                stride_result.combine_weight;
        end
    end

    if num_successful_strides == 0
        recon.failure_reason = "No stride produced a successful pure-CPD reconstruction.";
        recon.per_stride = stride_results;
        return;
    end

    for r = 1:rank_k
        valid = component_weights(:, r) > 0;
        if any(valid)
            est_components(valid, r) = est_components(valid, r) ./ component_weights(valid, r);
        end
    end

    if any(sum(component_weights, 2) == 0)
        recon.failure_reason = "Pure-CPD reconstruction left uncovered time indices.";
        recon.per_stride = stride_results;
        return;
    end

    residual_full = x - sum(est_components, 2);

    recon.success = true;
    recon.failure_reason = "";
    recon.support_indices = (1:signal_length).';
    recon.support_length = signal_length;
    recon.reference_stride = stride_results(find([stride_results.success], 1, 'first')).stride;
    recon.reference_slice_index = 1;
    recon.support_counts = sum(component_weights, 2);
    recon.rank = rank_k;
    recon.est_components_support = est_components;
    recon.est_components_full = est_components;
    recon.residual_support = residual_full;
    recon.residual_full = residual_full;
    recon.slice_components = {};
    recon.per_stride = stride_results;
    recon.grouping = struct( ...
        'mode', "cpd", ...
        'num_successful_strides', num_successful_strides, ...
        'per_stride', stride_results);
end

function stride_result = reconstruct_stride(stride_idx, stride_meta, block_info, x, A, B, block_factor, lambda)
    stride_result = init_stride_result();
    stride_result.stride = stride_meta.stride;
    stride_result.sample_indices = stride_meta.sample_indices(:);
    stride_result.decimated_length = stride_meta.decimated_length;
    stride_result.num_blocks = stride_meta.num_blocks;

    stride_block_ids = find([block_info.stride_index] == stride_idx);
    if isempty(stride_block_ids)
        stride_result.failure_reason = "No blocks found for stride.";
        return;
    end

    dec_len = stride_meta.decimated_length;
    rank_k = numel(lambda);
    count_matrix = zeros(stride_meta.raw_rows, stride_meta.raw_cols);
    hankel_terms = zeros(stride_meta.raw_rows, stride_meta.raw_cols, rank_k);

    for bi = 1:numel(stride_block_ids)
        block_id = stride_block_ids(bi);
        row_range = block_info(block_id).row_start:block_info(block_id).row_end;
        col_range = block_info(block_id).col_start:block_info(block_id).col_end;
        count_matrix(row_range, col_range) = count_matrix(row_range, col_range) + 1;

        coeff_all = lambda(:).' .* block_factor(block_id, :);
        for r = 1:rank_k
            H_r = coeff_all(r) * (A(:, r) * B(:, r).');
            hankel_terms(row_range, col_range, r) = hankel_terms(row_range, col_range, r) + H_r;
        end
    end

    if any(count_matrix(:) == 0)
        stride_result.failure_reason = "Stride Hankel reassembly left uncovered entries.";
        return;
    end

    term_components = zeros(dec_len, rank_k);
    valid = count_matrix > 0;
    for r = 1:rank_k
        H_r = hankel_terms(:, :, r);
        H_r(valid) = H_r(valid) ./ count_matrix(valid);
        term_components(:, r) = diagonal_average(H_r, dec_len);
    end
    total_decimated = sum(term_components, 2);

    x_dec = x(stride_result.sample_indices);
    signal_rel_error = norm(x_dec - total_decimated) / max(norm(x_dec), eps);

    stride_result.success = true;
    stride_result.failure_reason = "";
    stride_result.signal_rel_error = signal_rel_error;
    stride_result.combine_weight = max(1 - signal_rel_error, 0);
    stride_result.coverage_fraction = mean(valid(:));
    stride_result.support_weights = count_matrix;
    stride_result.total_decimated_reconstruction = total_decimated;
    stride_result.term_components_decimated = term_components;
end

function stride_result = init_stride_result()
    stride_result = struct( ...
        'success', false, ...
        'failure_reason', "", ...
        'stride', nan, ...
        'sample_indices', [], ...
        'decimated_length', nan, ...
        'num_blocks', nan, ...
        'signal_rel_error', nan, ...
        'combine_weight', nan, ...
        'coverage_fraction', nan, ...
        'support_weights', [], ...
        'total_decimated_reconstruction', [], ...
        'term_components_decimated', []);
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
    recon.grouping = struct();
    recon.per_stride = struct([]);
end

