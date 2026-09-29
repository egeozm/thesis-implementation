function recon = tensor_reconstruct_components(cpd_result, x)
%TENSOR_RECONSTRUCT_COMPONENTS Reconstruct full-length 1D terms from CPD blocks.
%
%   recon = tensor_reconstruct_components(cpd_result, x)
%
%   Uses the build metadata stored in tensor_run_cpd output to map each CP
%   rank-1 term back to time-domain samples by inverse-mapping each tensor
%   block and overlap-adding those contributions onto the original grid.

    x = x(:);
    recon = init_recon();

    if ~isfield(cpd_result, 'success') || ~cpd_result.success
        recon.failure_reason = "CPD result is not successful; cannot reconstruct components.";
        return;
    end
    if ~isfield(cpd_result, 'build_metadata') || ~isstruct(cpd_result.build_metadata) || ...
            ~isfield(cpd_result.build_metadata, 'block_info')
        recon.failure_reason = "Missing tensor build metadata in CPD result.";
        return;
    end

    build_metadata = cpd_result.build_metadata;
    block_info = build_metadata.block_info;
    if isempty(block_info)
        recon.failure_reason = "Tensor build metadata does not contain any blocks.";
        return;
    end

    U = cpd_result.factors;
    lambda = cpd_result.lambda(:);
    if numel(U) ~= 3
        recon.failure_reason = "Current reconstruction expects a third-order CPD.";
        return;
    end

    rank_k = numel(lambda);
    signal_length = numel(x);
    est_components = zeros(signal_length, rank_k);
    sample_weights = zeros(signal_length, 1);
    block_weights = block_weight_vector(build_metadata.params.recon_weighting, block_info(1));

    for b = 1:numel(block_info)
        idx = block_info(b).retained_original_indices(:);
        if any(idx < 1) || any(idx > signal_length)
            recon.failure_reason = "Block metadata contains indices outside the original signal range.";
            return;
        end
        sample_weights(idx) = sample_weights(idx) + block_weights;
    end

    if any(sample_weights == 0)
        recon.failure_reason = "Full-length reconstruction left uncovered time indices.";
        return;
    end

    A = U{1};
    B = U{2};
    block_factor = U{3};
    for r = 1:rank_k
        accum = zeros(signal_length, 1);
        per_block = cell(numel(block_info), 1);
        for b = 1:numel(block_info)
            coeff = lambda(r) * block_factor(b, r);
            H_rs = coeff * (A(:, r) * B(:, r).');
            y_slice = diagonal_average(H_rs, size(H_rs, 1) + size(H_rs, 2) - 1);
            idx = block_info(b).retained_original_indices(:);
            weighted_vals = y_slice .* block_weights;
            accum(idx) = accum(idx) + weighted_vals;
            per_block{b} = struct('indices', idx, 'values', y_slice);
        end
        est_components(:, r) = accum ./ sample_weights;
        slice_components{r} = per_block; %#ok<AGROW>
    end

    residual_full = x - sum(est_components, 2);

    recon.success = true;
    recon.failure_reason = "";
    recon.support_indices = (1:signal_length).';
    recon.support_length = signal_length;
    recon.reference_stride = 1;
    recon.reference_slice_index = 1;
    recon.support_counts = sample_weights;
    recon.rank = rank_k;
    recon.est_components_support = est_components;
    recon.est_components_full = est_components;
    recon.residual_support = residual_full;
    recon.residual_full = residual_full;
    recon.slice_components = slice_components;
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
end

function w = block_weight_vector(mode_name, block)
    support_len = numel(block.retained_original_indices);
    switch lower(string(mode_name))
        case "uniform"
            w = ones(support_len, 1);
        otherwise
            error('tensor_reconstruct_components:UnknownWeighting', ...
                'Unknown reconstruction weighting mode %s.', mode_name);
    end
end
