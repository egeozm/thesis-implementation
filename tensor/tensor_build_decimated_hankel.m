function out = tensor_build_decimated_hankel(x, params)
%TENSOR_BUILD_DECIMATED_HANKEL Build a tensor from decimated Hankel slices.
%
%   out = tensor_build_decimated_hankel(x, params)
%
%   The output struct contains:
%     tensor         - L x K x B numeric tensor of Hankel blocks
%     tensor_size    - [L K B]
%     common_shape   - struct with rows, cols, num_blocks, num_strides
%     valid_strides  - strides used in the block tensor
%     stride_info    - per-stride block metadata
%     block_info     - per-block metadata for full-length reconstruction
%     skipped_slices - invalid stride settings with reasons

    params = tensor_merge_params(params);
    x = x(:);
    validation = tensor_validate_params(x, params);

    block_rows = validation.block_num_rows;
    block_cols = validation.block_num_cols;
    stride_info = validation.candidate_slices;
    num_blocks = sum([stride_info.num_blocks]);
    T = zeros(block_rows, block_cols, num_blocks);
    block_info = repmat(init_block_info(), num_blocks, 1);

    block_id = 0;
    for i = 1:numel(stride_info)
        candidate = stride_info(i);
        x_dec = x(candidate.sample_indices);
        H = embed_hankel(x_dec, candidate.raw_rows);
        support_len = block_rows + block_cols - 1;
        block_counter = 0;

        for rr = 1:numel(candidate.row_starts)
            row_start = candidate.row_starts(rr);
            row_end = row_start + block_rows - 1;
            for cc = 1:numel(candidate.col_starts)
                col_start = candidate.col_starts(cc);
                col_end = col_start + block_cols - 1;
                block_counter = block_counter + 1;
                block_id = block_id + 1;

                H_block = H(row_start:row_end, col_start:col_end);
                decimated_support = (row_start + col_start - 1):(row_start + col_start + support_len - 2);
                original_indices = candidate.sample_indices(decimated_support);

                T(:, :, block_id) = H_block;
                block_info(block_id) = struct( ...
                    'block_id', block_id, ...
                    'stride', candidate.stride, ...
                    'stride_index', i, ...
                    'block_index_within_stride', block_counter, ...
                    'row_block_index', rr, ...
                    'col_block_index', cc, ...
                    'row_start', row_start, ...
                    'row_end', row_end, ...
                    'col_start', col_start, ...
                    'col_end', col_end, ...
                    'retained_decimated_support', decimated_support(:), ...
                    'retained_original_indices', original_indices(:), ...
                    'first_original_index', original_indices(1), ...
                    'last_original_index', original_indices(end));
            end
        end
    end

    coverage = build_coverage(block_info, numel(x));
    out = struct();
    out.tensor = T;
    out.tensor_size = size(T);
    out.params = params;
    out.signal_length = numel(x);
    out.valid_strides = validation.valid_strides;
    out.num_slices = num_blocks;
    out.num_strides = numel(stride_info);
    out.num_blocks = num_blocks;
    out.common_shape = struct('rows', block_rows, 'cols', block_cols, ...
        'num_blocks', num_blocks, 'num_strides', numel(stride_info));
    out.slice_info = stride_info;
    out.stride_info = stride_info;
    out.block_info = block_info;
    out.skipped_slices = validation.skipped_slices;
    out.tensor_norm = norm(T(:));
    out.coverage = coverage;
end

function info = init_block_info()
    info = struct( ...
        'block_id', nan, ...
        'stride', nan, ...
        'stride_index', nan, ...
        'block_index_within_stride', nan, ...
        'row_block_index', nan, ...
        'col_block_index', nan, ...
        'row_start', nan, ...
        'row_end', nan, ...
        'col_start', nan, ...
        'col_end', nan, ...
        'retained_decimated_support', [], ...
        'retained_original_indices', [], ...
        'first_original_index', nan, ...
        'last_original_index', nan);
end

function coverage = build_coverage(block_info, signal_length)
    coverage = struct();
    coverage.num_blocks = numel(block_info);
    coverage.min_first_original_index = min([block_info.first_original_index]);
    coverage.max_last_original_index = max([block_info.last_original_index]);
    coverage.support_lengths = arrayfun(@(b) numel(b.retained_original_indices), block_info).';
    counts = zeros(signal_length, 1);
    for i = 1:numel(block_info)
        counts(block_info(i).retained_original_indices) = counts(block_info(i).retained_original_indices) + 1;
    end
    coverage.covered_counts = counts;
    coverage.uncovered_indices = find(counts == 0);
    coverage.coverage_fraction = mean(counts > 0);
end
