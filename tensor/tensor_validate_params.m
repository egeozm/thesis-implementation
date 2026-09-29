function info = tensor_validate_params(x, params)
%TENSOR_VALIDATE_PARAMS Validate decimated Hankel tensor settings.
%
%   info = tensor_validate_params(x, params)
%
%   Returns a struct describing valid slices, skipped slices, and the
%   full-coverage block geometry that will be used for tensor stacking.

    params = tensor_merge_params(params);
    x = x(:);
    N = numel(x);

    if N < 4
        error('tensor_validate_params:ShortSignal', ...
            'Signal must have at least 4 samples (got N=%d).', N);
    end
    if ~isscalar(params.base_window_length) || params.base_window_length < 2
        error('tensor_validate_params:BadWindow', ...
            'base_window_length must be a scalar >= 2.');
    end
    if ~isscalar(params.start_index) || params.start_index < 1 || params.start_index > N || ...
            params.start_index ~= floor(params.start_index)
        error('tensor_validate_params:BadStartIndex', ...
            'start_index must be an integer in [1, N].');
    end
    if ~isscalar(params.cpd_rank) || params.cpd_rank < 1 || params.cpd_rank ~= floor(params.cpd_rank)
        error('tensor_validate_params:BadRank', ...
            'cpd_rank must be a positive integer.');
    end
    if ~isscalar(params.min_num_slices) || params.min_num_slices < 1 || ...
            params.min_num_slices ~= floor(params.min_num_slices)
        error('tensor_validate_params:BadMinSlices', ...
            'min_num_slices must be a positive integer.');
    end
    if ~isscalar(params.min_blocks_per_stride) || params.min_blocks_per_stride < 1 || ...
            params.min_blocks_per_stride ~= floor(params.min_blocks_per_stride)
        error('tensor_validate_params:BadMinBlocks', ...
            'min_blocks_per_stride must be a positive integer.');
    end
    reconstruction_mode = lower(string(params.reconstruction_mode));
    if reconstruction_mode ~= "svd" && reconstruction_mode ~= "blockwise" && reconstruction_mode ~= "cpd"
        error('tensor_validate_params:BadReconstructionMode', ...
            'reconstruction_mode must be one of "cpd", "svd", or "blockwise".');
    end

    strides = unique(params.strides(:).');
    if isempty(strides)
        error('tensor_validate_params:EmptyStrides', 'params.strides must not be empty.');
    end
    if any(~isfinite(strides)) || any(strides < 1) || any(strides ~= floor(strides))
        error('tensor_validate_params:BadStrides', ...
            'params.strides must contain positive integers only.');
    end

    candidate_slices = repmat(init_slice_info(), 0, 1);
    skipped_slices = repmat(init_skipped_slice(), 0, 1);

    for i = 1:numel(strides)
        stride = strides(i);
        sample_indices = (params.start_index:stride:N).';
        decimated_length = numel(sample_indices);
        raw_rows = floor((params.base_window_length - 1) / stride) + 1;
        raw_cols = decimated_length - raw_rows + 1;

        if raw_rows < 2
            skipped_slices(end + 1) = skipped_record(stride, decimated_length, ...
                raw_rows, max(raw_cols, 0), ...
                'window too short after decimation; need at least 2 rows');
            continue;
        end
        if raw_rows >= decimated_length
            skipped_slices(end + 1) = skipped_record(stride, decimated_length, ...
                raw_rows, max(raw_cols, 0), ...
                'decimated signal too short for requested window');
            continue;
        end
        if raw_cols < 2
            skipped_slices(end + 1) = skipped_record(stride, decimated_length, ...
                raw_rows, raw_cols, ...
                'need at least 2 columns in the Hankel slice');
            continue;
        end

        info_i = init_slice_info();
        info_i.stride = stride;
        info_i.decimated_length = decimated_length;
        info_i.raw_rows = raw_rows;
        info_i.raw_cols = raw_cols;
        info_i.sample_indices = sample_indices;
        candidate_slices(end + 1) = info_i; %#ok<AGROW>
    end

    if numel(candidate_slices) < params.min_num_slices
        error('tensor_validate_params:TooFewSlices', ...
            ['Only %d valid slices remain after validation; need at least %d. ', ...
             'Adjust strides or base_window_length.'], ...
            numel(candidate_slices), params.min_num_slices);
    end

    max_common_rows = min([candidate_slices.raw_rows]);
    max_common_cols = min([candidate_slices.raw_cols]);
    block_num_rows = pick_block_size(params.block_num_rows, max_common_rows, 'rows', max_common_rows);
    block_num_cols = pick_block_size(params.block_num_cols, ...
        min(max_common_cols, max(block_num_rows, 2)), 'cols', ...
        min(max_common_cols, max(block_num_rows, 2)));
    if block_num_rows < params.min_common_rows || block_num_cols < params.min_common_cols
        error('tensor_validate_params:CommonShapeTooSmall', ...
            ['Requested block shape is %d x %d, below the configured minimum ', ...
             'of %d x %d.'], ...
            block_num_rows, block_num_cols, params.min_common_rows, params.min_common_cols);
    end

    block_row_step = pick_step(params.block_row_step, block_num_rows);
    block_col_step = pick_step(params.block_col_step, block_num_cols);
    filtered_slices = repmat(init_slice_info(), 0, 1);
    for i = 1:numel(candidate_slices)
        info_i = candidate_slices(i);
        info_i.block_num_rows = block_num_rows;
        info_i.block_num_cols = block_num_cols;
        info_i.block_row_step = block_row_step;
        info_i.block_col_step = block_col_step;
        info_i.row_starts = make_starts(info_i.raw_rows, block_num_rows, block_row_step);
        info_i.col_starts = make_starts(info_i.raw_cols, block_num_cols, block_col_step);
        info_i.num_row_blocks = numel(info_i.row_starts);
        info_i.num_col_blocks = numel(info_i.col_starts);
        info_i.num_blocks = info_i.num_row_blocks * info_i.num_col_blocks;
        if info_i.num_blocks < params.min_blocks_per_stride
            skipped_slices(end + 1) = skipped_record(info_i.stride, info_i.decimated_length, ...
                info_i.raw_rows, info_i.raw_cols, ...
                sprintf('only %d blocks available; need at least %d', ...
                info_i.num_blocks, params.min_blocks_per_stride));
            continue;
        end
        filtered_slices(end + 1) = info_i; %#ok<AGROW>
    end

    if numel(filtered_slices) < params.min_num_slices
        error('tensor_validate_params:TooFewBlockSlices', ...
            ['Only %d valid strides remain after block extraction; need at least %d. ', ...
             'Adjust block size, step, or strides.'], ...
            numel(filtered_slices), params.min_num_slices);
    end

    coverage_counts = zeros(N, 1);
    for i = 1:numel(filtered_slices)
        info_i = filtered_slices(i);
        support_len = block_num_rows + block_num_cols - 1;
        for rr = 1:numel(info_i.row_starts)
            row_start = info_i.row_starts(rr);
            for cc = 1:numel(info_i.col_starts)
                col_start = info_i.col_starts(cc);
                local_pos = (row_start + col_start - 1):(row_start + col_start + support_len - 2);
                coverage_counts(info_i.sample_indices(local_pos)) = ...
                    coverage_counts(info_i.sample_indices(local_pos)) + 1;
            end
        end
    end

    info = struct();
    info.signal_length = N;
    info.start_index = params.start_index;
    info.base_window_length = params.base_window_length;
    info.valid_strides = [filtered_slices.stride];
    info.num_valid_slices = numel(filtered_slices);
    info.common_rows = block_num_rows;
    info.common_cols = block_num_cols;
    info.block_num_rows = block_num_rows;
    info.block_num_cols = block_num_cols;
    info.block_row_step = block_row_step;
    info.block_col_step = block_col_step;
    info.tensor_size = [block_num_rows, block_num_cols, sum([filtered_slices.num_blocks])];
    info.candidate_slices = filtered_slices;
    info.skipped_slices = skipped_slices;
    info.coverage_counts = coverage_counts;
    info.uncovered_indices = find(coverage_counts == 0);
    info.coverage_fraction = mean(coverage_counts > 0);
end

function info = init_slice_info()
    info = struct( ...
        'stride', nan, ...
        'decimated_length', nan, ...
        'raw_rows', nan, ...
        'raw_cols', nan, ...
        'sample_indices', [], ...
        'block_num_rows', nan, ...
        'block_num_cols', nan, ...
        'block_row_step', nan, ...
        'block_col_step', nan, ...
        'row_starts', [], ...
        'col_starts', [], ...
        'num_row_blocks', nan, ...
        'num_col_blocks', nan, ...
        'num_blocks', nan);
end

function info = init_skipped_slice()
    info = struct( ...
        'stride', nan, ...
        'decimated_length', nan, ...
        'raw_rows', nan, ...
        'raw_cols', nan, ...
        'reason', "");
end

function info = skipped_record(stride, decimated_length, raw_rows, raw_cols, reason)
    info = init_skipped_slice();
    info.stride = stride;
    info.decimated_length = decimated_length;
    info.raw_rows = raw_rows;
    info.raw_cols = raw_cols;
    info.reason = string(reason);
end

function block_size = pick_block_size(value, max_value, label, default_value)
    if isempty(value)
        block_size = default_value;
    else
        block_size = value;
    end
    if ~isscalar(block_size) || block_size < 2 || block_size ~= floor(block_size)
        error('tensor_validate_params:BadBlockSize', ...
            'block_num_%s must be an integer >= 2.', label);
    end
    if block_size > max_value
        error('tensor_validate_params:BlockTooLarge', ...
            'block_num_%s=%d exceeds the common admissible maximum of %d.', ...
            label, block_size, max_value);
    end
end

function step = pick_step(value, block_size)
    if isempty(value)
        step = max(1, floor(block_size / 2));
    else
        step = value;
    end
    if ~isscalar(step) || step < 1 || step ~= floor(step)
        error('tensor_validate_params:BadBlockStep', ...
            'Block step parameters must be positive integers.');
    end
end

function starts = make_starts(raw_size, block_size, step)
    max_start = raw_size - block_size + 1;
    starts = 1:step:max_start;
    if isempty(starts)
        starts = 1;
        return;
    end
    if starts(end) ~= max_start
        starts = [starts, max_start];
    end
    starts = unique(starts);
end
