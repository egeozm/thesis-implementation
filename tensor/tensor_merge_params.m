function params = tensor_merge_params(base_params, varargin)
%TENSOR_MERGE_PARAMS Merge tensor defaults with user overrides.
%
%   params = tensor_merge_params()
%   params = tensor_merge_params(struct('strides', [1 2 4], 'cpd_rank', 3))

    if nargin < 1 || isempty(base_params)
        params = struct();
    else
        params = base_params;
    end

    defaults = struct();
    defaults.fs = 1;
    defaults.strides = [1, 2, 4, 8];
    defaults.base_window_length = 400;
    defaults.start_index = 1;
    defaults.min_common_rows = 8;
    defaults.min_common_cols = 8;
    defaults.min_num_slices = 2;
    defaults.block_num_rows = [];
    defaults.block_num_cols = [];
    defaults.block_row_step = [];
    defaults.block_col_step = [];
    defaults.min_blocks_per_stride = 2;
    defaults.recon_weighting = 'uniform';
    defaults.reconstruction_mode = 'cpd';
    defaults.cpd_rank = 3;
    defaults.cpd_method = 'cpd';
    defaults.cpd_options = struct();
    defaults.cpd_initialization = [];
    defaults.grouping_mode = 'benchmark';
    defaults.group_r_max = 12;
    defaults.group_sv_ratio_thresh = 0.6;
    defaults.group_freq_tol = 0.02;
    defaults.group_rho_merge_thresh = 0.3;
    defaults.group_min_osc_freq = 0.005;
    defaults.group_oracle_min_rho = 0.05;
    defaults.benchmark_freq_tol = 0.01;
    defaults.benchmark_trend_freq_max = 0.005;
    defaults.benchmark_assign_min_score = 0.05;

    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(params, name) || isempty(params.(name))
            params.(name) = defaults.(name);
        end
    end

    for a = 1:numel(varargin)
        override = varargin{a};
        if isempty(override)
            continue;
        end
        names = fieldnames(override);
        for i = 1:numel(names)
            params.(names{i}) = override.(names{i});
        end
    end
end
