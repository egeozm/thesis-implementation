function params = ssd_merge_params(base_params, varargin)
%SSD_MERGE_PARAMS Merge frozen SSA grouping params with SSD-specific defaults.
%
%   params = ssd_merge_params(best_params)
%   params = ssd_merge_params(best_params, struct('max_modes', 8))
%
%   base_params must contain: r_max, sv_ratio_thresh, freq_tol, min_osc_freq
%   Optional overrides via additional struct(s): L_list, max_modes,
%   min_mode_energy_frac, min_modes_before_score_stop,
%   min_score_ratio, min_score_drop_ratio, fs, extract_trend_first, trend_L

    params = base_params;
    if ~isfield(params, 'fs') || isempty(params.fs)
        params.fs = 1;
    end
    if ~isfield(params, 'L_list') || isempty(params.L_list)
        params.L_list = [200, 400, 600, 1000];
    end
    if ~isfield(params, 'max_modes') || isempty(params.max_modes)
        params.max_modes = 4;
    end
    if ~isfield(params, 'min_mode_energy_frac') || isempty(params.min_mode_energy_frac)
        params.min_mode_energy_frac = 1e-2;
    end
    if ~isfield(params, 'min_modes_before_score_stop') || isempty(params.min_modes_before_score_stop)
        params.min_modes_before_score_stop = 2;
    end
    if ~isfield(params, 'min_score_ratio') || isempty(params.min_score_ratio)
        params.min_score_ratio = 0.1;
    end
    if ~isfield(params, 'min_score_drop_ratio') || isempty(params.min_score_drop_ratio)
        params.min_score_drop_ratio = 0.5;
    end
    if ~isfield(params, 'extract_trend_first') || isempty(params.extract_trend_first)
        params.extract_trend_first = false;
    end
    if ~isfield(params, 'trend_L')
        params.trend_L = [];
    end

    required = {'r_max', 'sv_ratio_thresh', 'freq_tol', 'min_osc_freq'};
    for k = 1:numel(required)
        if ~isfield(params, required{k})
            error('ssd_merge_params:MissingField', 'Missing field params.%s', required{k});
        end
    end

    for a = 1:numel(varargin)
        o = varargin{a};
        if isempty(o)
            continue;
        end
        f = fieldnames(o);
        for i = 1:numel(f)
            params.(f{i}) = o.(f{i});
        end
    end
end
