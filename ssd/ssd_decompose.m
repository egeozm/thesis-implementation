function out = ssd_decompose(x, params)
%SSD_DECOMPOSE Singular Spectrum Decomposition (iterative SSA on residual).
%
%   out = ssd_decompose(x, params)
%
%   params: output of ssd_merge_params (includes L_list, max_modes, r_max, fs, ...).
%
%   out fields:
%     modes          N x M extracted oscillatory modes (columns)
%     trend          N x 1 optional trend extracted before oscillatory modes
%     residual       N x 1 final residual
%     L_history      1 x M window lengths used per extraction
%     pair_history   1 x M cell of pair index rows
%     scores         1 x M selection scores
%     num_iterations scalar M

    x = x(:);
    params = ssd_merge_params(params);

    modes = [];
    L_hist = [];
    pair_hist = {};
    score_hist = [];
    first_score = nan;
    prev_score = nan;

    r = x;
    trend = zeros(size(x));
    trend_L = nan;
    trend_idx = [];
    trend_extracted = false;

    if params.extract_trend_first
        [trend, trend_L, trend_idx, trend_extracted] = ssd_extract_initial_trend(x, params);
        if trend_extracted
            r = r - trend;
        end
    end

    iter = 0;

    while iter < params.max_modes
        nr_before = norm(r);
        if nr_before <= eps
            break;
        end

        [L_best, feat_best, pair, sc] = ssd_select_window(r, params);
        if isempty(pair) || isnan(sc)
            break;
        end

        y = ssd_extract_mode_at_features(feat_best, pair);
        next_iter = iter + 1;
        if ssd_stop_rule(y, nr_before, next_iter, sc, first_score, prev_score, params)
            break;
        end

        iter = next_iter;
        modes = [modes, y]; %#ok<AGROW>
        L_hist(end + 1) = L_best; %#ok<AGROW>
        pair_hist{end + 1} = pair; %#ok<AGROW>
        score_hist(end + 1) = sc; %#ok<AGROW>
        if isnan(first_score)
            first_score = sc;
        end
        prev_score = sc;

        r = r - y;
    end

    out = struct();
    out.trend = trend;
    out.trend_L = trend_L;
    out.trend_idx = trend_idx;
    out.trend_extracted = trend_extracted;
    out.modes = modes;
    out.residual = r;
    out.L_history = L_hist;
    out.pair_history = pair_hist;
    out.scores = score_hist;
    out.num_iterations = size(modes, 2);
end

function [trend, trend_L, trend_idx, extracted] = ssd_extract_initial_trend(x, params)
    x = x(:);
    trend = zeros(size(x));
    trend_L = nan;
    trend_idx = [];
    extracted = false;

    L_list = params.L_list(:).';
    if isfield(params, 'trend_L') && ~isempty(params.trend_L)
        L_list = params.trend_L(:).';
    end
    L_list = sort(unique(L_list), 'descend');

    for L = L_list
        if ~(L >= 2 && L < numel(x))
            continue;
        end

        features = ssa_precompute(x, L, params.r_max, params.fs);
        grouping = auto_group_ssa(features, params);
        idx = grouping.trend_idx(:).';
        if isempty(idx)
            idx = leading_low_frequency_idx(features, params);
        end
        if isempty(idx)
            continue;
        end

        trend = ssa_reconstruct_groups(features, {idx});
        trend = trend(:);
        trend_L = L;
        trend_idx = idx;
        extracted = norm(trend) > eps;
        if extracted
            return;
        end
    end
end

function idx = leading_low_frequency_idx(features, params)
    r_use = min(params.r_max, features.r_max);
    idx = [];
    for r = 1:r_use
        if features.dominant_freqs(r) < params.min_osc_freq
            idx(end + 1) = r; %#ok<AGROW>
        elseif ~isempty(idx)
            break;
        end
    end
end
