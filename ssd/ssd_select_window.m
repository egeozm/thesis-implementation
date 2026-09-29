function [L_best, feat_best, pair, best_score] = ssd_select_window(residual, params)
%SSD_SELECT_WINDOW Choose L and best single oscillatory pair for current residual.
%
%   [L_best, feat_best, pair, best_score] = ssd_select_window(residual, params)
%
%   params: ssd_merge_params output (L_list, r_max, fs, ...).
%   If no pair is feasible at any L, L_best and feat_best are empty, pair is [].
%
%   Selection rule:
%   1. For each candidate L, rank all feasible adjacent pairs by score.
%   2. Keep only the top-scoring pair for that L.
%   3. Choose the single best pair across all L by comparing those per-L maxima.
%
%   This is equivalent to choosing the global maximum over all feasible adjacent
%   pairs across all L, since any non-maximal pair within a fixed L cannot beat the
%   maximal pair from the same L.

    residual = residual(:);
    N = numel(residual);
    L_list = params.L_list(:).';

    best_score = -inf;
    L_best = [];
    feat_best = [];
    pair = [];

    for L = L_list
        if ~(L >= 2 && L < N)
            continue;
        end
        features = ssa_precompute(residual, L, params.r_max, params.fs);
        [pairs, scores] = ssd_list_oscillatory_pairs(features, params);
        if isempty(pairs)
            continue;
        end
        sc = scores(1);
        pr = pairs{1};
        if sc > best_score || (sc == best_score && better_tie(L, pr, L_best, pair))
            best_score = sc;
            L_best = L;
            feat_best = features;
            pair = pr;
        end
    end

    if isempty(L_best)
        best_score = nan;
    end
end

function yes = better_tie(L_new, pair_new, L_old, pair_old)
    if isempty(L_old)
        yes = true;
        return;
    end
    if L_new > L_old
        yes = true;
    elseif L_new < L_old
        yes = false;
    else
        yes = lex_less(pair_new, pair_old);
    end
end

function yes = lex_less(a, b)
    yes = a(1) < b(1) || (a(1) == b(1) && a(2) < b(2));
end
