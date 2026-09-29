function [pairs, scores] = ssd_list_oscillatory_pairs(features, params)
%SSD_LIST_OSCILLATORY_PAIRS Rank adjacent oscillatory eigentriple pairs (SSD helper).
%
%   Uses the same feasibility rules as auto_group_ssa for a single residual SSA view.
%   pairs: cell array of 1x2 index rows, scores: descending.
%
%   For an adjacent pair (r, r+1), define:
%     sigma_mass   = sigma_r + sigma_{r+1}
%     sv_ratio     = min(sigma_r, sigma_{r+1}) / max(sigma_r, sigma_{r+1})
%     freq_diff    = |f_r - f_{r+1}|
%     freq_penalty = 1 + freq_diff / max(freq_tol, eps)
%
%   The pair score is the multiplicative heuristic
%
%     score(r, r+1) = sigma_mass * sv_ratio / freq_penalty
%
%   i.e. no separate weights are used. Larger singular-value mass, better balance
%   between the two singular values, and tighter dominant-frequency agreement all
%   increase the score.

    r_use = min(params.r_max, features.r_max);
    singular_values = features.singular_values(1:r_use);
    dominant_freqs = features.dominant_freqs(1:r_use);
    peak_values = features.peak_values(1:r_use);

    candidate_pairs = {};
    candidate_scores = [];

    for r = 1:(r_use - 1)
        sv_ratio = min(singular_values(r), singular_values(r + 1)) / ...
            max(singular_values(r), singular_values(r + 1));
        freq_diff = abs(dominant_freqs(r) - dominant_freqs(r + 1));
        mean_freq = mean(dominant_freqs([r, r + 1]));
        has_energy = max(peak_values([r, r + 1])) > eps;

        if has_energy && sv_ratio >= params.sv_ratio_thresh && ...
                freq_diff <= params.freq_tol && mean_freq >= params.min_osc_freq
            candidate_pairs{end + 1} = [r, r + 1]; %#ok<AGROW>
            candidate_scores(end + 1) = pair_score( ...
                singular_values([r, r + 1]), sv_ratio, freq_diff, params.freq_tol); %#ok<AGROW>
        end
    end

    if isempty(candidate_pairs)
        pairs = {};
        scores = [];
        return;
    end

    [scores, order] = sort(candidate_scores, 'descend');
    pairs = candidate_pairs(order);
end

function score = pair_score(sigmas, sv_ratio, freq_diff, freq_tol)
    sigma_mass = sum(sigmas);
    freq_penalty = 1 + (freq_diff / max(freq_tol, eps));
    score = sigma_mass * sv_ratio / freq_penalty;
end
