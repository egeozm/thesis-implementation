function grouping = auto_group_ssa(features, params)
%AUTO_GROUP_SSA Fixed automatic grouping rule for leading SSA eigentriples.
%
%   grouping = auto_group_ssa(features, params)
%
%   Required params fields:
%     r_max
%     sv_ratio_thresh
%     freq_tol
%     min_osc_freq

    params = validate_params(params);
    r_use = min(params.r_max, features.r_max);

    singular_values = features.singular_values(1:r_use);
    dominant_freqs = features.dominant_freqs(1:r_use);
    peak_values = features.peak_values(1:r_use);

    candidate_pairs = {};
    candidate_frequencies = [];
    candidate_ratios = [];
    candidate_freq_diffs = [];
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
            candidate_frequencies(end + 1) = mean_freq; %#ok<AGROW>
            candidate_ratios(end + 1) = sv_ratio; %#ok<AGROW>
            candidate_freq_diffs(end + 1) = freq_diff; %#ok<AGROW>
            candidate_scores(end + 1) = pair_score( ...
                singular_values([r, r + 1]), sv_ratio, freq_diff, params.freq_tol); %#ok<AGROW>
        end
    end

    [pairs, pair_frequencies, pair_ratios, pair_freq_diffs, pair_scores] = ...
        select_top_pairs(candidate_pairs, candidate_frequencies, candidate_ratios, ...
        candidate_freq_diffs, candidate_scores, 2);

    used = false(1, r_use);
    selected_idx = flatten_pairs(pairs);
    used(selected_idx) = true;

    unpaired_head_idx = find(~used);
    low_freq_unpaired_idx = unpaired_head_idx(dominant_freqs(unpaired_head_idx) < params.min_osc_freq);

    if isempty(low_freq_unpaired_idx) && ~isempty(unpaired_head_idx)
        if isempty(selected_idx)
            low_freq_unpaired_idx = unpaired_head_idx(1);
        else
            head_prefix = unpaired_head_idx(unpaired_head_idx < min(selected_idx));
            if isempty(head_prefix)
                low_freq_unpaired_idx = unpaired_head_idx(1);
            else
                low_freq_unpaired_idx = head_prefix;
            end
        end
    end

    trend_idx = low_freq_unpaired_idx;
    grouped_signal_idx = sort(unique([trend_idx, selected_idx]));
    residual_idx = setdiff(1:features.rank, grouped_signal_idx);

    grouping = struct();
    grouping.params = params;
    grouping.r_max_used = r_use;
    grouping.unpaired_head_idx = unpaired_head_idx;
    grouping.low_freq_unpaired_idx = low_freq_unpaired_idx;
    grouping.trend_idx = trend_idx;
    grouping.oscillatory_pairs = pairs;
    grouping.pair_frequencies = pair_frequencies;
    grouping.pair_ratios = pair_ratios;
    grouping.pair_freq_diffs = pair_freq_diffs;
    grouping.pair_scores = pair_scores;
    grouping.grouped_signal_idx = grouped_signal_idx;
    grouping.residual_idx = residual_idx;
end

function params = validate_params(params)
    required_fields = {'r_max', 'sv_ratio_thresh', 'freq_tol', 'min_osc_freq'};
    for k = 1:numel(required_fields)
        if ~isfield(params, required_fields{k})
            error('auto_group_ssa:MissingParam', 'Missing params.%s', required_fields{k});
        end
    end
end

function score = pair_score(sigmas, sv_ratio, freq_diff, freq_tol)
    sigma_mass = sum(sigmas);
    freq_penalty = 1 + (freq_diff / max(freq_tol, eps));
    score = sigma_mass * sv_ratio / freq_penalty;
end

function [pairs, pair_frequencies, pair_ratios, pair_freq_diffs, pair_scores] = ...
        select_top_pairs(candidate_pairs, candidate_frequencies, candidate_ratios, ...
        candidate_freq_diffs, candidate_scores, max_pairs)
    pairs = {};
    pair_frequencies = [];
    pair_ratios = [];
    pair_freq_diffs = [];
    pair_scores = [];

    if isempty(candidate_pairs)
        return;
    end

    [~, order] = sort(candidate_scores, 'descend');
    used = false(1, max(flatten_pairs(candidate_pairs)));

    for idx = order
        pair = candidate_pairs{idx};
        if any(used(pair))
            continue;
        end
        pairs{end + 1} = pair; %#ok<AGROW>
        pair_frequencies(end + 1) = candidate_frequencies(idx); %#ok<AGROW>
        pair_ratios(end + 1) = candidate_ratios(idx); %#ok<AGROW>
        pair_freq_diffs(end + 1) = candidate_freq_diffs(idx); %#ok<AGROW>
        pair_scores(end + 1) = candidate_scores(idx); %#ok<AGROW>
        used(pair) = true;

        if numel(pairs) >= max_pairs
            break;
        end
    end

    if ~isempty(pair_frequencies)
        [pair_frequencies, freq_order] = sort(pair_frequencies);
        pairs = pairs(freq_order);
        pair_ratios = pair_ratios(freq_order);
        pair_freq_diffs = pair_freq_diffs(freq_order);
        pair_scores = pair_scores(freq_order);
    end
end

function idx = flatten_pairs(pairs)
    idx = [];
    for k = 1:numel(pairs)
        idx = [idx, pairs{k}]; %#ok<AGROW>
    end
end
