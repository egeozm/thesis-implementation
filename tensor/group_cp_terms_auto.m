function [grouped_est, grouping] = group_cp_terms_auto(est_rank1, lambda, params)
%GROUP_CP_TERMS_AUTO Group CP rank-1 terms by frequency and correlation.
%
%   [grouped_est, grouping] = group_cp_terms_auto(est_rank1, lambda, params)
%
%   Clusters reconstructed rank-1 terms into higher-level components using:
%   - dominant-frequency proximity,
%   - pairwise waveform correlation,
%   - total CP weight for ordering the final groups.

    if nargin < 3
        params = struct();
    end
    params = tensor_merge_params(params);

    est_rank1 = est_rank1(:,:);
    lambda = lambda(:);
    [N, rank_k] = size(est_rank1);
    if numel(lambda) ~= rank_k
        error('group_cp_terms_auto:SizeMismatch', ...
            'lambda must have one entry per rank-1 term.');
    end

    dominant_freqs = nan(rank_k, 1);
    pair_corr = zeros(rank_k, rank_k);
    abs_pair_corr = zeros(rank_k, rank_k);
    for r = 1:rank_k
        dominant_freqs(r) = estimate_dominant_frequency(est_rank1(:, r), params.fs);
        for c = r:rank_k
            rho = pearson_raw(est_rank1(:, r), est_rank1(:, c));
            pair_corr(r, c) = rho;
            pair_corr(c, r) = rho;
            abs_pair_corr(r, c) = abs(rho);
            abs_pair_corr(c, r) = abs(rho);
        end
    end

    [~, order] = sort(abs(lambda), 'descend');
    clusters = {};
    cluster_info = empty_cluster_info();

    for order_idx = 1:numel(order)
        term_idx = order(order_idx);
        freq = dominant_freqs(term_idx);
        candidate_cluster = choose_cluster(term_idx, freq, clusters, cluster_info, dominant_freqs, ...
            abs_pair_corr, params.group_freq_tol, params.group_rho_merge_thresh, ...
            params.group_min_osc_freq);

        if isnan(candidate_cluster)
            aligned_wave = est_rank1(:, term_idx);
            clusters{end + 1} = aligned_wave; %#ok<AGROW>
            cluster_info(end + 1) = init_cluster( ...
                term_idx, aligned_wave, freq, lambda(term_idx), params.group_min_osc_freq); %#ok<AGROW>
            continue;
        end

        [aligned_wave, sign_flip, rho_to_cluster] = align_to_cluster(est_rank1(:, term_idx), clusters{candidate_cluster});
        clusters{candidate_cluster} = clusters{candidate_cluster} + aligned_wave;
        cluster_info(candidate_cluster).members(end + 1) = term_idx;
        cluster_info(candidate_cluster).member_signs(end + 1) = sign_flip;
        cluster_info(candidate_cluster).member_corr_to_cluster(end + 1) = rho_to_cluster;
        cluster_info(candidate_cluster).total_abs_lambda = ...
            cluster_info(candidate_cluster).total_abs_lambda + abs(lambda(term_idx));
        cluster_info(candidate_cluster).mean_frequency = mean(dominant_freqs(cluster_info(candidate_cluster).members), 'omitnan');
        cluster_info(candidate_cluster).is_trend = cluster_info(candidate_cluster).is_trend && ...
            is_low_frequency(freq, params.group_min_osc_freq);
    end

    if isempty(clusters)
        grouped_est = zeros(N, 0);
        grouping = empty_grouping(params, dominant_freqs, pair_corr, abs_pair_corr, lambda);
        return;
    end

    total_abs_lambda = arrayfun(@(s) s.total_abs_lambda, cluster_info);
    [~, cluster_order] = sort(total_abs_lambda, 'descend');
    clusters = clusters(cluster_order);
    cluster_info = cluster_info(cluster_order);

    grouped_est = zeros(N, numel(clusters));
    for k = 1:numel(clusters)
        grouped_est(:, k) = clusters{k};
        cluster_info(k).label = "group" + string(k);
    end

    grouping = struct();
    grouping.mode = "auto";
    grouping.params = params;
    grouping.order = order;
    grouping.dominant_freqs = dominant_freqs;
    grouping.lambda = lambda;
    grouping.pair_corr = pair_corr;
    grouping.abs_pair_corr = abs_pair_corr;
    grouping.clusters = cluster_info;
    grouping.num_grouped_components = size(grouped_est, 2);
end

function idx = choose_cluster(term_idx, freq, clusters, cluster_info, dominant_freqs, abs_pair_corr, ...
        freq_tol, rho_thresh, min_osc_freq)
    idx = nan;
    best_score = -inf;

    for k = 1:numel(clusters)
        members = cluster_info(k).members;
        member_freqs = dominant_freqs(members);
        freq_diff = abs(freq - mean(member_freqs, 'omitnan'));
        corr_score = max(abs_pair_corr(term_idx, members));

        if is_low_frequency(freq, min_osc_freq) && cluster_info(k).is_trend
            score = corr_score;
        elseif freq_diff <= freq_tol && corr_score >= rho_thresh
            score = corr_score - freq_diff / max(freq_tol, eps);
        else
            continue;
        end

        if score > best_score
            best_score = score;
            idx = k;
        end
    end
end

function cluster = init_cluster(term_idx, aligned_wave, freq, lambda_val, min_osc_freq)
    cluster = struct();
    cluster.members = term_idx;
    cluster.member_signs = 1;
    cluster.member_corr_to_cluster = 1;
    cluster.total_abs_lambda = abs(lambda_val);
    cluster.mean_frequency = freq;
    cluster.is_trend = is_low_frequency(freq, min_osc_freq);
    cluster.reference_energy = norm(aligned_wave);
    cluster.label = "";
end

function tf = is_low_frequency(freq, min_osc_freq)
    tf = ~isnan(freq) && freq < min_osc_freq;
end

function [aligned_wave, sign_flip, rho] = align_to_cluster(wave, cluster_sum)
    rho = pearson_raw(wave, cluster_sum);
    sign_flip = 1;
    if rho < 0
        sign_flip = -1;
    end
    aligned_wave = sign_flip * wave;
end

function grouping = empty_grouping(params, dominant_freqs, pair_corr, abs_pair_corr, lambda)
    grouping = struct();
    grouping.mode = "auto";
    grouping.params = params;
    grouping.order = [];
    grouping.dominant_freqs = dominant_freqs;
    grouping.lambda = lambda;
    grouping.pair_corr = pair_corr;
    grouping.abs_pair_corr = abs_pair_corr;
    grouping.clusters = empty_cluster_info();
    grouping.num_grouped_components = 0;
end

function s = empty_cluster_info()
    s = struct( ...
        'members', {}, ...
        'member_signs', {}, ...
        'member_corr_to_cluster', {}, ...
        'total_abs_lambda', {}, ...
        'mean_frequency', {}, ...
        'is_trend', {}, ...
        'reference_energy', {}, ...
        'label', {});
end

function c = pearson_raw(a, b)
    a = a(:) - mean(a);
    b = b(:) - mean(b);
    den = norm(a) * norm(b);
    if den == 0
        c = 0;
    else
        c = (a.' * b) / den;
    end
end
