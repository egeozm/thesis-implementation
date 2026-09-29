function [grouped_est, grouping] = group_cp_terms_oracle(est_rank1, true_components, opts)
%GROUP_CP_TERMS_ORACLE Group CP rank-1 terms using synthetic truth.
%
%   [grouped_est, grouping] = group_cp_terms_oracle(est_rank1, true_components, opts)
%
%   Assigns each reconstructed rank-1 term to the truth component with the
%   largest absolute Pearson correlation. Terms below opts.group_oracle_min_rho
%   are left out of the grouped estimates and tracked as residual members.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'group_oracle_min_rho') || isempty(opts.group_oracle_min_rho)
        opts.group_oracle_min_rho = 0.05;
    end

    est_rank1 = est_rank1(:,:);
    true_components = true_components(:,:);
    [N, rank_k] = size(est_rank1);
    [N_truth, num_truth] = size(true_components);

    if N_truth ~= N
        error('group_cp_terms_oracle:SizeMismatch', ...
            'est_rank1 and true_components must have the same number of rows.');
    end

    grouped_est = zeros(N, num_truth);
    corr_signed = zeros(num_truth, rank_k);
    corr_abs = zeros(num_truth, rank_k);
    assignments = zeros(rank_k, 1);
    sign_flips = ones(rank_k, 1);
    members = cell(num_truth, 1);

    for i = 1:num_truth
        for r = 1:rank_k
            c = pearson_raw(true_components(:, i), est_rank1(:, r));
            corr_signed(i, r) = c;
            corr_abs(i, r) = abs(c);
        end
    end

    [best_abs_rho, best_truth_idx] = max(corr_abs, [], 1);
    for r = 1:rank_k
        if best_abs_rho(r) < opts.group_oracle_min_rho
            continue;
        end

        i = best_truth_idx(r);
        assignments(r) = i;
        c = corr_signed(i, r);
        if c < 0
            sign_flips(r) = -1;
        end
        grouped_est(:, i) = grouped_est(:, i) + sign_flips(r) * est_rank1(:, r);
        members{i}(end + 1) = r; %#ok<AGROW>
    end

    residual_terms = find(assignments == 0);
    grouping = struct();
    grouping.mode = "oracle";
    grouping.params = opts;
    grouping.corr_matrix = corr_signed;
    grouping.abs_corr_matrix = corr_abs;
    grouping.best_abs_rho = best_abs_rho(:);
    grouping.best_truth_idx = best_truth_idx(:);
    grouping.assignments = assignments;
    grouping.sign_flips = sign_flips;
    grouping.component_members = members;
    grouping.residual_terms = residual_terms(:);
    grouping.num_grouped_components = num_truth;
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
