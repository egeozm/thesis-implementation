function [assignment, corr_matrix] = match_components(true_components, est_components)
%MATCH_COMPONENTS Greedy one-to-one matching by maximum absolute Pearson correlation.
%
%   [assignment, corr_matrix] = match_components(true_components, est_components)
%
%   true_components: N x Nt matrix (columns = reference components).
%   est_components:  N x Ne matrix (columns = estimated components).
%
%   assignment: Nt x 1 — assignment(i) is the column index in est_components matched
%               to true column i. If Nt > Ne, extra entries are NaN (unmatched).
%   corr_matrix: Nt x Ne matrix of absolute Pearson correlations.

    Nt = size(true_components, 2);
    Ne = size(est_components, 2);
    corr_matrix = zeros(Nt, Ne);

    for i = 1:Nt
        for j = 1:Ne
            c = pearson_corr(true_components(:, i), est_components(:, j));
            corr_matrix(i, j) = abs(c);
        end
    end

    assignment = nan(Nt, 1);
    C = corr_matrix;
    for k = 1:min(Nt, Ne)
        [val, linear_ix] = max(C(:));
        if isnan(val) || val == -Inf
            break;
        end
        [i_row, j_col] = ind2sub(size(C), linear_ix);
        assignment(i_row) = j_col;
        C(i_row, :) = -Inf;
        C(:, j_col) = -Inf;
    end
end

function c = pearson_corr(a, b)
    a = a(:) - mean(a);
    b = b(:) - mean(b);
    den = norm(a) * norm(b);
    if den == 0
        c = 0;
    else
        c = (a.' * b) / den;
    end
end
