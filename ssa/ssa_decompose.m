function [reconstructed, singular_values, U, S, V] = ssa_decompose(x, L, groups)
%SSA_DECOMPOSE Classical SSA: Hankel embed, SVD, group elementary matrices, diagonal average.
%
%   [reconstructed, singular_values, U, S, V] = ssa_decompose(x, L, groups)
%
%   x: N x 1 column vector.
%   L: window length (1 < L < N).
%   groups: cell array; groups{g} is a vector of singular triplet indices (1-based),
%           each index in 1:min(L, K). Each group yields one reconstructed series.
%
%   reconstructed: N x numel(groups). Column g is the Hankel-averaged series for
%                  the sum of elementary matrices in groups{g}.

    x = x(:);
    N = numel(x);
    H = embed_hankel(x, L);
    [U, S, V] = svd(H, 'econ');
    singular_values = diag(S);
    R = numel(singular_values);

    num_groups = numel(groups);
    reconstructed = zeros(N, num_groups);

    for g = 1:num_groups
        idx = groups{g}(:).';
        idx = idx(idx >= 1 & idx <= R);
        if isempty(idx)
            continue;
        end
        HG = zeros(size(H));
        for r = idx
            HG = HG + S(r, r) * (U(:, r) * V(:, r).');
        end
        reconstructed(:, g) = diagonal_average(HG, N);
    end
end
