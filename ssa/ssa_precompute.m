function features = ssa_precompute(x, L, r_max, fs)
%SSA_PRECOMPUTE Compute SSA decomposition and leading elementary reconstructions once.
%
%   features = ssa_precompute(x, L, r_max, fs)

    if nargin < 4 || isempty(fs)
        fs = 1;
    end

    x = x(:);
    N = numel(x);
    H = embed_hankel(x, L);
    [U, S, V] = svd(H, 'econ');
    singular_values = diag(S);
    rank_H = numel(singular_values);
    r_use = min(r_max, rank_H);

    elementary = zeros(N, r_use);
    dominant_freqs = nan(1, r_use);
    peak_values = nan(1, r_use);

    for r = 1:r_use
        Hr = S(r, r) * (U(:, r) * V(:, r).');
        elementary(:, r) = diagonal_average(Hr, N);
        [dominant_freqs(r), ~, peak_values(r)] = estimate_dominant_frequency(elementary(:, r), fs);
    end

    features = struct();
    features.x = x;
    features.N = N;
    features.L = L;
    features.K = size(H, 2);
    features.rank = rank_H;
    features.U = U;
    features.S = S;
    features.V = V;
    features.singular_values = singular_values;
    features.r_max = r_use;
    features.elementary = elementary;
    features.dominant_freqs = dominant_freqs;
    features.peak_values = peak_values;
    features.fs = fs;
end
