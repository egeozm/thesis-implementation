function H = embed_hankel(x, L)
%EMBED_HANKEL Trajectory (Hankel) matrix for SSA.
%
%   H = embed_hankel(x, L)
%
%   x: column vector length N, L window length with 1 < L < N.
%   H(i,j) = x(i+j-1), size L x K, K = N - L + 1.

    x = x(:);
    N = numel(x);
    if ~(L > 1 && L < N)
        error('embed_hankel:BadL', 'Require 1 < L < N (N=%d, L=%d).', N, L);
    end
    K = N - L + 1;
    H = hankel(x(1:L), x(L:N));
    if size(H, 1) ~= L || size(H, 2) ~= K
        error('embed_hankel:Size', 'Unexpected Hankel size.');
    end
end
