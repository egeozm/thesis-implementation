function y = diagonal_average(X, N)
%DIAGONAL_AVERAGING SSA Hankelization: average anti-diagonals i+j-1 = k.
%
%   y = diagonal_average(X, N)
%
%   X: L x K matrix; N must equal L + K - 1 (length of recovered series).

    [L, K] = size(X);
    if L + K - 1 ~= N
        error('diagonal_average:BadN', ...
            'Need N = L + K - 1 (got L=%d, K=%d, N=%d).', L, K, N);
    end

    y = zeros(N, 1);
    for k = 1:N
        s = 0;
        c = 0;
        for i = 1:L
            j = k - i + 1;
            if j >= 1 && j <= K
                s = s + X(i, j);
                c = c + 1;
            end
        end
        y(k) = s / c;
    end
end
