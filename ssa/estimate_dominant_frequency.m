function [f_peak, bin_ix, peak_value] = estimate_dominant_frequency(x, fs)
%ESTIMATE_DOMINANT_FREQUENCY Dominant non-DC FFT peak in cycles per sample.
%
%   [f_peak, bin_ix, peak_value] = estimate_dominant_frequency(x, fs)

    if nargin < 2 || isempty(fs)
        fs = 1;
    end

    x = x(:);
    N = numel(x);
    if N < 8
        f_peak = nan;
        bin_ix = nan;
        peak_value = nan;
        return;
    end

    Xf = fft(x);
    n_half = floor(N / 2) + 1;
    P = abs(Xf(1:n_half));
    P(1) = 0; % ignore DC

    [peak_value, bin_ix] = max(P);
    if isempty(bin_ix) || peak_value <= eps
        f_peak = 0;
        bin_ix = 1;
        peak_value = 0;
        return;
    end

    f_peak = (bin_ix - 1) * (fs / N);
end
