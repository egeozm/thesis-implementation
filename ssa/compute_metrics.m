function metrics = compute_metrics(true_components, est_components, assignment, fs)
%COMPUTE_METRICS NMSE, Pearson correlation, dominant peak frequency error, leakage.
%
%   metrics = compute_metrics(true_components, est_components, assignment, fs)
%
%   true_components, est_components: N x Nt / N x Ne matrices.
%   assignment: Nt x 1 from match_components (NaN if unmatched).
%   fs: sampling rate (Hz); use 1 for normalized time (cycles per sample for FFT bins).
%
%   metrics fields (vectors length Nt, NaN where unmatched or undefined):
%     nmse, pearson, peak_freq_error, leakage
%   Also: true_idx, est_idx (same as assignment), dominant_freq_true, dominant_freq_est

    if nargin < 4 || isempty(fs)
        fs = 1;
    end

    [N, Nt] = size(true_components);
    Ne = size(est_components, 2);
    if size(est_components, 1) ~= N
        error('compute_metrics:SizeMismatch', ...
            'true_components and est_components must have the same number of rows.');
    end

    metrics.nmse = nan(Nt, 1);
    metrics.pearson = nan(Nt, 1);
    metrics.peak_freq_error = nan(Nt, 1);
    metrics.leakage = nan(Ne, 1);
    metrics.dominant_freq_true = nan(Nt, 1);
    metrics.dominant_freq_est = nan(Nt, 1);
    metrics.true_idx = (1:Nt).';
    metrics.est_idx = assignment(:);

    for i = 1:Nt
        s = true_components(:, i);
        j = assignment(i);
        if isnan(j) || j < 1 || j > Ne
            continue;
        end
        shat = est_components(:, j);
        s_energy = sum(s.^2);

        % A zero-energy reference component (e.g. noiseless "noise") does not admit
        % meaningful NMSE/correlation metrics, so leave those entries as NaN.
        if s_energy <= eps
            continue;
        end

        metrics.nmse(i) = sum((s - shat).^2) / s_energy;
        metrics.pearson(i) = pearson_raw(s, shat);

        [ft, ~] = dominant_peak_freq(s, fs);
        [fe, ~] = dominant_peak_freq(shat, fs);
        metrics.dominant_freq_true(i) = ft;
        metrics.dominant_freq_est(i) = fe;
        if ~isnan(ft) && ~isnan(fe)
            metrics.peak_freq_error(i) = abs(ft - fe);
        end
    end

    % Leakage per estimated column: correlation energy on non-matched true components
    for j = 1:Ne
        matched_true = find(assignment == j, 1);
        ej = est_components(:, j);
        leak = 0;
        for k = 1:Nt
            if ~isempty(matched_true) && k == matched_true
                continue;
            end
            c = pearson_raw(ej, true_components(:, k));
            leak = leak + c^2;
        end
        metrics.leakage(j) = leak;
    end
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

function [f_peak, bin_ix] = dominant_peak_freq(x, fs)
    x = x(:);
    N = numel(x);
    if N < 8
        f_peak = nan;
        bin_ix = nan;
        return;
    end
    Xf = fft(x);
    n_half = floor(N / 2) + 1;
    P = abs(Xf(1:n_half));
    P(1) = 0; % ignore DC
    if numel(P) < 2
        f_peak = nan;
        bin_ix = nan;
        return;
    end
    [~, bin_ix] = max(P);
    % Bin bin_ix corresponds to frequency (bin_ix-1)*fs/N in Hz
    f_peak = (bin_ix - 1) * (fs / N);
end
