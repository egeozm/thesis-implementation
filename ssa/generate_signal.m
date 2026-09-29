function [x, components, t] = generate_signal(type, N, snr_db, seed)
%GENERATE_SIGNAL Synthetic signals from thesis Section 3.1.
%
%   [x, components, t] = generate_signal(type, N, snr_db, seed)
%
%   type:
%     'multiscale' or 1 — trend + slow sine + fast sine + noise
%     'close'      or 2 — two close frequencies + noise
%
%   snr_db: SNR in dB; use Inf for noiseless (sigma = 0).
%   seed:   optional; if nonempty, rng(seed) is applied.
%
%   components: struct with fields depending on type (column vectors N x 1).
%   t: time indices 1:N as column vector.

    if nargin < 4
        seed = [];
    end
    if ~isempty(seed)
        rng(seed);
    end

    t = (1:N).';
    n = t;

    if ischar(type) || isstring(type)
        ts = lower(char(type));
        if strcmp(ts, 'multiscale')
            type = 1;
        elseif strcmp(ts, 'close')
            type = 2;
        else
            error('generate_signal:UnknownType', 'Unknown type "%s".', type);
        end
    end

    switch type
        case 1
            [x, components] = multiscale_signal(n, N, snr_db);
        case 2
            [x, components] = close_frequency_signal(n, N, snr_db);
        otherwise
            error('generate_signal:UnknownType', 'type must be 1, 2, ''multiscale'', or ''close''.');
    end
end

function [x, c] = multiscale_signal(n, N, snr_db)
    Tn = 0.002 * (n / N).^2 - 0.15 * (n / N);
    As = 1;
    fs = 0.01;
    Af = 0.5;
    ff = 0.08;
    phi = 2 * pi * rand();
    slow = As * sin(2 * pi * fs * n);
    fast = Af * sin(2 * pi * ff * n + phi);

    x_clean = Tn + slow + fast;
    sigma = noise_sigma(x_clean, snr_db);
    noise = sigma * randn(N, 1);
    x = x_clean + noise;

    c = struct('trend', Tn, 'slow', slow, 'fast', fast, 'noise', noise, ...
        'x_clean', x_clean);
end

function [x, c] = close_frequency_signal(n, N, snr_db)
    f1 = 0.05;
    f2 = 0.055;
    phi = 2 * pi * rand();
    comp1 = sin(2 * pi * f1 * n);
    comp2 = sin(2 * pi * f2 * n + phi);
    x_clean = comp1 + comp2;
    sigma = noise_sigma(x_clean, snr_db);
    noise = sigma * randn(N, 1);
    x = x_clean + noise;

    c = struct('comp1', comp1, 'comp2', comp2, 'noise', noise, ...
        'x_clean', x_clean);
end

function sigma = noise_sigma(x_clean, snr_db)
    if isinf(snr_db) || isnan(snr_db)
        sigma = 0;
        return
    end
    v = var(x_clean, 0, 'all');
    if v == 0
        sigma = 0;
    else
        sigma = sqrt(v / (10^(snr_db / 10)));
    end
end
