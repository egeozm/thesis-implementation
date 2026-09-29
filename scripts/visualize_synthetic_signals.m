%VISUALIZE_SYNTHETIC_SIGNALS Plot the thesis synthetic signals and components.
%
%   From project root:
%       run('scripts/visualize_synthetic_signals.m')
%
%   This script visualizes both synthetic signal families defined in
%   `ssa/generate_signal.m`:
%     1. multiscale
%     2. close
%
%   For each case it shows:
%     - noisy observation vs clean signal
%     - each ground-truth component
%     - single-sided FFT magnitude
%
%   Figures are also exported to the `results/` directory.

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));

resultsDir = fullfile(projRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Parameters
N = 2000;
snr_db = 20;      % use Inf for noiseless
seed = 42;
signalTypes = {'multiscale', 'close'};

for i = 1:numel(signalTypes)
    signalType = signalTypes{i};
    [x, components, t] = generate_signal(signalType, N, snr_db, seed);
    make_signal_figure(signalType, x, components, t, snr_db, resultsDir);
end

fprintf('\nSynthetic signal visualization finished.\n');

%% Local helpers
function make_signal_figure(signalType, x, components, t, snr_db, resultsDir)
    [componentMatrix, componentNames] = unpack_components(signalType, components);
    xClean = components.x_clean(:);
    x = x(:);
    t = t(:);

    fig = figure( ...
        'Name', sprintf('Synthetic signal: %s', signalType), ...
        'NumberTitle', 'off', ...
        'Color', 'w');
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    plot(ax1, t, x, 'Color', [0.00, 0.45, 0.74], 'LineWidth', 1.0, ...
        'DisplayName', 'observed');
    hold(ax1, 'on');
    plot(ax1, t, xClean, 'k--', 'LineWidth', 1.1, 'DisplayName', 'clean');
    title(ax1, sprintf('%s signal (SNR = %s dB)', pretty_name(signalType), snr_label(snr_db)));
    xlabel(ax1, 'Sample index');
    ylabel(ax1, 'Amplitude');
    legend(ax1, 'Location', 'best');
    grid(ax1, 'on');

    ax2 = nexttile;
    hold(ax2, 'on');
    colors = lines(size(componentMatrix, 2));
    for k = 1:size(componentMatrix, 2)
        plot(ax2, t, componentMatrix(:, k), 'LineWidth', 1.1, ...
            'Color', colors(k, :), 'DisplayName', componentNames{k});
    end
    title(ax2, 'Ground-truth components');
    xlabel(ax2, 'Sample index');
    ylabel(ax2, 'Amplitude');
    legend(ax2, 'Location', 'best');
    grid(ax2, 'on');

    ax3 = nexttile;
    [f, magObserved] = single_sided_fft(x);
    [~, magClean] = single_sided_fft(xClean);
    plot(ax3, f, magObserved, 'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.0, ...
        'DisplayName', 'observed');
    hold(ax3, 'on');
    plot(ax3, f, magClean, 'k--', 'LineWidth', 1.1, 'DisplayName', 'clean');
    title(ax3, 'Single-sided FFT magnitude');
    xlabel(ax3, 'Normalized frequency (cycles/sample)');
    ylabel(ax3, '|X(f)|');
    xlim(ax3, [0, 0.15]);
    legend(ax3, 'Location', 'best');
    grid(ax3, 'on');

    exportName = sprintf('synthetic_signal_%s.png', signalType);
    exportgraphics(fig, fullfile(resultsDir, exportName), 'Resolution', 300);
end

function [componentMatrix, componentNames] = unpack_components(signalType, components)
    switch lower(signalType)
        case 'multiscale'
            componentMatrix = [components.trend, components.slow, components.fast, components.noise];
            componentNames = {'trend', 'slow', 'fast', 'noise'};
        case 'close'
            componentMatrix = [components.comp1, components.comp2, components.noise];
            componentNames = {'comp1', 'comp2', 'noise'};
        otherwise
            error('Unknown signal type "%s".', signalType);
    end
end

function [f, magnitude] = single_sided_fft(x)
    x = x(:);
    N = numel(x);
    Y = fft(x);
    P2 = abs(Y / N);
    P1 = P2(1:floor(N / 2) + 1);
    if numel(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end
    f = (0:floor(N / 2))' / N;
    magnitude = P1;
end

function label = snr_label(snr_db)
    if isinf(snr_db)
        label = 'Inf';
    else
        label = num2str(snr_db);
    end
end

function name = pretty_name(signalType)
    switch lower(signalType)
        case 'multiscale'
            name = 'Multiscale';
        case 'close'
            name = 'Close-frequency';
        otherwise
            name = signalType;
    end
end
