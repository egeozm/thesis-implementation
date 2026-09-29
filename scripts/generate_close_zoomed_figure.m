%GENERATE_CLOSE_ZOOMED_FIGURE Generate zoomed close-frequency thesis figure.
%
%   Exports:
%     results/representative_close_inf_zoomed_4x1.png

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));
addpath(fullfile(projRoot, 'ssd'));
addpath(fullfile(projRoot, 'tensor'));
if exist(fullfile(projRoot, 'tensorlab'), 'dir')
    addpath(genpath(fullfile(projRoot, 'tensorlab')));
end

resultsDir = fullfile(projRoot, 'results');
loaded = load(fullfile(resultsDir, 'ssa_calibration_results.mat'), 'best_params');
best_params = loaded.best_params;
ssa_summary = readtable(fullfile(resultsDir, 'ssa_frozen_summary.csv'), 'TextType', 'string');
tensor_best = readtable(fullfile(resultsDir, 'tensor_benchmark_best.csv'), 'TextType', 'string');

N = 2000;
fs = 1;
seed = 11;
snr_db = Inf;
signal_type = "close";
zoom_idx = 1:200;

[x, components, t] = generate_signal(signal_type, N, snr_db, seed);
truth = [components.comp1, components.comp2];

ssa_L = select_best_ssa_L(ssa_summary, signal_type, snr_db);
ssa_features = ssa_precompute(x, ssa_L, best_params.r_max, fs);
ssa_out = evaluate_auto_group_ssa(signal_type, ssa_features, components, best_params, fs);
ssa_components = matched_components(ssa_out, 2);

ssd_params = ssd_merge_params(best_params, struct('fs', fs));
ssd_out = evaluate_ssd(signal_type, x, components, ssd_params, fs);
ssd_components = matched_components(ssd_out, 2);

tensor_params = tensor_params_from_best(tensor_best, signal_type, snr_db, fs);
tensor_out = evaluate_tensor_cpd(signal_type, x, components, tensor_params, fs);
tensor_components = matched_components(tensor_out, 2);

fig = figure('Name', 'Close-frequency zoomed 4x1', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1350 1450]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'loose');

ax = nexttile;
plot(ax, t(zoom_idx), x(zoom_idx), 'Color', [0.75 0.75 0.75], ...
    'LineWidth', 1.2, 'DisplayName', 'mixture');
hold(ax, 'on');
plot(ax, t(zoom_idx), truth(zoom_idx, 1), 'k-', ...
    'LineWidth', 1.7, 'DisplayName', 'true comp1');
plot(ax, t(zoom_idx), truth(zoom_idx, 2), '-', 'Color', [0 0.25 0.75], ...
    'LineWidth', 1.7, 'DisplayName', 'true comp2');
title(ax, 'Ground truth: two nearby-frequency sinusoids');
format_axis(ax, 'amplitude');
format_legend(legend(ax, 'Location', 'northeast', 'FontSize', 16));

plot_method_panel(nexttile, t, zoom_idx, truth, ssa_components, ...
    sprintf('SSA reconstruction (L = %d): mode mixing remains', ssa_L), true);
plot_method_panel(nexttile, t, zoom_idx, truth, ssd_components, ...
    'SSD reconstruction: improved separation with amplitude distortion', false);
plot_method_panel(nexttile, t, zoom_idx, truth, tensor_components, ...
    sprintf('TensorCPD reconstruction (L = %d, R = %d): closest visual match', ...
    tensor_params.base_window_length, tensor_params.cpd_rank), false);

fig_title = sgtitle('Zoomed close-frequency reconstruction, SNR = Inf dB, seed = 11', ...
    'FontSize', 21, 'FontWeight', 'bold');
fig_title.Color = 'k';

axes_handles = findall(fig, 'Type', 'axes');
for i = 1:numel(axes_handles)
    if isprop(axes_handles(i), 'Toolbar') && ~isempty(axes_handles(i).Toolbar)
        axes_handles(i).Toolbar.Visible = 'off';
    end
end

outFile = fullfile(resultsDir, 'representative_close_inf_zoomed_4x1.png');
exportgraphics(fig, outFile, 'Resolution', 300);
close(fig);
fprintf('Saved %s\n', outFile);

function plot_method_panel(ax, t, zoom_idx, truth, estimate, panel_title, show_legend)
    plot(ax, t(zoom_idx), truth(zoom_idx, 1), 'k-', ...
        'LineWidth', 1.5, 'DisplayName', 'true comp1');
    hold(ax, 'on');
    plot(ax, t(zoom_idx), truth(zoom_idx, 2), '-', 'Color', [0 0.25 0.75], ...
        'LineWidth', 1.5, 'DisplayName', 'true comp2');
    plot(ax, t(zoom_idx), estimate(zoom_idx, 1), 'r--', ...
        'LineWidth', 1.5, 'DisplayName', 'recon comp1');
    plot(ax, t(zoom_idx), estimate(zoom_idx, 2), '--', 'Color', [0.85 0.45 0], ...
        'LineWidth', 1.5, 'DisplayName', 'recon comp2');
    title(ax, panel_title);
    format_axis(ax, 'amplitude');
    if show_legend
        format_legend(legend(ax, 'Location', 'northeast', 'FontSize', 16));
    end
end

function format_axis(ax, y_label)
    xlabel(ax, 'sample index n');
    ylabel(ax, y_label);
    grid(ax, 'on');
    xlim(ax, [1 200]);
    set(ax, 'FontSize', 17, 'LineWidth', 1.1, ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.75 0.75 0.75], 'GridAlpha', 0.45);
    ax.Title.FontSize = 18;
    ax.Title.FontWeight = 'bold';
    ax.Title.Color = 'k';
    ax.XLabel.FontSize = 17;
    ax.XLabel.Color = 'k';
    ax.YLabel.FontSize = 17;
    ax.YLabel.Color = 'k';
end

function format_legend(lgd)
    lgd.Color = 'w';
    lgd.TextColor = 'k';
    lgd.EdgeColor = [0.8 0.8 0.8];
end

function L_best = select_best_ssa_L(ssa_summary, signal_type, snr_db)
    mask = string(ssa_summary.signal_type) == string(signal_type) & snr_equal(ssa_summary.snr_db, snr_db);
    subset = ssa_summary(mask, :);
    if isempty(subset)
        error('generate_close_zoomed_figure:MissingSSASetting', ...
            'No SSA setting found for signal=%s.', signal_type);
    end
    subset = sortrows(subset, {'mean_nmse', 'mean_abs_rho', 'L'}, {'ascend', 'descend', 'ascend'});
    L_best = subset.L(1);
end

function params = tensor_params_from_best(tensor_best, signal_type, snr_db, fs)
    mask = string(tensor_best.signal_type) == string(signal_type) & snr_equal(tensor_best.snr_db, snr_db);
    subset = tensor_best(mask, :);
    if isempty(subset)
        error('generate_close_zoomed_figure:MissingTensorSetting', ...
            'No tensor best setting found for signal=%s.', signal_type);
    end
    row = subset(1, :);
    params = tensor_merge_params(struct( ...
        'strides', parse_stride_label(row.stride_label(1)), ...
        'base_window_length', row.base_window_length, ...
        'cpd_rank', row.cpd_rank, ...
        'grouping_mode', 'benchmark', ...
        'reconstruction_mode', 'cpd', ...
        'cpd_method', 'cpd', ...
        'cpd_options', struct(), ...
        'fs', fs));
end

function strides = parse_stride_label(label)
    parts = split(string(label), "_");
    strides = zeros(1, numel(parts));
    for i = 1:numel(parts)
        strides(i) = str2double(erase(parts(i), "S"));
    end
end

function matched = matched_components(outcome, num_truth)
    matched = zeros(size(outcome.est_components, 1), num_truth);
    for k = 1:num_truth
        idx = outcome.metrics.est_idx(k);
        if isfinite(idx) && idx >= 1 && idx <= size(outcome.est_components, 2)
            matched(:, k) = outcome.est_components(:, idx);
        end
    end
end

function tf = snr_equal(values, target)
    if isinf(target)
        tf = isinf(values);
    else
        tf = values == target;
    end
end
