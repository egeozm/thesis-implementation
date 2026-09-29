%RUN_REPRESENTATIVE_RECONSTRUCTIONS Visualize representative SSA/SSD/Tensor runs.
%
%   Generates a small set of representative signal-space figures for the three
%   methods used in the thesis comparison:
%     - SSA (frozen automatic grouping)
%     - SSD (adaptive residual extraction)
%     - Tensor CPD (finalized benchmark path)
%
%   For each selected case, the script exports:
%     - a component-by-component comparison figure
%     - a signal-flow figure showing input / internal signal objects / outputs

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));
addpath(fullfile(projRoot, 'ssd'));
addpath(fullfile(projRoot, 'tensor'));
if exist(fullfile(projRoot, 'tensorlab'), 'dir')
    addpath(genpath(fullfile(projRoot, 'tensorlab')));
end

resultsDir = fullfile(projRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

calibrationMat = fullfile(resultsDir, 'ssa_calibration_results.mat');
ssaSummaryCsv = fullfile(resultsDir, 'ssa_frozen_summary.csv');
tensorBestCsv = fullfile(resultsDir, 'tensor_benchmark_best.csv');

require_file(calibrationMat, 'scripts/run_ssa_calibration.m');
require_file(ssaSummaryCsv, 'scripts/run_ssa_frozen_experiments.m');
require_file(tensorBestCsv, 'scripts/run_tensor_benchmark_experiments.m');

loaded = load(calibrationMat, 'best_params');
best_params = loaded.best_params;
ssa_summary = readtable(ssaSummaryCsv);
tensor_best = readtable(tensorBestCsv);

config = struct();
config.N = 2000;
config.fs = 1;
config.cases = [ ...
    case_struct("close", Inf, 11, "close_inf"), ...
    case_struct("close", 5, 11, "close_5db"), ...
    case_struct("multiscale", 20, 11, "multiscale_20db"), ...
    case_struct("multiscale", 5, 11, "multiscale_5db") ...
    ];

for i = 1:numel(config.cases)
    case_i = config.cases(i);
    [x, components, t] = generate_signal(case_i.signal_type, config.N, case_i.snr_db, case_i.seed);
    [true_components, component_names] = truth_components(case_i.signal_type, components);

    ssa_L = select_best_ssa_L(ssa_summary, case_i.signal_type, case_i.snr_db);
    ssa_features = ssa_precompute(x, ssa_L, best_params.r_max, config.fs);
    ssa_out = evaluate_auto_group_ssa(case_i.signal_type, ssa_features, components, best_params, config.fs);

    ssd_params = ssd_merge_params(best_params, struct('fs', config.fs));
    ssd_out = evaluate_ssd(case_i.signal_type, x, components, ssd_params, config.fs);

    tensor_params = tensor_params_from_best(tensor_best, case_i.signal_type, case_i.snr_db, config.fs);
    tensor_out = evaluate_tensor_cpd(case_i.signal_type, x, components, tensor_params, config.fs);

    methods = build_method_structs(x, true_components, component_names, t, ...
        ssa_out, ssa_features, ssa_L, ssd_out, tensor_out);

    comparison_fig = make_component_comparison_figure(case_i, methods, true_components, component_names, t);
    flow_fig = make_signal_flow_figure(case_i, methods, x, components, t);

    comparison_file = fullfile(resultsDir, sprintf('representative_%s_components.png', case_i.export_tag));
    flow_file = fullfile(resultsDir, sprintf('representative_%s_signal_flow.png', case_i.export_tag));
    exportgraphics(comparison_fig, comparison_file, 'Resolution', 300);
    exportgraphics(flow_fig, flow_file, 'Resolution', 300);
end

fprintf('\nRepresentative reconstruction visualization finished.\n');

%% --- Local helpers ---
function require_file(path_str, script_name)
    if ~exist(path_str, 'file')
        error('run_representative_reconstructions:MissingFile', ...
            'Missing %s. Run %s first.', path_str, script_name);
    end
end

function c = case_struct(signal_type, snr_db, seed, export_tag)
    c = struct( ...
        'signal_type', string(signal_type), ...
        'snr_db', snr_db, ...
        'seed', seed, ...
        'export_tag', string(export_tag));
end

function [true_components, component_names] = truth_components(signal_type, components)
    switch lower(string(signal_type))
        case "multiscale"
            true_components = [components.trend, components.slow, components.fast];
            component_names = {"trend", "slow", "fast"};
        case "close"
            true_components = [components.comp1, components.comp2];
            component_names = {"comp1", "comp2"};
        otherwise
            error('run_representative_reconstructions:UnknownSignal', ...
                'Unknown signal type %s.', signal_type);
    end
end

function L_best = select_best_ssa_L(ssa_summary, signal_type, snr_db)
    signal_mask = ssa_summary.signal_type == string(signal_type);
    snr_mask = snr_equal_table(ssa_summary.snr_db, snr_db);
    subset = ssa_summary(signal_mask & snr_mask, :);
    if isempty(subset)
        error('run_representative_reconstructions:MissingSSASetting', ...
            'No SSA summary row found for signal=%s, snr=%s.', signal_type, snr_label(snr_db));
    end
    subset = sortrows(subset, {'mean_nmse', 'mean_abs_rho', 'L'}, {'ascend', 'descend', 'ascend'});
    L_best = subset.L(1);
end

function params = tensor_params_from_best(tensor_best, signal_type, snr_db, fs)
    signal_mask = tensor_best.signal_type == string(signal_type);
    snr_mask = snr_equal_table(tensor_best.snr_db, snr_db);
    subset = tensor_best(signal_mask & snr_mask, :);
    if isempty(subset)
        error('run_representative_reconstructions:MissingTensorSetting', ...
            'No tensor best row found for signal=%s, snr=%s.', signal_type, snr_label(snr_db));
    end
    row = subset(1, :);
    params = tensor_merge_params(struct( ...
        'strides', parse_stride_label(row.stride_label{1}), ...
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

function tf = snr_equal_table(values, target)
    if isinf(target)
        tf = isinf(values);
    else
        tf = values == target;
    end
end

function methods = build_method_structs(x, true_components, component_names, t, ...
        ssa_out, ssa_features, ssa_L, ssd_out, tensor_out)
    methods = repmat(init_method_struct(), 3, 1);
    methods(1) = build_ssa_method(x, true_components, component_names, t, ssa_out, ssa_features, ssa_L);
    methods(2) = build_ssd_method(x, true_components, component_names, t, ssd_out);
    methods(3) = build_tensor_method(x, true_components, component_names, t, tensor_out);
end

function method = init_method_struct()
    method = struct( ...
        'name', "", ...
        'input_signal', [], ...
        'input_label', "", ...
        'internal_matrix', [], ...
        'internal_names', {{}}, ...
        'internal_title', "", ...
        'internal_aux', struct(), ...
        'output_components', [], ...
        'output_names', {{}}, ...
        'output_residual', [], ...
        'subtitle', "");
end

function method = build_ssa_method(x, true_components, component_names, t, outcome, features, L_best)
    method = init_method_struct();
    method.name = "SSA";
    method.input_signal = x(:);
    method.input_label = "observed";
    num_truth = size(true_components, 2);
    method.output_components = matched_components(outcome, num_truth);
    method.output_names = component_names;
    method.output_residual = x(:) - sum(method.output_components, 2);
    r_show = min(6, size(features.elementary, 2));
    method.internal_matrix = features.elementary(:, 1:r_show);
    method.internal_names = arrayfun(@(r) sprintf('elem%d', r), 1:r_show, 'UniformOutput', false);
    method.internal_title = sprintf('Leading elementary signals (L = %d)', L_best);
    method.internal_aux = struct('t', t(:));
    method.subtitle = sprintf('L = %d, mean NMSE = %.4g', L_best, outcome.mean_nmse);
end

function method = build_ssd_method(x, true_components, component_names, t, outcome)
    method = init_method_struct();
    method.name = "SSD (oscillatory-only)";
    method.input_signal = x(:);
    method.input_label = "observed";
    num_truth = size(true_components, 2);
    method.output_components = matched_components(outcome, num_truth);
    method.output_names = component_names;
    method.output_residual = outcome.decomposition.residual(:);
    num_modes_show = min(3, size(outcome.decomposition.modes, 2));
    method.internal_matrix = [outcome.decomposition.modes(:, 1:num_modes_show), outcome.decomposition.residual(:)];
    method.internal_names = [ ...
        arrayfun(@(r) sprintf('mode%d', r), 1:num_modes_show, 'UniformOutput', false), ...
        {'residual'}];
    method.internal_title = 'Extracted modes and final residual';
    method.internal_aux = struct('t', t(:), 'L_history', outcome.decomposition.L_history);
    method.subtitle = sprintf('%d extracted modes, mean NMSE = %.4g', ...
        outcome.decomposition.num_iterations, outcome.mean_nmse);
end

function method = build_tensor_method(x, true_components, component_names, t, outcome)
    method = init_method_struct();
    method.name = "Tensor CPD";
    method.input_signal = x(:);
    method.input_label = "observed";
    num_truth = size(true_components, 2);
    method.output_components = matched_components(outcome, num_truth);
    method.output_names = component_names;
    method.output_residual = x(:) - sum(method.output_components, 2);
    if isfield(outcome.grouping, 'total_reconstruction')
        total_recon = outcome.grouping.total_reconstruction(:);
    else
        total_recon = sum(method.output_components, 2);
    end
    method.internal_matrix = total_recon;
    method.internal_names = {'x_hat'};
    method.internal_title = 'CPD reconstruction and per-stride signals';
    method.internal_aux = struct('t', t(:), 'per_stride', outcome.reconstruction.per_stride);
    method.subtitle = sprintf('fit = %.3f, mean NMSE = %.4g', outcome.fit, outcome.mean_nmse);
end

function matched = matched_components(outcome, num_truth)
    matched = zeros(size(outcome.est_components, 1), num_truth);
    for k = 1:num_truth
        idx = outcome.metrics.est_idx(k);
        if isnan(idx) || idx < 1 || idx > size(outcome.est_components, 2)
            continue;
        end
        matched(:, k) = outcome.est_components(:, idx);
    end
end

function fig = make_component_comparison_figure(case_i, methods, true_components, component_names, t)
    num_components = size(true_components, 2);
    fig = figure('Name', sprintf('Representative components: %s', case_i.export_tag), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [50 50 1450 900]);
    tiledlayout(num_components, numel(methods), 'Padding', 'compact', 'TileSpacing', 'compact');
    for c = 1:num_components
        for m = 1:numel(methods)
            ax = nexttile;
            plot(ax, t, true_components(:, c), 'k-', 'LineWidth', 1.0, 'DisplayName', 'truth');
            hold(ax, 'on');
            plot(ax, t, methods(m).output_components(:, c), 'r--', 'LineWidth', 1.0, 'DisplayName', char(methods(m).name));
            title(ax, sprintf('%s | %s', char(methods(m).name), component_names{c}));
            xlabel(ax, 'n');
            ylabel(ax, 'amplitude');
            grid(ax, 'on');
            if c == 1 && m == numel(methods)
                legend(ax, 'Location', 'best');
            end
        end
    end
    sgtitle(sprintf('Representative component reconstructions: %s, SNR = %s dB, seed = %d', ...
        pretty_name(case_i.signal_type), snr_label(case_i.snr_db), case_i.seed));
end

function fig = make_signal_flow_figure(case_i, methods, x, components, t)
    fig = figure('Name', sprintf('Representative flow: %s', case_i.export_tag), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1550 950]);
    tiledlayout(numel(methods), 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    clean_signal = components.x_clean(:);

    for m = 1:numel(methods)
        ax_in = nexttile;
        plot(ax_in, t, x(:), 'Color', [0.20 0.20 0.20], 'LineWidth', 1.0, 'DisplayName', methods(m).input_label);
        hold(ax_in, 'on');
        plot(ax_in, t, clean_signal, 'b--', 'LineWidth', 1.0, 'DisplayName', 'clean');
        title(ax_in, sprintf('%s | input', char(methods(m).name)));
        xlabel(ax_in, 'n');
        ylabel(ax_in, 'amplitude');
        grid(ax_in, 'on');
        if m == 1
            legend(ax_in, 'Location', 'best');
        end

        ax_mid = nexttile;
        plot_internal_panel(ax_mid, methods(m), t, x);

        ax_out = nexttile;
        plot_output_panel(ax_out, methods(m), t);
    end

    sgtitle(sprintf('Signal-flow view: %s, SNR = %s dB, seed = %d', ...
        pretty_name(case_i.signal_type), snr_label(case_i.snr_db), case_i.seed));
end

function plot_internal_panel(ax, method, t, x)
    switch char(method.name)
        case 'SSA'
            plot_offset_series(ax, t, method.internal_matrix, method.internal_names);
            title(ax, sprintf('SSA inside: %s', method.internal_title));
            xlabel(ax, 'n');
            ylabel(ax, 'offset signals');
            grid(ax, 'on');
        case 'SSD (oscillatory-only)'
            plot_offset_series(ax, t, method.internal_matrix, method.internal_names);
            title(ax, 'Oscillatory-only SSD inside: extracted modes -> residual');
            xlabel(ax, 'n');
            ylabel(ax, 'offset signals');
            grid(ax, 'on');
        otherwise
            plot(ax, t, x(:), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'x');
            hold(ax, 'on');
            plot(ax, t, method.internal_matrix(:), 'r-', 'LineWidth', 1.1, 'DisplayName', 'xhat');
            per_stride = method.internal_aux.per_stride;
            colors = lines(max(numel(per_stride), 1));
            for i = 1:numel(per_stride)
                if ~per_stride(i).success
                    continue;
                end
                plot(ax, per_stride(i).sample_indices(:), per_stride(i).total_decimated_reconstruction(:), ...
                    '.', 'Color', colors(i, :), 'MarkerSize', 8, ...
                    'DisplayName', sprintf('stride %d', per_stride(i).stride));
            end
            title(ax, 'Tensor inside: CPD total and stride reconstructions');
            xlabel(ax, 'n');
            ylabel(ax, 'amplitude');
            legend(ax, 'Location', 'best');
            grid(ax, 'on');
    end
end

function plot_output_panel(ax, method, t)
    Y = [method.output_components, method.output_residual];
    names = [method.output_names, {'residual'}];
    plot_offset_series(ax, t, Y, names);
    title(ax, sprintf('%s output | %s', char(method.name), method.subtitle));
    xlabel(ax, 'n');
    ylabel(ax, 'offset signals');
    grid(ax, 'on');
end

function plot_offset_series(ax, t, Y, names)
    Y = Y(:,:);
    n_series = size(Y, 2);
    if n_series == 0
        title(ax, 'No signals to show');
        axis(ax, 'off');
        return;
    end
    colors = lines(n_series);
    amp = max(abs(Y(:)));
    if amp <= eps
        amp = 1;
    end
    spacing = 2.5 * amp;
    hold(ax, 'on');
    for k = 1:n_series
        offset = (n_series - k) * spacing;
        plot(ax, t, Y(:, k) + offset, 'LineWidth', 1.0, 'Color', colors(k, :));
        text(ax, t(1), offset, ['  ' names{k}], 'Color', colors(k, :), ...
            'VerticalAlignment', 'bottom', 'Interpreter', 'none');
    end
    yticks(ax, []);
end

function label = snr_label(v)
    if isinf(v)
        label = 'Inf';
    else
        label = num2str(v);
    end
end

function name = pretty_name(signal_type)
    switch lower(string(signal_type))
        case "multiscale"
            name = 'Multiscale';
        case "close"
            name = 'Close-frequency';
        otherwise
            name = char(signal_type);
    end
end
