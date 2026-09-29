%RUN_SLOW_COMPONENT_RECOVERY Analyze multiscale slow-component recovery.
%
%   This focused results script answers the slow-component question for the
%   multiscale synthetic benchmark:
%     - true slow component vs extracted slow component for SSA/SSD/TensorCPD
%     - slow-component correlation and RMSE by SNR
%     - a representative visual comparison of recovered slow components
%
%   Writes:
%     results/slow_component_recovery_trials.csv
%     results/slow_component_recovery_summary.csv
%     results/slow_component_recovery_plot_data.csv
%     results/slow_component_recovery_multiscale_<snr>db_seed<seed>.png

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

config = struct();
config.N = 2000;
config.fs = 1;
config.plot_snr_db = 20;
config.plot_seed = 11;

paths = struct();
paths.calibration_mat = fullfile(resultsDir, 'ssa_calibration_results.mat');
paths.ssa_trials_csv = fullfile(resultsDir, 'ssa_frozen_trials.csv');
paths.ssa_summary_csv = fullfile(resultsDir, 'ssa_frozen_summary.csv');
paths.ssd_trials_csv = fullfile(resultsDir, 'ssd_trials.csv');
paths.tensor_best_csv = fullfile(resultsDir, 'tensor_benchmark_best.csv');
paths.tensor_component_trials_csv = fullfile(resultsDir, 'tensor_benchmark_component_trials.csv');

require_file(paths.calibration_mat, 'scripts/run_ssa_calibration.m');
require_file(paths.ssa_trials_csv, 'scripts/run_ssa_frozen_experiments.m');
require_file(paths.ssa_summary_csv, 'scripts/run_ssa_frozen_experiments.m');
require_file(paths.ssd_trials_csv, 'scripts/run_ssd_experiments.m');
require_file(paths.tensor_best_csv, 'scripts/run_tensor_benchmark_experiments.m');
require_file(paths.tensor_component_trials_csv, 'scripts/run_tensor_benchmark_component_metrics.m');

[~, reference_components] = generate_signal('multiscale', config.N, Inf, 11);
slow_mean_square = mean(reference_components.slow(:).^2);

ssa_trials = readtable(paths.ssa_trials_csv, 'TextType', 'string');
ssa_summary = readtable(paths.ssa_summary_csv, 'TextType', 'string');
ssd_trials = readtable(paths.ssd_trials_csv, 'TextType', 'string');
tensor_trials = readtable(paths.tensor_component_trials_csv, 'TextType', 'string');

trial_table = build_slow_metric_table(ssa_trials, ssa_summary, ssd_trials, tensor_trials, slow_mean_square);
summary_table = summarize_slow_metrics(trial_table);

trialCsv = fullfile(resultsDir, 'slow_component_recovery_trials.csv');
summaryCsv = fullfile(resultsDir, 'slow_component_recovery_summary.csv');
writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);

[plot_table, plot_file] = make_recovery_plot(paths, config, resultsDir);
plotDataCsv = fullfile(resultsDir, 'slow_component_recovery_plot_data.csv');
writetable(plot_table, plotDataCsv);

fprintf('\n=== Slow-component recovery ===\n');
disp(summary_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', trialCsv, summaryCsv, plotDataCsv, plot_file);

%% --- Local helpers ---
function require_file(path_str, script_name)
    if ~exist(path_str, 'file')
        error('run_slow_component_recovery:MissingFile', ...
            'Missing %s. Run %s first.', path_str, script_name);
    end
end

function trial_table = build_slow_metric_table(ssa_trials, ssa_summary, ssd_trials, tensor_trials, slow_mean_square)
    ssa_slow = select_ssa_slow_rows(ssa_trials, ssa_summary);
    ssd_slow = ssd_trials(ssd_trials.signal_type == "multiscale", :);
    tensor_slow = tensor_trials(tensor_trials.signal_type == "multiscale" & ...
        tensor_trials.component_name == "slow", :);

    rows = [
        rows_from_wide_metrics(ssa_slow, "SSA", slow_mean_square); ...
        rows_from_wide_metrics(ssd_slow, "SSD (oscillatory-only)", slow_mean_square); ...
        rows_from_tensor_metrics(tensor_slow, slow_mean_square) ...
        ];
    trial_table = sortrows(rows, {'snr_db', 'method', 'seed'}, {'descend', 'ascend', 'ascend'});
end

function ssa_slow = select_ssa_slow_rows(ssa_trials, ssa_summary)
    ssa_trials = ssa_trials(ssa_trials.signal_type == "multiscale", :);
    ssa_summary = ssa_summary(ssa_summary.signal_type == "multiscale", :);
    keep = false(height(ssa_trials), 1);

    snrs = unique(ssa_trials.snr_db);
    for i = 1:numel(snrs)
        snr_db = snrs(i);
        summary_subset = ssa_summary(snr_equal(ssa_summary.snr_db, snr_db), :);
        if isempty(summary_subset)
            continue;
        end
        summary_subset = sortrows(summary_subset, ...
            {'mean_nmse', 'mean_abs_rho', 'L'}, {'ascend', 'descend', 'ascend'});
        best_L = summary_subset.L(1);
        keep = keep | (snr_equal(ssa_trials.snr_db, snr_db) & ssa_trials.L == best_L);
    end

    ssa_slow = ssa_trials(keep, :);
end

function rows = rows_from_wide_metrics(source, method_name, slow_mean_square)
    n = height(source);
    rows = init_trial_table(n);
    rows.method(:) = method_name;
    rows.signal_type(:) = "multiscale";
    rows.snr_db = source.snr_db;
    rows.seed = source.seed;
    rows.success = logical(source.success);
    rows.rho = source.comp2_rho;
    rows.abs_rho = abs(source.comp2_rho);
    rows.nmse = source.comp2_nmse;
    rows.rmse = sqrt(source.comp2_nmse .* slow_mean_square);
    if method_name == "SSA"
        rows.setting = "L=" + string(source.L);
    else
        rows.setting(:) = "adaptive L";
    end
end

function rows = rows_from_tensor_metrics(source, slow_mean_square)
    n = height(source);
    rows = init_trial_table(n);
    rows.method(:) = "TensorCPD";
    rows.signal_type(:) = "multiscale";
    rows.snr_db = source.snr_db;
    rows.seed = source.seed;
    rows.success = logical(source.success);
    rows.rho = source.rho;
    rows.abs_rho = source.abs_rho;
    rows.nmse = source.nmse;
    rows.rmse = sqrt(source.nmse .* slow_mean_square);
    rows.setting = source.setting_label;
end

function rows = init_trial_table(n)
    rows = table();
    rows.method = strings(n, 1);
    rows.signal_type = strings(n, 1);
    rows.snr_db = nan(n, 1);
    rows.seed = nan(n, 1);
    rows.setting = strings(n, 1);
    rows.success = false(n, 1);
    rows.rho = nan(n, 1);
    rows.abs_rho = nan(n, 1);
    rows.nmse = nan(n, 1);
    rows.rmse = nan(n, 1);
end

function summary_table = summarize_slow_metrics(trial_table)
    keys = unique(trial_table(:, {'method', 'signal_type', 'snr_db'}), 'rows');
    rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.method == keys.method(i) & ...
            trial_table.signal_type == keys.signal_type(i) & ...
            snr_equal(trial_table.snr_db, keys.snr_db(i));
        subset = trial_table(mask, :);
        rows(i) = pack_summary_row(keys(i, :), subset);
    end
    summary_table = struct2table(rows);
    summary_table = sortrows(summary_table, {'signal_type', 'snr_db', 'method'}, ...
        {'ascend', 'descend', 'ascend'});
end

function row = init_summary_row()
    row = struct( ...
        'method', "", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_rho', nan, ...
        'mean_abs_rho', nan, ...
        'mean_rmse', nan, ...
        'median_rmse', nan, ...
        'std_rmse', nan, ...
        'mean_nmse', nan);
end

function row = pack_summary_row(key_row, subset)
    row = init_summary_row();
    row.method = key_row.method;
    row.signal_type = key_row.signal_type;
    row.snr_db = key_row.snr_db;
    row.num_trials = height(subset);
    row.failure_rate = mean(~subset.success);
    row.mean_rho = mean(subset.rho, 'omitnan');
    row.mean_abs_rho = mean(subset.abs_rho, 'omitnan');
    row.mean_rmse = mean(subset.rmse, 'omitnan');
    row.median_rmse = median(subset.rmse, 'omitnan');
    row.std_rmse = std(subset.rmse, 'omitnan');
    row.mean_nmse = mean(subset.nmse, 'omitnan');
end

function [plot_table, plot_file] = make_recovery_plot(paths, config, resultsDir)
    loaded = load(paths.calibration_mat, 'best_params');
    best_params = loaded.best_params;
    ssa_summary = readtable(paths.ssa_summary_csv, 'TextType', 'string');
    tensor_best = readtable(paths.tensor_best_csv, 'TextType', 'string');

    [x, components, t] = generate_signal('multiscale', config.N, config.plot_snr_db, config.plot_seed);
    true_slow = components.slow(:);

    ssa_L = select_best_ssa_L(ssa_summary, config.plot_snr_db);
    ssa_features = ssa_precompute(x, ssa_L, best_params.r_max, config.fs);
    ssa_out = evaluate_auto_group_ssa('multiscale', ssa_features, components, best_params, config.fs);

    ssd_params = ssd_merge_params(best_params, struct('fs', config.fs));
    ssd_out = evaluate_ssd('multiscale', x, components, ssd_params, config.fs);

    tensor_params = tensor_params_from_best(tensor_best, config.plot_snr_db, config.fs);
    tensor_out = evaluate_tensor_cpd('multiscale', x, components, tensor_params, config.fs);

    methods = ["SSA"; "SSD (oscillatory-only)"; "TensorCPD"];
    estimates = [
        matched_component(ssa_out, 2, numel(true_slow)), ...
        matched_component(ssd_out, 2, numel(true_slow)), ...
        matched_component(tensor_out, 2, numel(true_slow)) ...
        ];
    rho = nan(numel(methods), 1);
    rmse = nan(numel(methods), 1);
    for i = 1:numel(methods)
        rho(i) = pearson_raw(true_slow, estimates(:, i));
        rmse(i) = sqrt(mean((true_slow - estimates(:, i)).^2, 'omitnan'));
    end

    plot_table = table();
    plot_table.n = t(:);
    plot_table.true_slow = true_slow;
    plot_table.ssa_slow = estimates(:, 1);
    plot_table.ssd_slow = estimates(:, 2);
    plot_table.tensorcpd_slow = estimates(:, 3);

    fig = figure('Name', 'Slow-component recovery', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [80 80 1300 850]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:numel(methods)
        ax = nexttile;
        plot(ax, t, true_slow, 'k-', 'LineWidth', 1.0, 'DisplayName', 'true slow');
        hold(ax, 'on');
        plot(ax, t, estimates(:, i), 'r--', 'LineWidth', 1.0, ...
            'DisplayName', char(methods(i) + " extracted slow"));
        title(ax, sprintf('%s slow recovery | rho = %.4f, RMSE = %.4g', ...
            char(methods(i)), rho(i), rmse(i)));
        xlabel(ax, 'n');
        ylabel(ax, 'amplitude');
        grid(ax, 'on');
        legend(ax, 'Location', 'best');
    end
    sgtitle(sprintf('True vs extracted slow component: multiscale, SNR = %s dB, seed = %d', ...
        snr_label(config.plot_snr_db), config.plot_seed));

    plot_file = fullfile(resultsDir, sprintf('slow_component_recovery_multiscale_%sdb_seed%d.png', ...
        lower(snr_label(config.plot_snr_db)), config.plot_seed));
    exportgraphics(fig, plot_file, 'Resolution', 300);
end

function L_best = select_best_ssa_L(ssa_summary, snr_db)
    subset = ssa_summary(ssa_summary.signal_type == "multiscale" & ...
        snr_equal(ssa_summary.snr_db, snr_db), :);
    if isempty(subset)
        error('run_slow_component_recovery:MissingSSASetting', ...
            'No SSA summary row found for multiscale SNR %s.', snr_label(snr_db));
    end
    subset = sortrows(subset, {'mean_nmse', 'mean_abs_rho', 'L'}, ...
        {'ascend', 'descend', 'ascend'});
    L_best = subset.L(1);
end

function params = tensor_params_from_best(tensor_best, snr_db, fs)
    subset = tensor_best(tensor_best.signal_type == "multiscale" & ...
        snr_equal(tensor_best.snr_db, snr_db), :);
    if isempty(subset)
        error('run_slow_component_recovery:MissingTensorSetting', ...
            'No TensorCPD best row found for multiscale SNR %s.', snr_label(snr_db));
    end
    row = subset(1, :);
    params = tensor_merge_params(struct( ...
        'strides', strides_from_label(row.stride_label), ...
        'base_window_length', row.base_window_length, ...
        'cpd_rank', row.cpd_rank, ...
        'grouping_mode', 'benchmark', ...
        'reconstruction_mode', 'cpd', ...
        'cpd_method', 'cpd', ...
        'cpd_options', struct(), ...
        'fs', fs));
end

function strides = strides_from_label(label)
    parts = split(erase(string(label), "S"), "_");
    strides = str2double(parts).';
    if any(isnan(strides))
        error('run_slow_component_recovery:InvalidStrideLabel', ...
            'Could not parse stride label "%s".', label);
    end
end

function y = matched_component(outcome, true_idx, n)
    y = nan(n, 1);
    if ~outcome.success || true_idx > numel(outcome.metrics.est_idx)
        return;
    end
    est_idx = outcome.metrics.est_idx(true_idx);
    if isnan(est_idx) || est_idx < 1 || est_idx > size(outcome.est_components, 2)
        return;
    end
    y = outcome.est_components(:, est_idx);
end

function c = pearson_raw(a, b)
    a = a(:);
    b = b(:);
    mask = ~(isnan(a) | isnan(b));
    a = a(mask) - mean(a(mask));
    b = b(mask) - mean(b(mask));
    den = norm(a) * norm(b);
    if den == 0
        c = nan;
    else
        c = (a.' * b) / den;
    end
end

function mask = snr_equal(v, target)
    if isinf(target)
        mask = isinf(v);
    else
        mask = v == target;
    end
end

function label = snr_label(v)
    if isinf(v)
        label = 'Inf';
    else
        label = num2str(v);
    end
end
