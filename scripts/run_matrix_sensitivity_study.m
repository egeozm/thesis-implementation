%RUN_MATRIX_SENSITIVITY_STUDY Controlled sensitivity study for matrix methods only.
%
%   Goal:
%   - show that SSA depends strongly on window length L
%   - show that SSD reduces some of that dependence
%   - keep both methods within single-resolution matrix Hankel embeddings
%
%   Requires results/ssa_calibration_results.mat from run_ssa_calibration.m.
%
%   Writes:
%     results/matrix_sensitivity_trials.csv
%     results/matrix_sensitivity_summary.csv
%     results/matrix_sensitivity_stability.csv
%     results/matrix_sensitivity_results.mat
%     results/matrix_sensitivity_nmse_by_setting.png
%     results/matrix_sensitivity_rho_by_setting.png
%     results/matrix_sensitivity_nmse_span.png

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));
addpath(fullfile(projRoot, 'ssd'));

resultsDir = fullfile(projRoot, 'results');
calibrationMat = fullfile(resultsDir, 'ssa_calibration_results.mat');
if ~exist(calibrationMat, 'file')
    error(['Missing calibration results. Run scripts/run_ssa_calibration.m first ', ...
        'to select and freeze the automatic grouping parameters.']);
end

loaded = load(calibrationMat, 'best_params');
best_params = loaded.best_params;

config = struct();
config.N = 2000;
config.fs = 1;
config.signal_types = {'multiscale', 'close'};
config.snr_list = [Inf, 20, 10, 5];
config.seed_list = 11:50;
config.ssa_L_list = [200, 400, 600, 1000];
config.ssd_setting_labels = {'L200', 'L400', 'L600', 'L1000', 'adaptive'};
config.ssd_L_sets = {
    [200], ...
    [400], ...
    [600], ...
    [1000], ...
    [200, 400, 600, 1000]
    };

fprintf('\n=== Matrix-method sensitivity study ===\n');
fprintf('Signal types : %s\n', strjoin(config.signal_types, ', '));
fprintf('SNR values   : %s\n', numvec_to_str(config.snr_list));
fprintf('Seeds        : %d:%d\n', config.seed_list(1), config.seed_list(end));
fprintf('SSA L values : %s\n', numvec_to_str(config.ssa_L_list));
fprintf('SSD settings : %s\n', strjoin(config.ssd_setting_labels, ', '));

num_trials = numel(config.signal_types) * numel(config.snr_list) * ...
    numel(config.seed_list) * (numel(config.ssa_L_list) + numel(config.ssd_L_sets));
trial_rows = repmat(init_trial_row(), num_trials, 1);
row = 0;

for st = 1:numel(config.signal_types)
    signal_type = config.signal_types{st};
    for si = 1:numel(config.snr_list)
        snr_db = config.snr_list(si);
        for seed = config.seed_list
            [x, components] = generate_signal(signal_type, config.N, snr_db, seed);

            for Li = 1:numel(config.ssa_L_list)
                L = config.ssa_L_list(Li);
                features = ssa_precompute(x, L, best_params.r_max, config.fs);
                outcome = evaluate_auto_group_ssa(signal_type, features, components, best_params, config.fs);
                row = row + 1;
                trial_rows(row) = build_ssa_trial_row(signal_type, snr_db, seed, L, outcome);
            end

            for ci = 1:numel(config.ssd_L_sets)
                ssd_params = ssd_merge_params(best_params, struct('L_list', config.ssd_L_sets{ci}));
                outcome = evaluate_ssd(signal_type, x, components, ssd_params, config.fs);
                row = row + 1;
                trial_rows(row) = build_ssd_trial_row(signal_type, snr_db, seed, ...
                    config.ssd_setting_labels{ci}, config.ssd_L_sets{ci}, outcome);
            end
        end
    end
end

trial_table = struct2table(trial_rows(1:row));
summary_table = summarize_by_condition(trial_table);
summary_table = sortrows(summary_table, {'signal_type', 'method', 'snr_db', 'setting_order'}, ...
    {'ascend', 'ascend', 'descend', 'ascend'});

stability_table = build_stability_table(summary_table);

trialCsv = fullfile(resultsDir, 'matrix_sensitivity_trials.csv');
summaryCsv = fullfile(resultsDir, 'matrix_sensitivity_summary.csv');
stabilityCsv = fullfile(resultsDir, 'matrix_sensitivity_stability.csv');
matFile = fullfile(resultsDir, 'matrix_sensitivity_results.mat');

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
writetable(stability_table, stabilityCsv);
save(matFile, 'config', 'best_params', 'trial_table', 'summary_table', 'stability_table');

plotFiles = plot_sensitivity_figures(summary_table, stability_table, resultsDir, config);

fprintf('\nCondition summary:\n');
disp(summary_table);
fprintf('\nStability summary (fixed settings only):\n');
disp(stability_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    trialCsv, summaryCsv, stabilityCsv, matFile, plotFiles{1}, plotFiles{2}, plotFiles{3});

%% --- Local helpers ---
function row = init_trial_row()
    row = struct( ...
        'method', "", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'seed', nan, ...
        'setting_label', "", ...
        'setting_kind', "", ...
        'setting_order', nan, ...
        'L', nan, ...
        'num_candidate_L', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'num_pairs', nan, ...
        'trend_group_size', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan);
end

function row = build_ssa_trial_row(signal_type, snr_db, seed, L, outcome)
    row = init_trial_row();
    row.method = "SSA";
    row.signal_type = string(signal_type);
    row.snr_db = snr_db;
    row.seed = seed;
    row.setting_label = "L" + string(L);
    row.setting_kind = "fixed";
    row.setting_order = setting_order_from_label(row.setting_label);
    row.L = L;
    row.num_candidate_L = 1;
    row.success = outcome.success;
    row.failure_reason = outcome.failure_reason;
    row.num_pairs = outcome.num_pairs;
    row.trend_group_size = outcome.trend_group_size;
    row.mean_nmse = outcome.mean_nmse;
    row.mean_abs_rho = outcome.mean_abs_rho;
    row.mean_peak_error = outcome.mean_peak_error;
end

function row = build_ssd_trial_row(signal_type, snr_db, seed, label, L_set, outcome)
    row = init_trial_row();
    row.method = "SSD";
    row.signal_type = string(signal_type);
    row.snr_db = snr_db;
    row.seed = seed;
    row.setting_label = string(label);
    row.setting_kind = ternary(numel(L_set) == 1, "fixed", "adaptive");
    row.setting_order = setting_order_from_label(row.setting_label);
    row.L = ternary(numel(L_set) == 1, L_set(1), nan);
    row.num_candidate_L = numel(L_set);
    row.success = outcome.success;
    row.failure_reason = outcome.failure_reason;
    row.num_pairs = outcome.num_ssd_modes;
    row.trend_group_size = nan;
    row.mean_nmse = outcome.mean_nmse;
    row.mean_abs_rho = outcome.mean_abs_rho;
    row.mean_peak_error = outcome.mean_peak_error;
end

function summary_table = summarize_by_condition(trial_table)
    keys = unique(trial_table(:, {'method', 'signal_type', 'snr_db', 'setting_label', 'setting_kind', 'setting_order'}), 'rows');
    rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.method == keys.method(i) & ...
            trial_table.signal_type == keys.signal_type(i) & ...
            trial_table.snr_db == keys.snr_db(i) & ...
            trial_table.setting_label == keys.setting_label(i);
        subset = trial_table(mask, :);
        rows(i) = pack_summary_row(keys(i, :), subset);
    end
    summary_table = struct2table(rows);
end

function row = init_summary_row()
    row = struct( ...
        'method', "", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'setting_label', "", ...
        'setting_kind', "", ...
        'setting_order', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_pairs', nan, ...
        'mean_trend_group_size', nan);
end

function row = pack_summary_row(key_row, subset)
    row = init_summary_row();
    row.method = key_row.method;
    row.signal_type = key_row.signal_type;
    row.snr_db = key_row.snr_db;
    row.setting_label = key_row.setting_label;
    row.setting_kind = key_row.setting_kind;
    row.setting_order = key_row.setting_order;
    row.num_trials = height(subset);
    row.failure_rate = mean(~subset.success);
    row.mean_nmse = mean(subset.mean_nmse, 'omitnan');
    row.mean_abs_rho = mean(subset.mean_abs_rho, 'omitnan');
    row.mean_peak_error = mean(subset.mean_peak_error, 'omitnan');
    row.mean_pairs = mean(subset.num_pairs, 'omitnan');
    row.mean_trend_group_size = mean(subset.trend_group_size, 'omitnan');
end

function stability_table = build_stability_table(summary_table)
    fixed_only = summary_table(summary_table.setting_kind == "fixed", :);
    keys = unique(fixed_only(:, {'method', 'signal_type', 'snr_db'}), 'rows');
    rows = repmat(init_stability_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = fixed_only.method == keys.method(i) & ...
            fixed_only.signal_type == keys.signal_type(i) & ...
            fixed_only.snr_db == keys.snr_db(i);
        subset = fixed_only(mask, :);
        rows(i) = pack_stability_row(keys.method(i), keys.signal_type(i), keys.snr_db(i), subset);
    end
    stability_table = struct2table(rows);
    stability_table = sortrows(stability_table, {'signal_type', 'snr_db', 'method'}, {'ascend', 'descend', 'ascend'});
end

function row = init_stability_row()
    row = struct( ...
        'method', "", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'num_settings', nan, ...
        'nmse_min', nan, ...
        'nmse_max', nan, ...
        'nmse_span', nan, ...
        'rho_min', nan, ...
        'rho_max', nan, ...
        'rho_span', nan);
end

function row = pack_stability_row(method, signal_type, snr_db, subset)
    row = init_stability_row();
    row.method = method;
    row.signal_type = signal_type;
    row.snr_db = snr_db;
    row.num_settings = height(subset);
    row.nmse_min = min(subset.mean_nmse, [], 'omitnan');
    row.nmse_max = max(subset.mean_nmse, [], 'omitnan');
    row.nmse_span = row.nmse_max - row.nmse_min;
    row.rho_min = min(subset.mean_abs_rho, [], 'omitnan');
    row.rho_max = max(subset.mean_abs_rho, [], 'omitnan');
    row.rho_span = row.rho_max - row.rho_min;
end

function plotFiles = plot_sensitivity_figures(summary_table, stability_table, resultsDir, config)
    signal_types = string(config.signal_types);
    snr_vals = sorted_snr(config.snr_list);

    nmseFig = figure('Name', 'Matrix sensitivity: mean NMSE by setting', 'NumberTitle', 'off');
    tiledlayout(numel(signal_types), 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signal_types)
        make_metric_panel(summary_table, "SSA", signal_types(si), snr_vals, 'mean_nmse', 'Mean NMSE', true);
        make_metric_panel(summary_table, "SSD", signal_types(si), snr_vals, 'mean_nmse', 'Mean NMSE', true);
    end
    sgtitle('Matrix methods: NMSE sensitivity to window configuration');
    nmsePng = fullfile(resultsDir, 'matrix_sensitivity_nmse_by_setting.png');
    exportgraphics(nmseFig, nmsePng, 'Resolution', 300);

    rhoFig = figure('Name', 'Matrix sensitivity: mean abs rho by setting', 'NumberTitle', 'off');
    tiledlayout(numel(signal_types), 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signal_types)
        make_metric_panel(summary_table, "SSA", signal_types(si), snr_vals, 'mean_abs_rho', 'Mean |Pearson \rho|', false);
        make_metric_panel(summary_table, "SSD", signal_types(si), snr_vals, 'mean_abs_rho', 'Mean |Pearson \rho|', false);
    end
    sgtitle('Matrix methods: correlation sensitivity to window configuration');
    rhoPng = fullfile(resultsDir, 'matrix_sensitivity_rho_by_setting.png');
    exportgraphics(rhoFig, rhoPng, 'Resolution', 300);

    spanFig = figure('Name', 'Matrix sensitivity: NMSE span across fixed settings', 'NumberTitle', 'off');
    tiledlayout(1, numel(signal_types), 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signal_types)
        ax = nexttile;
        hold(ax, 'on');
        subset = stability_table(stability_table.signal_type == signal_types(si), :);
        ssa_subset = subset(subset.method == "SSA", :);
        ssd_subset = subset(subset.method == "SSD", :);
        x = 1:numel(snr_vals);
        ssa_y = values_by_snr(ssa_subset, snr_vals, 'nmse_span');
        ssd_y = values_by_snr(ssd_subset, snr_vals, 'nmse_span');
        bar(ax, x - 0.18, ssa_y, 0.35, 'DisplayName', 'SSA');
        bar(ax, x + 0.18, ssd_y, 0.35, 'DisplayName', 'SSD (oscillatory-only)');
        set(ax, 'YScale', 'log');
        xticks(ax, x);
        xticklabels(ax, snr_labels(snr_vals));
        xlabel(ax, 'SNR (dB)');
        ylabel(ax, 'NMSE span across fixed settings');
        title(ax, char(signal_types(si)));
        grid(ax, 'on');
        legend(ax, 'Location', 'best');
    end
    sgtitle('Matrix methods: fixed-setting sensitivity span');
    spanPng = fullfile(resultsDir, 'matrix_sensitivity_nmse_span.png');
    exportgraphics(spanFig, spanPng, 'Resolution', 300);

    plotFiles = {nmsePng, rhoPng, spanPng};
end

function make_metric_panel(summary_table, method, signal_type, snr_vals, field_name, y_label, use_log)
    ax = nexttile;
    hold(ax, 'on');
    subset = summary_table(summary_table.method == method & summary_table.signal_type == signal_type, :);
    subset = sortrows(subset, {'setting_order'}, {'ascend'});
    setting_labels = unique(subset(:, {'setting_label', 'setting_order'}), 'rows');
    setting_labels = sortrows(setting_labels, 'setting_order');
    x = 1:height(setting_labels);
    for si = 1:numel(snr_vals)
        y = nan(size(x));
        for k = 1:numel(x)
            rowMask = subset.setting_label == setting_labels.setting_label(k) & snr_equal(subset.snr_db, snr_vals(si));
            if any(rowMask)
                y(k) = subset.(field_name)(find(rowMask, 1));
            end
        end
        plot(ax, x, y, 'o-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('SNR=%s', snr_label_scalar(snr_vals(si))));
    end
    if use_log
        set(ax, 'YScale', 'log');
    else
        ylim(ax, [0, 1.05]);
    end
    xticks(ax, x);
    xticklabels(ax, cellstr(setting_labels.setting_label));
    xlabel(ax, ternary(method == "SSA", 'Fixed L', 'Oscillatory-only SSD L setting'));
    ylabel(ax, y_label);
    title(ax, sprintf('%s - %s', signal_type, display_method_label(method)));
    grid(ax, 'on');
    legend(ax, 'Location', 'best');
end

function y = values_by_snr(subset, snr_vals, field_name)
    y = nan(size(snr_vals));
    for i = 1:numel(snr_vals)
        rowMask = snr_equal(subset.snr_db, snr_vals(i));
        if any(rowMask)
            y(i) = subset.(field_name)(find(rowMask, 1));
        end
    end
end

function label = display_method_label(method)
    if method == "SSD"
        label = "SSD (oscillatory-only)";
    else
        label = method;
    end
end

function order = setting_order_from_label(label)
    switch string(label)
        case "L200"
            order = 1;
        case "L400"
            order = 2;
        case "L600"
            order = 3;
        case "L1000"
            order = 4;
        case "adaptive"
            order = 5;
        otherwise
            error('run_matrix_sensitivity_study:UnknownSetting', 'Unknown setting label %s', label);
    end
end

function vals = sorted_snr(v)
    v = unique(v(:).');
    finiteVals = sort(v(isfinite(v)));
    if any(isinf(v))
        vals = [finiteVals, Inf];
    else
        vals = finiteVals;
    end
end

function labels = snr_labels(v)
    labels = cell(size(v));
    for i = 1:numel(v)
        labels{i} = snr_label_scalar(v(i));
    end
end

function label = snr_label_scalar(v)
    if isinf(v)
        label = 'Inf';
    else
        label = num2str(v);
    end
end

function mask = snr_equal(v, target)
    if isinf(target)
        mask = isinf(v);
    else
        mask = v == target;
    end
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end

function s = numvec_to_str(v)
    parts = cell(size(v));
    for i = 1:numel(v)
        parts{i} = snr_label_scalar(v(i));
    end
    s = strjoin(parts, ', ');
end
