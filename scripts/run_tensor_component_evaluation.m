%RUN_TENSOR_COMPONENT_EVALUATION Evaluate candidate tensor CPD settings on synthetics.
%
%   Focused follow-up to the tensor fit sweep:
%   - reconstruct CP rank-1 terms back to 1D on the full signal grid
%   - compare grouping strategies for turning rank-1 terms into components
%   - evaluate against synthetic truth using existing matching/metric helpers
%
%   Writes:
%     results/tensor_component_eval_trials.csv
%     results/tensor_component_eval_summary.csv
%     results/tensor_component_eval_best.csv
%     results/tensor_component_eval_grouping_compare.csv
%     results/tensor_component_eval_results.mat

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));
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
config.signal_types = {'multiscale', 'close'};
config.snr_list = [Inf, 20, 10, 5];
config.seed_list = 11:15;
config.grouping_modes = {'none', 'benchmark', 'oracle', 'auto'};
config.candidate_settings = build_candidate_settings();

fprintf('\n=== Tensor component evaluation ===\n');
fprintf('Signal types : %s\n', strjoin(config.signal_types, ', '));
fprintf('SNR values   : %s\n', numvec_to_str(config.snr_list));
fprintf('Seeds        : %d:%d\n', config.seed_list(1), config.seed_list(end));
fprintf('Grouping     : %s\n', strjoin(config.grouping_modes, ', '));
fprintf('Candidate settings per signal are fixed in this script.\n');

num_trials = estimate_num_trials(config);
trial_rows = repmat(init_trial_row(), num_trials, 1);
row = 0;

for st = 1:numel(config.signal_types)
    signal_type = string(config.signal_types{st});
    settings = config.candidate_settings.(char(signal_type));
    for si = 1:numel(config.snr_list)
        snr_db = config.snr_list(si);
        for seed = config.seed_list
            [x, components] = generate_signal(signal_type, config.N, snr_db, seed);
            for ci = 1:numel(settings)
                for gi = 1:numel(config.grouping_modes)
                    row = row + 1;
                    trial_rows(row) = run_single_trial( ...
                        signal_type, snr_db, seed, x, components, settings(ci), ...
                        config.grouping_modes{gi}, config.fs);
                end
            end
        end
    end
end

trial_table = struct2table(trial_rows(1:row));
summary_table = summarize_by_condition(trial_table);
summary_table = sortrows(summary_table, {'signal_type', 'snr_db', 'grouping_mode', 'setting_order'}, {'ascend', 'descend', 'ascend', 'ascend'});
best_table = build_best_table(summary_table);
comparison_table = build_grouping_comparison_table(summary_table);

trialCsv = fullfile(resultsDir, 'tensor_component_eval_trials.csv');
summaryCsv = fullfile(resultsDir, 'tensor_component_eval_summary.csv');
bestCsv = fullfile(resultsDir, 'tensor_component_eval_best.csv');
comparisonCsv = fullfile(resultsDir, 'tensor_component_eval_grouping_compare.csv');
matFile = fullfile(resultsDir, 'tensor_component_eval_results.mat');

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
writetable(best_table, bestCsv);
writetable(comparison_table, comparisonCsv);
save(matFile, 'config', 'trial_table', 'summary_table', 'best_table', 'comparison_table');

fprintf('\nCondition summary:\n');
disp(summary_table);
fprintf('\nBest settings by signal/SNR:\n');
disp(best_table);
fprintf('\nGrouping comparison:\n');
disp(comparison_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    trialCsv, summaryCsv, bestCsv, comparisonCsv, matFile);

%% --- Local helpers ---
function n = estimate_num_trials(config)
    n = 0;
    for st = 1:numel(config.signal_types)
        signal_type = char(config.signal_types{st});
        n = n + numel(config.snr_list) * numel(config.seed_list) * ...
            numel(config.candidate_settings.(signal_type)) * numel(config.grouping_modes);
    end
end

function settings = build_candidate_settings()
    settings = struct();
    L_values = [200, 400, 600, 1000];
    R_values = [4, 6, 8];
    settings.multiscale = grid_settings("ms", L_values, R_values);
    settings.close = grid_settings("cl", L_values, R_values);
end

function settings = grid_settings(prefix, L_values, R_values)
    settings = repmat(setting_struct("", nan, [1, 2], nan, nan), 1, numel(L_values) * numel(R_values));
    order = 0;
    for i = 1:numel(L_values)
        for j = 1:numel(R_values)
            order = order + 1;
            L = L_values(i);
            R = R_values(j);
            settings(order) = setting_struct(sprintf("%s_L%d_R%d", prefix, L, R), order, [1, 2], L, R);
        end
    end
end

function s = setting_struct(label, order, strides, base_window_length, cpd_rank)
    s = struct( ...
        'label', string(label), ...
        'order', order, ...
        'strides', strides, ...
        'base_window_length', base_window_length, ...
        'cpd_rank', cpd_rank);
end

function row = init_trial_row()
    row = struct( ...
        'method', "TensorCPD", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'seed', nan, ...
        'setting_label', "", ...
        'setting_order', nan, ...
        'grouping_mode', "", ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'num_est_components', nan, ...
        'support_length', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'fit', nan, ...
        'rel_error', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan);
end

function row = run_single_trial(signal_type, snr_db, seed, x, components, setting, grouping_mode, fs)
    row = init_trial_row();
    row.signal_type = signal_type;
    row.snr_db = snr_db;
    row.seed = seed;
    row.setting_label = setting.label;
    row.setting_order = setting.order;
    row.grouping_mode = string(grouping_mode);
    row.stride_label = stride_label_from_values(setting.strides);
    row.base_window_length = setting.base_window_length;
    row.cpd_rank = setting.cpd_rank;

    reconstruction_mode = "cpd";

    params = tensor_merge_params(struct( ...
        'strides', setting.strides, ...
        'base_window_length', setting.base_window_length, ...
        'cpd_rank', setting.cpd_rank, ...
        'grouping_mode', grouping_mode, ...
        'reconstruction_mode', reconstruction_mode, ...
        'cpd_method', 'cpd', ...
        'cpd_options', struct()));

    try
        outcome = evaluate_tensor_cpd(signal_type, x, components, params, fs);
        row.success = outcome.success;
        row.failure_reason = string(outcome.failure_reason);
        row.support_length = outcome.support_length;
        row.num_est_components = size(outcome.est_components, 2);
        row.fit = outcome.fit;
        row.rel_error = outcome.rel_error;
        row.mean_nmse = outcome.mean_nmse;
        row.mean_abs_rho = outcome.mean_abs_rho;
        row.mean_peak_error = outcome.mean_peak_error;
    catch ME
        row.success = false;
        row.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
    end
end

function summary_table = summarize_by_condition(trial_table)
    keys = unique(trial_table(:, ...
        {'signal_type', 'snr_db', 'setting_label', 'setting_order', 'grouping_mode', 'stride_label', 'base_window_length', 'cpd_rank'}), ...
        'rows');
    rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            snr_equal(trial_table.snr_db, keys.snr_db(i)) & ...
            trial_table.setting_label == keys.setting_label(i) & ...
            trial_table.grouping_mode == keys.grouping_mode(i);
        subset = trial_table(mask, :);
        rows(i) = pack_summary_row(keys(i, :), subset);
    end
    summary_table = struct2table(rows);
end

function row = init_summary_row()
    row = struct( ...
        'method', "TensorCPD", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'setting_label', "", ...
        'setting_order', nan, ...
        'grouping_mode', "", ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_num_est_components', nan, ...
        'mean_fit', nan, ...
        'mean_rel_error', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_support_length', nan);
end

function row = pack_summary_row(key_row, subset)
    row = init_summary_row();
    row.signal_type = key_row.signal_type;
    row.snr_db = key_row.snr_db;
    row.setting_label = key_row.setting_label;
    row.setting_order = key_row.setting_order;
    row.grouping_mode = key_row.grouping_mode;
    row.stride_label = key_row.stride_label;
    row.base_window_length = key_row.base_window_length;
    row.cpd_rank = key_row.cpd_rank;
    row.num_trials = height(subset);
    row.failure_rate = mean(~subset.success);
    row.mean_num_est_components = mean(subset.num_est_components, 'omitnan');
    row.mean_fit = mean(subset.fit, 'omitnan');
    row.mean_rel_error = mean(subset.rel_error, 'omitnan');
    row.mean_nmse = mean(subset.mean_nmse, 'omitnan');
    row.mean_abs_rho = mean(subset.mean_abs_rho, 'omitnan');
    row.mean_peak_error = mean(subset.mean_peak_error, 'omitnan');
    row.mean_support_length = mean(subset.support_length, 'omitnan');
end

function best_table = build_best_table(summary_table)
    keys = unique(summary_table(:, {'signal_type', 'snr_db', 'grouping_mode'}), 'rows');
    rows = repmat(init_best_row(), height(keys), 1);
    for i = 1:height(keys)
        subset = summary_table(summary_table.signal_type == keys.signal_type(i) & ...
            snr_equal(summary_table.snr_db, keys.snr_db(i)) & ...
            summary_table.grouping_mode == keys.grouping_mode(i), :);
        rows(i) = choose_best_row(subset, keys.signal_type(i), keys.snr_db(i), keys.grouping_mode(i));
    end
    best_table = struct2table(rows);
    best_table = sortrows(best_table, {'signal_type', 'snr_db', 'grouping_mode'}, {'ascend', 'descend', 'ascend'});
end

function row = init_best_row()
    row = struct( ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'grouping_mode', "", ...
        'setting_label', "", ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'failure_rate', nan, ...
        'mean_num_est_components', nan, ...
        'mean_fit', nan, ...
        'mean_rel_error', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_support_length', nan);
end

function row = choose_best_row(subset, signal_type, snr_db, grouping_mode)
    row = init_best_row();
    row.signal_type = signal_type;
    row.snr_db = snr_db;
    row.grouping_mode = grouping_mode;
    if isempty(subset)
        return;
    end

    sort_subset = sortrows(subset, {'mean_nmse', 'mean_abs_rho', 'failure_rate', 'mean_fit'}, ...
        {'ascend', 'descend', 'ascend', 'descend'});
    best = sort_subset(1, :);
    row.setting_label = best.setting_label;
    row.stride_label = best.stride_label;
    row.base_window_length = best.base_window_length;
    row.cpd_rank = best.cpd_rank;
    row.failure_rate = best.failure_rate;
    row.mean_num_est_components = best.mean_num_est_components;
    row.mean_fit = best.mean_fit;
    row.mean_rel_error = best.mean_rel_error;
    row.mean_nmse = best.mean_nmse;
    row.mean_abs_rho = best.mean_abs_rho;
    row.mean_peak_error = best.mean_peak_error;
    row.mean_support_length = best.mean_support_length;
end

function comparison_table = build_grouping_comparison_table(summary_table)
    keys = unique(summary_table(:, ...
        {'signal_type', 'snr_db', 'setting_label', 'setting_order', 'stride_label', 'base_window_length', 'cpd_rank'}), ...
        'rows');
    rows = repmat(init_comparison_row(), height(keys), 1);
    modes = ["none", "benchmark", "oracle", "auto"];

    for i = 1:height(keys)
        rows(i).signal_type = keys.signal_type(i);
        rows(i).snr_db = keys.snr_db(i);
        rows(i).setting_label = keys.setting_label(i);
        rows(i).setting_order = keys.setting_order(i);
        rows(i).stride_label = keys.stride_label(i);
        rows(i).base_window_length = keys.base_window_length(i);
        rows(i).cpd_rank = keys.cpd_rank(i);

        subset = summary_table(summary_table.signal_type == keys.signal_type(i) & ...
            snr_equal(summary_table.snr_db, keys.snr_db(i)) & ...
            summary_table.setting_label == keys.setting_label(i), :);

        for m = 1:numel(modes)
            mode_name = modes(m);
            mode_subset = subset(subset.grouping_mode == mode_name, :);
            if isempty(mode_subset)
                continue;
            end
            rows(i) = assign_mode_metrics(rows(i), mode_name, mode_subset(1, :));
        end
    end

    comparison_table = struct2table(rows);
    comparison_table = sortrows(comparison_table, {'signal_type', 'snr_db', 'setting_order'}, {'ascend', 'descend', 'ascend'});
end

function row = init_comparison_row()
    row = struct( ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'setting_label', "", ...
        'setting_order', nan, ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'none_mean_nmse', nan, ...
        'benchmark_mean_nmse', nan, ...
        'oracle_mean_nmse', nan, ...
        'auto_mean_nmse', nan, ...
        'none_mean_abs_rho', nan, ...
        'benchmark_mean_abs_rho', nan, ...
        'oracle_mean_abs_rho', nan, ...
        'auto_mean_abs_rho', nan, ...
        'none_mean_fit', nan, ...
        'benchmark_mean_fit', nan, ...
        'oracle_mean_fit', nan, ...
        'auto_mean_fit', nan);
end

function row = assign_mode_metrics(row, mode_name, summary_row)
    prefix = char(mode_name);
    row.([prefix '_mean_nmse']) = summary_row.mean_nmse;
    row.([prefix '_mean_abs_rho']) = summary_row.mean_abs_rho;
    row.([prefix '_mean_fit']) = summary_row.mean_fit;
end

function label = stride_label_from_values(v)
    label = "S" + join(string(v), "_S");
end

function mask = snr_equal(v, target)
    if isinf(target)
        mask = isinf(v);
    else
        mask = v == target;
    end
end

function label = snr_label_scalar(v)
    if isinf(v)
        label = 'Inf';
    else
        label = num2str(v);
    end
end

function s = numvec_to_str(v)
    parts = cell(size(v));
    for i = 1:numel(v)
        parts{i} = snr_label_scalar(v(i));
    end
    s = strjoin(parts, ', ');
end
