%RUN_TENSOR_BENCHMARK_COMPONENT_METRICS Summarize per-component TensorCPD metrics.
%
%   Replays the best TensorCPD setting selected for each benchmark condition
%   and writes one component-level table for the RQ2 discussion:
%     results/tensor_benchmark_component_trials.csv
%     results/tensor_benchmark_component_summary.csv
%     results/tensor_benchmark_component_metrics.mat

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));
addpath(fullfile(projRoot, 'tensor'));
if exist(fullfile(projRoot, 'tensorlab'), 'dir')
    addpath(genpath(fullfile(projRoot, 'tensorlab')));
end

resultsDir = fullfile(projRoot, 'results');
bestCsv = fullfile(resultsDir, 'tensor_benchmark_best.csv');
if ~exist(bestCsv, 'file')
    error('run_tensor_benchmark_component_metrics:MissingBestTable', ...
        'Run scripts/run_tensor_benchmark_experiments.m before this script.');
end

config = struct();
config.N = 2000;
config.fs = 1;
config.seed_list = 11:50;
config.grouping_mode = 'benchmark';
config.reconstruction_mode = 'cpd';

best_table = readtable(bestCsv, 'TextType', 'string');

fprintf('\n=== Tensor benchmark per-component metrics ===\n');
fprintf('Best settings : %s\n', bestCsv);
fprintf('Seeds         : %d:%d\n', config.seed_list(1), config.seed_list(end));

trial_rows = repmat(init_trial_row(), estimate_num_trials(best_table, config), 1);
row = 0;

for bi = 1:height(best_table)
    best = best_table(bi, :);
    signal_type = string(best.signal_type);
    snr_db = best.snr_db;
    params = params_from_best_row(best, config);

    for seed = config.seed_list
        [x, components] = generate_signal(signal_type, config.N, snr_db, seed);
        try
            outcome = evaluate_tensor_cpd(signal_type, x, components, params, config.fs);
            [trial_rows, row] = append_outcome_rows(trial_rows, row, best, seed, outcome);
        catch ME
            names = component_names_for_signal(signal_type);
            for ci = 1:numel(names)
                row = row + 1;
                trial_rows(row) = pack_failed_row(best, seed, ci, names{ci}, ...
                    sprintf('%s: %s', ME.identifier, ME.message));
            end
        end
    end
end

trial_table = struct2table(trial_rows(1:row));
summary_table = summarize_by_component(trial_table);
summary_table = sortrows(summary_table, ...
    {'signal_type', 'snr_db', 'component_order'}, {'ascend', 'descend', 'ascend'});

trialCsv = fullfile(resultsDir, 'tensor_benchmark_component_trials.csv');
summaryCsv = fullfile(resultsDir, 'tensor_benchmark_component_summary.csv');
matFile = fullfile(resultsDir, 'tensor_benchmark_component_metrics.mat');

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
save(matFile, 'config', 'best_table', 'trial_table', 'summary_table');

fprintf('\nPer-component summary:\n');
disp(summary_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n', trialCsv, summaryCsv, matFile);

%% --- Local helpers ---
function n = estimate_num_trials(best_table, config)
    n = 0;
    for i = 1:height(best_table)
        n = n + numel(config.seed_list) * numel(component_names_for_signal(best_table.signal_type(i)));
    end
end

function params = params_from_best_row(best, config)
    params = tensor_merge_params(struct( ...
        'strides', strides_from_label(best.stride_label), ...
        'base_window_length', best.base_window_length, ...
        'cpd_rank', best.cpd_rank, ...
        'grouping_mode', config.grouping_mode, ...
        'reconstruction_mode', config.reconstruction_mode, ...
        'cpd_method', 'cpd', ...
        'cpd_options', struct()));
end

function strides = strides_from_label(label)
    parts = split(erase(string(label), "S"), "_");
    strides = str2double(parts).';
    if any(isnan(strides))
        error('run_tensor_benchmark_component_metrics:InvalidStrideLabel', ...
            'Could not parse stride label "%s".', label);
    end
end

function names = component_names_for_signal(signal_type)
    switch lower(string(signal_type))
        case "multiscale"
            names = {'trend', 'slow', 'fast'};
        case "close"
            names = {'comp1', 'comp2'};
        otherwise
            error('run_tensor_benchmark_component_metrics:UnknownSignal', ...
                'Unknown signal type "%s".', signal_type);
    end
end

function row = init_trial_row()
    row = struct( ...
        'method', "TensorCPD", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'seed', nan, ...
        'setting_label', "", ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'component_order', nan, ...
        'component_name', "", ...
        'success', false, ...
        'failure_reason', "", ...
        'nmse', nan, ...
        'rho', nan, ...
        'abs_rho', nan, ...
        'peak_error', nan);
end

function row = pack_base_row(best, seed, component_order, component_name)
    row = init_trial_row();
    row.signal_type = string(best.signal_type);
    row.snr_db = best.snr_db;
    row.seed = seed;
    row.setting_label = string(best.setting_label);
    row.stride_label = string(best.stride_label);
    row.base_window_length = best.base_window_length;
    row.cpd_rank = best.cpd_rank;
    row.component_order = component_order;
    row.component_name = string(component_name);
end

function [trial_rows, row_index] = append_outcome_rows(trial_rows, row_index, best, seed, outcome)
    names = component_names_for_signal(best.signal_type);
    for ci = 1:numel(names)
        row_index = row_index + 1;
        value = pack_base_row(best, seed, ci, names{ci});
        value.success = outcome.success;
        value.failure_reason = string(outcome.failure_reason);
        if outcome.success && ci <= numel(outcome.metrics.nmse)
            value.nmse = outcome.metrics.nmse(ci);
            value.rho = outcome.metrics.pearson(ci);
            value.abs_rho = abs(outcome.metrics.pearson(ci));
            value.peak_error = outcome.metrics.peak_freq_error(ci);
        end
        trial_rows(row_index) = value;
    end
end

function row = pack_failed_row(best, seed, component_order, component_name, failure_reason)
    row = pack_base_row(best, seed, component_order, component_name);
    row.success = false;
    row.failure_reason = string(failure_reason);
end

function summary_table = summarize_by_component(trial_table)
    keys = unique(trial_table(:, ...
        {'signal_type', 'snr_db', 'setting_label', 'stride_label', ...
        'base_window_length', 'cpd_rank', 'component_order', 'component_name'}), ...
        'rows');
    rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            snr_equal(trial_table.snr_db, keys.snr_db(i)) & ...
            trial_table.component_order == keys.component_order(i);
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
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'component_order', nan, ...
        'component_name', "", ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_nmse', nan, ...
        'mean_rho', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan);
end

function row = pack_summary_row(key_row, subset)
    row = init_summary_row();
    row.signal_type = key_row.signal_type;
    row.snr_db = key_row.snr_db;
    row.setting_label = key_row.setting_label;
    row.stride_label = key_row.stride_label;
    row.base_window_length = key_row.base_window_length;
    row.cpd_rank = key_row.cpd_rank;
    row.component_order = key_row.component_order;
    row.component_name = key_row.component_name;
    row.num_trials = height(subset);
    row.failure_rate = mean(~subset.success);
    row.mean_nmse = mean(subset.nmse, 'omitnan');
    row.mean_rho = mean(subset.rho, 'omitnan');
    row.mean_abs_rho = mean(subset.abs_rho, 'omitnan');
    row.mean_peak_error = mean(subset.peak_error, 'omitnan');
end

function mask = snr_equal(v, target)
    if isinf(target)
        mask = isinf(v);
    else
        mask = v == target;
    end
end
