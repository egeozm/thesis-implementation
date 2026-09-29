%RUN_TENSOR_BENCHMARK_EXPERIMENTS Evaluate tensor CPD on thesis synthetics.
%
%   Uses the benchmark-oriented grouping mode so tensor outputs can be
%   compared directly against the SSA/SSD synthetic benchmark tables.
%
%   Writes:
%     results/tensor_benchmark_trials.csv
%     results/tensor_benchmark_summary.csv
%     results/tensor_benchmark_best.csv
%     results/tensor_benchmark_results.mat

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
config.grouping_mode = 'benchmark';
config.reconstruction_mode = 'cpd';
config.signal_types = {'multiscale', 'close'};
config.snr_list = [Inf, 20, 10, 5];
config.seed_list = 11:50;
config.candidate_settings = build_candidate_settings();

fprintf('\n=== Tensor benchmark experiments ===\n');
fprintf('Signal types : %s\n', strjoin(config.signal_types, ', '));
fprintf('SNR values   : %s\n', numvec_to_str(config.snr_list));
fprintf('Seeds        : %d:%d\n', config.seed_list(1), config.seed_list(end));
fprintf('Grouping     : %s\n', config.grouping_mode);
fprintf('Reconstruction: %s\n', config.reconstruction_mode);

trialCsv = fullfile(resultsDir, 'tensor_benchmark_trials.csv');
summaryCsv = fullfile(resultsDir, 'tensor_benchmark_summary.csv');
bestCsv = fullfile(resultsDir, 'tensor_benchmark_best.csv');
matFile = fullfile(resultsDir, 'tensor_benchmark_results.mat');

existing_table = load_existing_trials(trialCsv, config);
num_trials = estimate_num_trials(config);
trial_rows = repmat(init_trial_row(), num_trials, 1);
row = 0;
skipped_existing = 0;

for st = 1:numel(config.signal_types)
    signal_type = string(config.signal_types{st});
    settings = config.candidate_settings.(char(signal_type));
    for si = 1:numel(config.snr_list)
        snr_db = config.snr_list(si);
        for seed = config.seed_list
            [x, components] = generate_signal(signal_type, config.N, snr_db, seed);
            for ci = 1:numel(settings)
                if has_existing_trial(existing_table, signal_type, snr_db, seed, settings(ci))
                    skipped_existing = skipped_existing + 1;
                    continue;
                end
                row = row + 1;
                fprintf('Running %s SNR=%s seed=%d setting=%s\n', ...
                    signal_type, num2str(snr_db), seed, settings(ci).label);
                trial_rows(row) = run_single_trial(signal_type, snr_db, seed, x, components, settings(ci), config);
            end
        end
    end
end

new_table = struct2table(trial_rows(1:row));
trial_table = [existing_table; new_table];
trial_table = restrict_to_current_grid(trial_table, config);
trial_table = sortrows(trial_table, ...
    {'signal_type', 'snr_db', 'seed', 'setting_order'}, ...
    {'ascend', 'descend', 'ascend', 'ascend'});
summary_table = summarize_by_condition(trial_table);
summary_table = sortrows(summary_table, {'signal_type', 'snr_db', 'setting_order'}, {'ascend', 'descend', 'ascend'});
best_table = build_best_table(summary_table);

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
writetable(best_table, bestCsv);
save(matFile, 'config', 'trial_table', 'summary_table', 'best_table');

fprintf('\nSkipped existing trials: %d\n', skipped_existing);
fprintf('New trials run        : %d\n', row);
fprintf('\nCondition summary:\n');
disp(summary_table);
fprintf('\nBest settings by signal/SNR:\n');
disp(best_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', trialCsv, summaryCsv, bestCsv, matFile);

%% --- Local helpers ---
function n = estimate_num_trials(config)
    n = 0;
    for st = 1:numel(config.signal_types)
        signal_type = char(config.signal_types{st});
        n = n + numel(config.snr_list) * numel(config.seed_list) * numel(config.candidate_settings.(signal_type));
    end
end

function trial_table = load_existing_trials(trialCsv, config)
    if ~exist(trialCsv, 'file')
        trial_table = struct2table(repmat(init_trial_row(), 0, 1));
        return;
    end

    trial_table = readtable(trialCsv, 'TextType', 'string');
    trial_table = normalize_trial_table(trial_table, config);
end

function trial_table = normalize_trial_table(trial_table, config)
    string_vars = {'method', 'signal_type', 'setting_label', 'stride_label', 'failure_reason'};
    for i = 1:numel(string_vars)
        name = string_vars{i};
        if ismember(name, trial_table.Properties.VariableNames)
            trial_table.(name) = string(trial_table.(name));
        end
    end
    if ismember('success', trial_table.Properties.VariableNames)
        trial_table.success = logical(trial_table.success);
    end

    for i = 1:height(trial_table)
        setting = setting_for_row(trial_table.signal_type(i), ...
            trial_table.base_window_length(i), trial_table.cpd_rank(i), config);
        if ~isempty(setting)
            trial_table.setting_label(i) = setting.label;
            trial_table.setting_order(i) = setting.order;
            trial_table.stride_label(i) = stride_label_from_values(setting.strides);
        end
    end
end

function setting = setting_for_row(signal_type, base_window_length, cpd_rank, config)
    settings = config.candidate_settings.(char(signal_type));
    setting = [];
    for i = 1:numel(settings)
        if settings(i).base_window_length == base_window_length && settings(i).cpd_rank == cpd_rank
            setting = settings(i);
            return;
        end
    end
end

function yes = has_existing_trial(trial_table, signal_type, snr_db, seed, setting)
    if isempty(trial_table) || height(trial_table) == 0
        yes = false;
        return;
    end
    mask = trial_table.signal_type == signal_type & ...
        snr_equal(trial_table.snr_db, snr_db) & ...
        trial_table.seed == seed & ...
        trial_table.base_window_length == setting.base_window_length & ...
        trial_table.cpd_rank == setting.cpd_rank;
    yes = any(mask);
end

function trial_table = restrict_to_current_grid(trial_table, config)
    keep = false(height(trial_table), 1);
    for i = 1:height(trial_table)
        if ~ismember(char(trial_table.signal_type(i)), config.signal_types)
            continue;
        end
        if ~any(arrayfun(@(x) snr_equal(trial_table.snr_db(i), x), config.snr_list))
            continue;
        end
        if ~ismember(trial_table.seed(i), config.seed_list)
            continue;
        end
        setting = setting_for_row(trial_table.signal_type(i), ...
            trial_table.base_window_length(i), trial_table.cpd_rank(i), config);
        keep(i) = ~isempty(setting);
    end
    trial_table = trial_table(keep, :);
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

function row = run_single_trial(signal_type, snr_db, seed, x, components, setting, config)
    row = init_trial_row();
    row.signal_type = signal_type;
    row.snr_db = snr_db;
    row.seed = seed;
    row.setting_label = setting.label;
    row.setting_order = setting.order;
    row.stride_label = stride_label_from_values(setting.strides);
    row.base_window_length = setting.base_window_length;
    row.cpd_rank = setting.cpd_rank;

    params = tensor_merge_params(struct( ...
        'strides', setting.strides, ...
        'base_window_length', setting.base_window_length, ...
        'cpd_rank', setting.cpd_rank, ...
        'grouping_mode', config.grouping_mode, ...
        'reconstruction_mode', config.reconstruction_mode, ...
        'cpd_method', 'cpd', ...
        'cpd_options', struct()));

    try
        outcome = evaluate_tensor_cpd(signal_type, x, components, params, config.fs);
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
        {'signal_type', 'snr_db', 'setting_label', 'setting_order', 'stride_label', 'base_window_length', 'cpd_rank'}), ...
        'rows');
    rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            snr_equal(trial_table.snr_db, keys.snr_db(i)) & ...
            trial_table.setting_label == keys.setting_label(i);
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
    keys = unique(summary_table(:, {'signal_type', 'snr_db'}), 'rows');
    rows = repmat(init_best_row(), height(keys), 1);
    for i = 1:height(keys)
        subset = summary_table(summary_table.signal_type == keys.signal_type(i) & ...
            snr_equal(summary_table.snr_db, keys.snr_db(i)), :);
        rows(i) = choose_best_row(subset, keys.signal_type(i), keys.snr_db(i));
    end
    best_table = struct2table(rows);
    best_table = sortrows(best_table, {'signal_type', 'snr_db'}, {'ascend', 'descend'});
end

function row = init_best_row()
    row = struct( ...
        'signal_type', "", ...
        'snr_db', nan, ...
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

function row = choose_best_row(subset, signal_type, snr_db)
    row = init_best_row();
    row.signal_type = signal_type;
    row.snr_db = snr_db;
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
