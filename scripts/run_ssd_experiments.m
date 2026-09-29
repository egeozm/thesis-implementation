%RUN_SSD_EXPERIMENTS Run SSD experiments with frozen SSA grouping hyperparameters.
%
%   Requires results/ssa_calibration_results.mat (best_params) from run_ssa_calibration.m.
%   Window length is adaptive (candidate set in params.L_list); trial rows use L = NaN.

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

ssd_params = ssd_merge_params(best_params);

fprintf('\n=== SSD experiments (adaptive L, frozen grouping rule) ===\n');
fprintf(['Using frozen grouping params: r_max=%d, sv_ratio_thresh=%.3f, ', ...
    'freq_tol=%.4f, min_osc_freq=%.4f\n'], ...
    ssd_params.r_max, ssd_params.sv_ratio_thresh, ...
    ssd_params.freq_tol, ssd_params.min_osc_freq);
fprintf('L candidates : %s\n', numvec_to_str(ssd_params.L_list));

num_trials = numel(config.signal_types) * numel(config.snr_list) * numel(config.seed_list);
trial_rows = repmat(init_trial_row(), num_trials, 1);
row = 0;

for st = 1:numel(config.signal_types)
    signal_type = config.signal_types{st};
    for si = 1:numel(config.snr_list)
        snr_db = config.snr_list(si);
        for seed = config.seed_list
            [x, components] = generate_signal(signal_type, config.N, snr_db, seed);
            outcome = evaluate_ssd(signal_type, x, components, ssd_params, config.fs);
            row = row + 1;
            trial_rows(row) = build_trial_row(signal_type, snr_db, seed, outcome);
        end
    end
end

trial_table = struct2table(trial_rows(1:row));
summary_table = summarize_by_condition(trial_table);
summary_table = sortrows(summary_table, {'signal_type', 'snr_db'}, {'ascend', 'descend'});

overall_table = summarize_overall(trial_table);

trialCsv = fullfile(resultsDir, 'ssd_trials.csv');
summaryCsv = fullfile(resultsDir, 'ssd_summary.csv');
overallCsv = fullfile(resultsDir, 'ssd_overall.csv');
matFile = fullfile(resultsDir, 'ssd_results.mat');

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
writetable(overall_table, overallCsv);
save(matFile, 'config', 'ssd_params', 'trial_table', 'summary_table', 'overall_table');

fprintf('\nCondition summary:\n');
disp(summary_table);
fprintf('\nOverall summary:\n');
disp(overall_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', trialCsv, summaryCsv, overallCsv, matFile);

%% --- Local helpers ---
function row = init_trial_row()
    row = struct( ...
        'method', "SSD", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'seed', nan, ...
        'L', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'num_pairs', nan, ...
        'trend_group_size', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'comp1_nmse', nan, ...
        'comp2_nmse', nan, ...
        'comp3_nmse', nan, ...
        'comp1_rho', nan, ...
        'comp2_rho', nan, ...
        'comp3_rho', nan);
end

function row = build_trial_row(signal_type, snr_db, seed, outcome)
    row = init_trial_row();
    row.signal_type = string(signal_type);
    row.snr_db = snr_db;
    row.seed = seed;
    row.L = nan;
    row.success = outcome.success;
    row.failure_reason = outcome.failure_reason;
    row.num_pairs = outcome.num_ssd_modes;
    row.trend_group_size = nan;
    row.mean_nmse = outcome.mean_nmse;
    row.mean_abs_rho = outcome.mean_abs_rho;
    row.mean_peak_error = outcome.mean_peak_error;

    num_components = numel(outcome.component_names);
    for k = 1:min(3, num_components)
        row.(sprintf('comp%d_nmse', k)) = outcome.metrics.nmse(k);
        row.(sprintf('comp%d_rho', k)) = outcome.metrics.pearson(k);
    end
end

function summary_table = summarize_by_condition(trial_table)
    keys = unique(trial_table(:, {'signal_type', 'snr_db'}), 'rows');
    summary_rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            trial_table.snr_db == keys.snr_db(i);
        subset = trial_table(mask, :);
        summary_rows(i) = aggregate_subset(keys.signal_type(i), keys.snr_db(i), subset);
    end
    summary_table = struct2table(summary_rows);
end

function overall_table = summarize_overall(trial_table)
    signal_types = unique(trial_table.signal_type);
    overall_rows = repmat(init_summary_row(), numel(signal_types), 1);
    for i = 1:numel(signal_types)
        mask = trial_table.signal_type == signal_types(i);
        subset = trial_table(mask, :);
        overall_rows(i) = aggregate_subset(signal_types(i), nan, subset);
    end
    overall_table = struct2table(overall_rows);
end

function row = aggregate_subset(signal_type, snr_db, subset)
    row = init_summary_row();
    row.signal_type = signal_type;
    row.snr_db = snr_db;
    row.num_trials = height(subset);
    row.failure_rate = mean(~subset.success);
    row.mean_nmse = mean(subset.mean_nmse, 'omitnan');
    row.mean_abs_rho = mean(subset.mean_abs_rho, 'omitnan');
    row.mean_peak_error = mean(subset.mean_peak_error, 'omitnan');
    row.mean_pairs = mean(subset.num_pairs, 'omitnan');
    row.mean_trend_group_size = mean(subset.trend_group_size, 'omitnan');
end

function row = init_summary_row()
    row = struct( ...
        'method', "SSD", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_pairs', nan, ...
        'mean_trend_group_size', nan);
end

function s = numvec_to_str(v)
    s = strjoin(arrayfun(@num2str, v, 'UniformOutput', false), ', ');
end
