%RUN_SSA_CALIBRATION Calibrate fixed automatic SSA grouping parameters.
%
%   Uses only the two thesis synthetic signals and saves results under results/.

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));

resultsDir = fullfile(projRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Calibration design
config = struct();
config.N = 2000;
config.fs = 1;
config.signal_types = {'multiscale', 'close'};
config.snr_list = [Inf, 20, 10, 5];
config.seed_list = 1:10;
config.L_list = [200, 400, 600, 1000];
config.r_max_list = [8, 10, 12];
config.sv_ratio_thresh_list = [0.5, 0.6, 0.7, 0.8];
config.freq_tol_list = [0.005, 0.01, 0.02];
config.min_osc_freq_list = [0.005, 0.01, 0.015];
config.selection_margin = 0.05;
config.failure_penalty = 1e3;

param_grid = build_param_grid(config.r_max_list, ...
    config.sv_ratio_thresh_list, config.freq_tol_list, config.min_osc_freq_list);
max_rmax = max(config.r_max_list);

fprintf('\n=== SSA automatic grouping calibration ===\n');
fprintf('Signal types : %s\n', strjoin(config.signal_types, ', '));
fprintf('SNR values   : %s\n', numvec_to_str(config.snr_list));
fprintf('L values     : %s\n', numvec_to_str(config.L_list));
fprintf('Seeds        : %d:%d\n', config.seed_list(1), config.seed_list(end));
fprintf('Grid size    : %d configurations\n', numel(param_grid));

num_trials = numel(config.signal_types) * numel(config.snr_list) * ...
    numel(config.seed_list) * numel(config.L_list) * numel(param_grid);
trial_rows = repmat(init_trial_row(), num_trials, 1);
row = 0;

for st = 1:numel(config.signal_types)
    signal_type = config.signal_types{st};
    for si = 1:numel(config.snr_list)
        snr_db = config.snr_list(si);
        for seed = config.seed_list
            [x, components] = generate_signal(signal_type, config.N, snr_db, seed);
            for Li = 1:numel(config.L_list)
                L = config.L_list(Li);
                features = ssa_precompute(x, L, max_rmax, config.fs);

                for pi = 1:numel(param_grid)
                    outcome = evaluate_auto_group_ssa(signal_type, features, components, ...
                        param_grid(pi), config.fs);
                    row = row + 1;
                    trial_rows(row) = build_trial_row(signal_type, snr_db, seed, L, ...
                        pi, param_grid(pi), outcome, config.failure_penalty);
                end
            end
        end
    end
end

trial_table = struct2table(trial_rows(1:row));
summary_table = summarize_trials(trial_table, param_grid);
summary_table = sortrows(summary_table, ...
    {'score', 'failure_rate', 'mean_nmse', 'mean_abs_rho'}, ...
    {'ascend', 'ascend', 'ascend', 'descend'});

[best_params, best_row, eligible_rows] = select_ssa_params(summary_table, config.selection_margin);

trialCsv = fullfile(resultsDir, 'ssa_calibration_trials.csv');
summaryCsv = fullfile(resultsDir, 'ssa_calibration_summary.csv');
matFile = fullfile(resultsDir, 'ssa_calibration_results.mat');

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
save(matFile, 'config', 'param_grid', 'trial_table', 'summary_table', ...
    'best_params', 'best_row', 'eligible_rows');

fprintf('\nTop calibration settings:\n');
disp(summary_table(1:min(10, height(summary_table)), :));

fprintf('\nSelected fixed parameters (within %.0f%% of best score):\n', ...
    100 * config.selection_margin);
fprintf('  r_max = %d\n', best_params.r_max);
fprintf('  sv_ratio_thresh = %.3f\n', best_params.sv_ratio_thresh);
fprintf('  freq_tol = %.4f\n', best_params.freq_tol);
fprintf('  min_osc_freq = %.4f\n', best_params.min_osc_freq);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n', trialCsv, summaryCsv, matFile);

%% --- Local helpers ---
function grid = build_param_grid(r_max_list, sv_ratio_list, freq_tol_list, min_osc_freq_list)
    idx = 0;
    grid(numel(r_max_list) * numel(sv_ratio_list) * numel(freq_tol_list) * numel(min_osc_freq_list)) = struct();
    for r_max = r_max_list
        for sv_ratio = sv_ratio_list
            for freq_tol = freq_tol_list
                for min_osc_freq = min_osc_freq_list
                    idx = idx + 1;
                    grid(idx).r_max = r_max;
                    grid(idx).sv_ratio_thresh = sv_ratio;
                    grid(idx).freq_tol = freq_tol;
                    grid(idx).min_osc_freq = min_osc_freq;
                end
            end
        end
    end
end

function row = init_trial_row()
    row = struct( ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'seed', nan, ...
        'L', nan, ...
        'param_id', nan, ...
        'r_max', nan, ...
        'sv_ratio_thresh', nan, ...
        'freq_tol', nan, ...
        'min_osc_freq', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'num_pairs', nan, ...
        'trend_group_size', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'score', nan, ...
        'comp1_nmse', nan, ...
        'comp2_nmse', nan, ...
        'comp3_nmse', nan, ...
        'comp1_rho', nan, ...
        'comp2_rho', nan, ...
        'comp3_rho', nan);
end

function row = build_trial_row(signal_type, snr_db, seed, L, param_id, params, outcome, failure_penalty)
    row = init_trial_row();
    row.signal_type = string(signal_type);
    row.snr_db = snr_db;
    row.seed = seed;
    row.L = L;
    row.param_id = param_id;
    row.r_max = params.r_max;
    row.sv_ratio_thresh = params.sv_ratio_thresh;
    row.freq_tol = params.freq_tol;
    row.min_osc_freq = params.min_osc_freq;
    row.success = outcome.success;
    row.failure_reason = outcome.failure_reason;
    row.num_pairs = outcome.num_pairs;
    row.trend_group_size = outcome.trend_group_size;
    row.mean_nmse = outcome.mean_nmse;
    row.mean_abs_rho = outcome.mean_abs_rho;
    row.mean_peak_error = outcome.mean_peak_error;
    row.score = outcome.score;
    if ~row.success
        row.score = failure_penalty;
    end

    num_components = numel(outcome.component_names);
    for k = 1:min(3, num_components)
        row.(sprintf('comp%d_nmse', k)) = outcome.metrics.nmse(k);
        row.(sprintf('comp%d_rho', k)) = outcome.metrics.pearson(k);
    end
end

function summary_table = summarize_trials(trial_table, param_grid)
    summary_rows = repmat(init_summary_row(), numel(param_grid), 1);
    for pi = 1:numel(param_grid)
        mask = trial_table.param_id == pi;
        subset = trial_table(mask, :);
        summary_rows(pi) = init_summary_row();
        summary_rows(pi).param_id = pi;
        summary_rows(pi).r_max = param_grid(pi).r_max;
        summary_rows(pi).sv_ratio_thresh = param_grid(pi).sv_ratio_thresh;
        summary_rows(pi).freq_tol = param_grid(pi).freq_tol;
        summary_rows(pi).min_osc_freq = param_grid(pi).min_osc_freq;
        summary_rows(pi).num_trials = height(subset);
        summary_rows(pi).failure_rate = mean(~subset.success);
        summary_rows(pi).score = mean(subset.score, 'omitnan');
        summary_rows(pi).mean_nmse = mean(subset.mean_nmse, 'omitnan');
        summary_rows(pi).mean_abs_rho = mean(subset.mean_abs_rho, 'omitnan');
        summary_rows(pi).mean_peak_error = mean(subset.mean_peak_error, 'omitnan');
        summary_rows(pi).mean_pairs = mean(subset.num_pairs, 'omitnan');
        summary_rows(pi).mean_trend_group_size = mean(subset.trend_group_size, 'omitnan');
    end
    summary_table = struct2table(summary_rows);
end

function row = init_summary_row()
    row = struct( ...
        'param_id', nan, ...
        'r_max', nan, ...
        'sv_ratio_thresh', nan, ...
        'freq_tol', nan, ...
        'min_osc_freq', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'score', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_pairs', nan, ...
        'mean_trend_group_size', nan);
end

function out = numvec_to_str(v)
    parts = strings(size(v));
    for i = 1:numel(v)
        if isinf(v(i))
            parts(i) = "Inf";
        else
            parts(i) = string(v(i));
        end
    end
    out = strjoin(cellstr(parts), ', ');
end
