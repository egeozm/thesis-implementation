%RUN_COMPARE_SSA_SSD_TENSOR Aggregate SSA, SSD, and tensor benchmark results.
%
%   Requires:
%     results/ssa_frozen_results.mat
%     results/ssd_results.mat
%     results/tensor_benchmark_results.mat
%
%   Writes:
%     results/ssa_ssd_tensor_comparison_long.csv
%     results/ssa_ssd_tensor_comparison_by_condition.csv

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
resultsDir = fullfile(projRoot, 'results');

ssaMat = fullfile(resultsDir, 'ssa_frozen_results.mat');
ssdMat = fullfile(resultsDir, 'ssd_results.mat');
tensorMat = fullfile(resultsDir, 'tensor_benchmark_results.mat');

require_file(ssaMat, 'scripts/run_ssa_frozen_experiments.m');
require_file(ssdMat, 'scripts/run_ssd_experiments.m');
require_file(tensorMat, 'scripts/run_tensor_benchmark_experiments.m');

S = load(ssaMat, 'trial_table');
D = load(ssdMat, 'trial_table');
T = load(tensorMat, 'best_table');

ssa_agg = aggregate_ssa_or_ssd(S.trial_table, "SSA");
ssd_agg = aggregate_ssa_or_ssd(D.trial_table, "SSD");
tensor_agg = aggregate_tensor_best(T.best_table);

long_table = [ssa_agg; ssd_agg; tensor_agg];
long_table = sortrows(long_table, {'signal_type', 'snr_db', 'method'}, {'ascend', 'descend', 'ascend'});
longCsv = fullfile(resultsDir, 'ssa_ssd_tensor_comparison_long.csv');
writetable(long_table, longCsv);

wide_table = build_wide_table(ssa_agg, ssd_agg, tensor_agg);
wideCsv = fullfile(resultsDir, 'ssa_ssd_tensor_comparison_by_condition.csv');
writetable(wide_table, wideCsv);

fprintf('\n=== SSA vs SSD vs Tensor benchmark ===\n');
disp(wide_table);
fprintf('\nSaved:\n  %s\n  %s\n', wideCsv, longCsv);

%% --- Local helpers ---
function require_file(path_str, script_name)
    if ~exist(path_str, 'file')
        error('Missing %s. Run %s first.', path_str, script_name);
    end
end

function agg = aggregate_ssa_or_ssd(trial_table, method_name)
    keys = unique(trial_table(:, {'signal_type', 'snr_db'}), 'rows');
    rows = repmat(init_agg_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            trial_table.snr_db == keys.snr_db(i);
        sub = trial_table(mask, :);
        rows(i) = pack_agg_row(method_name, keys.signal_type(i), keys.snr_db(i), sub, mean(sub.num_pairs, 'omitnan'));
    end
    agg = struct2table(rows);
end

function agg = aggregate_tensor_best(best_table)
    rows = repmat(init_agg_row(), height(best_table), 1);
    for i = 1:height(best_table)
        extra = best_table.mean_num_est_components(i);
        rows(i) = pack_direct_row("TensorCPD", best_table.signal_type(i), best_table.snr_db(i), ...
            1, best_table.failure_rate(i), best_table.mean_nmse(i), best_table.mean_abs_rho(i), ...
            best_table.mean_peak_error(i), extra);
    end
    agg = struct2table(rows);
end

function r = init_agg_row()
    r = struct( ...
        'method', "", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_pairs', nan);
end

function r = pack_agg_row(method_name, signal_type, snr_db, sub, extra_value)
    r = init_agg_row();
    r.method = method_name;
    r.signal_type = signal_type;
    r.snr_db = snr_db;
    r.num_trials = height(sub);
    r.failure_rate = mean(~sub.success);
    r.mean_nmse = mean(sub.mean_nmse, 'omitnan');
    r.mean_abs_rho = mean(sub.mean_abs_rho, 'omitnan');
    r.mean_peak_error = mean(sub.mean_peak_error, 'omitnan');
    r.mean_pairs = extra_value;
end

function r = pack_direct_row(method_name, signal_type, snr_db, num_trials, failure_rate, mean_nmse, mean_abs_rho, mean_peak_error, extra_value)
    r = init_agg_row();
    r.method = method_name;
    r.signal_type = signal_type;
    r.snr_db = snr_db;
    r.num_trials = num_trials;
    r.failure_rate = failure_rate;
    r.mean_nmse = mean_nmse;
    r.mean_abs_rho = mean_abs_rho;
    r.mean_peak_error = mean_peak_error;
    r.mean_pairs = extra_value;
end

function wide = build_wide_table(ssa_agg, ssd_agg, tensor_agg)
    keys = unique([ ...
        ssa_agg(:, {'signal_type', 'snr_db'}); ...
        ssd_agg(:, {'signal_type', 'snr_db'}); ...
        tensor_agg(:, {'signal_type', 'snr_db'})], 'rows');

    wide = table();
    wide.signal_type = keys.signal_type;
    wide.snr_db = keys.snr_db;
    n = height(keys);
    nan_col = nan(n, 1);

    wide.ssa_mean_nmse = nan_col;
    wide.ssa_mean_abs_rho = nan_col;
    wide.ssa_failure_rate = nan_col;

    wide.ssd_mean_nmse = nan_col;
    wide.ssd_mean_abs_rho = nan_col;
    wide.ssd_failure_rate = nan_col;

    wide.tensor_mean_nmse = nan_col;
    wide.tensor_mean_abs_rho = nan_col;
    wide.tensor_failure_rate = nan_col;
    wide.tensor_mean_pairs = nan_col;

    wide.nmse_delta_tensor_minus_ssa = nan_col;
    wide.nmse_delta_tensor_minus_ssd = nan_col;

    for i = 1:n
        st = keys.signal_type(i);
        sn = keys.snr_db(i);
        ia = find(ssa_agg.signal_type == st & ssa_agg.snr_db == sn, 1);
        ib = find(ssd_agg.signal_type == st & ssd_agg.snr_db == sn, 1);
        ic = find(tensor_agg.signal_type == st & tensor_agg.snr_db == sn, 1);

        if ~isempty(ia)
            wide.ssa_mean_nmse(i) = ssa_agg.mean_nmse(ia);
            wide.ssa_mean_abs_rho(i) = ssa_agg.mean_abs_rho(ia);
            wide.ssa_failure_rate(i) = ssa_agg.failure_rate(ia);
        end
        if ~isempty(ib)
            wide.ssd_mean_nmse(i) = ssd_agg.mean_nmse(ib);
            wide.ssd_mean_abs_rho(i) = ssd_agg.mean_abs_rho(ib);
            wide.ssd_failure_rate(i) = ssd_agg.failure_rate(ib);
        end
        if ~isempty(ic)
            wide.tensor_mean_nmse(i) = tensor_agg.mean_nmse(ic);
            wide.tensor_mean_abs_rho(i) = tensor_agg.mean_abs_rho(ic);
            wide.tensor_failure_rate(i) = tensor_agg.failure_rate(ic);
            wide.tensor_mean_pairs(i) = tensor_agg.mean_pairs(ic);
        end
        if ~isempty(ia) && ~isempty(ic)
            wide.nmse_delta_tensor_minus_ssa(i) = tensor_agg.mean_nmse(ic) - ssa_agg.mean_nmse(ia);
        end
        if ~isempty(ib) && ~isempty(ic)
            wide.nmse_delta_tensor_minus_ssd(i) = tensor_agg.mean_nmse(ic) - ssd_agg.mean_nmse(ib);
        end
    end
end
