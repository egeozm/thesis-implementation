%RUN_SSA_CLOSE_TRUTH_GUIDED_DIAGNOSTIC
% Compare automatic vs truth-guided grouping on the close-frequency noiseless case.
%
%   This diagnostic isolates whether failure comes from the automatic grouping rule
%   or from SSA separability itself.

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));

resultsDir = fullfile(projRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

calibrationMat = fullfile(resultsDir, 'ssa_calibration_results.mat');
if ~exist(calibrationMat, 'file')
    error(['Missing calibration results. Run scripts/run_ssa_calibration.m first ', ...
        'to select and freeze the automatic grouping parameters.']);
end

loaded = load(calibrationMat, 'best_params');
best_params = loaded.best_params;

config = struct();
config.signal_type = 'close';
config.N = 2000;
config.snr_db = Inf;
config.seed = 42;
config.fs = 1;
config.L_list = [200, 400, 600, 1000];
config.inspect_r = max(12, best_params.r_max);

fprintf('\n=== Close-frequency noiseless diagnostic ===\n');
fprintf('Frozen automatic params: r_max=%d, sv_ratio_thresh=%.3f, freq_tol=%.4f, min_osc_freq=%.4f\n', ...
    best_params.r_max, best_params.sv_ratio_thresh, best_params.freq_tol, best_params.min_osc_freq);

[x, components, t] = generate_signal(config.signal_type, config.N, config.snr_db, config.seed);
true_components = [components.comp1, components.comp2];

rows = repmat(init_row(), numel(config.L_list), 1);

for i = 1:numel(config.L_list)
    L = config.L_list(i);
    features = ssa_precompute(x, L, config.inspect_r, config.fs);

    automatic = evaluate_auto_group_ssa(config.signal_type, features, components, best_params, config.fs);
    truth_guided = evaluate_close_truth_guided(features, components, config.fs);

    rows(i) = build_row(L, automatic, truth_guided);
end

summary_table = struct2table(rows);
csvFile = fullfile(resultsDir, 'ssa_close_truth_guided_diagnostic.csv');
matFile = fullfile(resultsDir, 'ssa_close_truth_guided_diagnostic.mat');
writetable(summary_table, csvFile);
save(matFile, 'config', 'best_params', 'summary_table');

fprintf('\nDiagnostic summary:\n');
disp(summary_table);
fprintf('\nSaved:\n  %s\n  %s\n', csvFile, matFile);

figure('Name', 'Close-frequency noiseless: automatic vs truth-guided NMSE', 'NumberTitle', 'off');
plot(summary_table.L, summary_table.auto_mean_nmse, 'o-', 'LineWidth', 1.2, 'DisplayName', 'automatic'); hold on;
plot(summary_table.L, summary_table.truth_mean_nmse, 's-', 'LineWidth', 1.2, 'DisplayName', 'truth-guided');
set(gca, 'YScale', 'log');
xlabel('L');
ylabel('Mean NMSE');
title('Close-frequency noiseless diagnostic');
legend('Location', 'best');
grid on;

bestL = summary_table.L(end);
features = ssa_precompute(x, bestL, config.inspect_r, config.fs);
truth_guided = evaluate_close_truth_guided(features, components, config.fs);

figure('Name', 'Close-frequency truth-guided reconstructions', 'NumberTitle', 'off');
sgtitle(sprintf('Close-frequency noiseless, L = %d (truth-guided)', bestL));
for k = 1:2
    subplot(1, 2, k);
    plot(t, true_components(:, k), 'k-', 'LineWidth', 1); hold on;
    plot(t, truth_guided.est_components(:, k), 'r--', 'LineWidth', 1);
    title(sprintf('Component %d', k));
    xlabel('n');
    ylabel('amplitude');
    legend('true', 'truth-guided SSA', 'Location', 'best');
    grid on;
end

%% --- Local helpers ---
function row = init_row()
    row = struct( ...
        'L', nan, ...
        'auto_success', false, ...
        'truth_success', false, ...
        'auto_mean_nmse', nan, ...
        'truth_mean_nmse', nan, ...
        'auto_mean_abs_rho', nan, ...
        'truth_mean_abs_rho', nan, ...
        'auto_pair_1', "", ...
        'auto_pair_2', "", ...
        'truth_pair_1', "", ...
        'truth_pair_2', "");
end

function row = build_row(L, automatic, truth_guided)
    row = init_row();
    row.L = L;
    row.auto_success = automatic.success;
    row.truth_success = truth_guided.success;
    row.auto_mean_nmse = automatic.mean_nmse;
    row.truth_mean_nmse = truth_guided.mean_nmse;
    row.auto_mean_abs_rho = automatic.mean_abs_rho;
    row.truth_mean_abs_rho = truth_guided.mean_abs_rho;
    row.auto_pair_1 = pair_to_string(automatic.grouping.oscillatory_pairs, 1);
    row.auto_pair_2 = pair_to_string(automatic.grouping.oscillatory_pairs, 2);
    row.truth_pair_1 = pair_to_string(truth_guided.groups, 1);
    row.truth_pair_2 = pair_to_string(truth_guided.groups, 2);
end

function outcome = evaluate_close_truth_guided(features, components, fs)
    true_components = [components.comp1, components.comp2];
    inspect_r = min(features.r_max, size(features.elementary, 2));
    candidate_pairs = adjacent_pairs(inspect_r);

    best_score = inf;
    best_groups = {};
    best_est = [];
    best_metrics = empty_metrics(2);

    for i = 1:numel(candidate_pairs)
        for j = (i + 1):numel(candidate_pairs)
            pair1 = candidate_pairs{i};
            pair2 = candidate_pairs{j};
            if any(ismember(pair1, pair2))
                continue;
            end

            est12 = ssa_reconstruct_groups(features, {pair1, pair2});
            m12 = compute_metrics(true_components, est12, [1; 2], fs);
            score12 = mean(m12.nmse, 'omitnan');

            est21 = ssa_reconstruct_groups(features, {pair2, pair1});
            m21 = compute_metrics(true_components, est21, [1; 2], fs);
            score21 = mean(m21.nmse, 'omitnan');

            if score12 <= score21
                score = score12;
                est = est12;
                metrics = m12;
                groups = {pair1, pair2};
            else
                score = score21;
                est = est21;
                metrics = m21;
                groups = {pair2, pair1};
            end

            if score < best_score
                best_score = score;
                best_groups = groups;
                best_est = est;
                best_metrics = metrics;
            end
        end
    end

    outcome = struct();
    outcome.success = ~isempty(best_groups);
    outcome.groups = best_groups;
    outcome.est_components = best_est;
    outcome.metrics = best_metrics;
    outcome.mean_nmse = mean(best_metrics.nmse, 'omitnan');
    outcome.mean_abs_rho = mean(abs(best_metrics.pearson), 'omitnan');
end

function pairs = adjacent_pairs(r_use)
    pairs = cell(1, max(r_use - 1, 0));
    for r = 1:(r_use - 1)
        pairs{r} = [r, r + 1];
    end
end

function metrics = empty_metrics(num_components)
    metrics = struct();
    metrics.nmse = nan(num_components, 1);
    metrics.pearson = nan(num_components, 1);
    metrics.peak_freq_error = nan(num_components, 1);
    metrics.leakage = nan(num_components, 1);
    metrics.dominant_freq_true = nan(num_components, 1);
    metrics.dominant_freq_est = nan(num_components, 1);
    metrics.true_idx = (1:num_components).';
    metrics.est_idx = nan(num_components, 1);
end

function s = pair_to_string(pairs, idx)
    if isempty(pairs) || numel(pairs) < idx || isempty(pairs{idx})
        s = "";
    else
        s = string(mat2str(pairs{idx}));
    end
end
