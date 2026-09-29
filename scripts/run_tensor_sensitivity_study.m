%RUN_TENSOR_SENSITIVITY_STUDY Parameter sweep for decimated Hankel tensor CPD.
%
%   Goal:
%   - assess how CP tensor fit changes with CP rank, stride configuration,
%     and base window length
%   - identify whether the low-fit verification result is due to an
%     underpowered setting or a broader modeling limitation
%
%   Writes:
%     results/tensor_sensitivity_trials.csv
%     results/tensor_sensitivity_summary.csv
%     results/tensor_sensitivity_best_by_condition.csv
%     results/tensor_sensitivity_best_by_stride.csv
%     results/tensor_sensitivity_results.mat
%     results/tensor_sensitivity_fit_heatmaps.png
%     results/tensor_sensitivity_best_by_stride.png

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
config.cpd_method = 'cpd';
config.signal_types = {'multiscale', 'close'};
config.snr_list = [Inf, 20, 10, 5];
% Keep the default sweep broad enough to compare settings, but small enough
% to remain practical for CPD experiments.
config.seed_list = 11:13;
config.rank_list = [2, 3, 4, 5, 6];
config.window_list = [200, 400, 600, 1000];
config.stride_setting_labels = {'S1_S2', 'S1_S2_S4', 'S1_S2_S4_S8', 'S1_S4_S8'};
config.stride_setting_values = {
    [1, 2], ...
    [1, 2, 4], ...
    [1, 2, 4, 8], ...
    [1, 4, 8]
    };
config.primary_stride_label = "S1_S2_S4_S8";

fprintf('\n=== Tensor sensitivity study ===\n');
fprintf('Signal types     : %s\n', strjoin(config.signal_types, ', '));
fprintf('SNR values       : %s\n', numvec_to_str(config.snr_list));
fprintf('Seeds            : %d:%d\n', config.seed_list(1), config.seed_list(end));
fprintf('CP ranks         : %s\n', numvec_to_str(config.rank_list));
fprintf('Window lengths   : %s\n', numvec_to_str(config.window_list));
fprintf('Stride settings  : %s\n', strjoin(config.stride_setting_labels, ', '));

num_trials = numel(config.signal_types) * numel(config.snr_list) * numel(config.seed_list) * ...
    numel(config.rank_list) * numel(config.window_list) * numel(config.stride_setting_values);
trial_rows = repmat(init_trial_row(), num_trials, 1);
row = 0;

for st = 1:numel(config.signal_types)
    signal_type = config.signal_types{st};
    for si = 1:numel(config.snr_list)
        snr_db = config.snr_list(si);
        for seed = config.seed_list
            [x, ~] = generate_signal(signal_type, config.N, snr_db, seed);

            for wi = 1:numel(config.window_list)
                base_window_length = config.window_list(wi);
                for ri = 1:numel(config.rank_list)
                    cpd_rank = config.rank_list(ri);
                    for gi = 1:numel(config.stride_setting_values)
                        stride_label = string(config.stride_setting_labels{gi});
                        strides = config.stride_setting_values{gi};
                        row = row + 1;
                        trial_rows(row) = run_single_trial( ...
                            x, signal_type, snr_db, seed, base_window_length, ...
                            cpd_rank, stride_label, strides, config);
                    end
                end
            end
        end
    end
end

trial_table = struct2table(trial_rows(1:row));
summary_table = summarize_by_condition(trial_table);
summary_table = sortrows(summary_table, ...
    {'signal_type', 'snr_db', 'stride_order', 'base_window_length', 'cpd_rank'}, ...
    {'ascend', 'descend', 'ascend', 'ascend', 'ascend'});

best_condition_table = build_best_condition_table(summary_table);
best_stride_table = build_best_stride_table(summary_table);

trialCsv = fullfile(resultsDir, 'tensor_sensitivity_trials.csv');
summaryCsv = fullfile(resultsDir, 'tensor_sensitivity_summary.csv');
bestConditionCsv = fullfile(resultsDir, 'tensor_sensitivity_best_by_condition.csv');
bestStrideCsv = fullfile(resultsDir, 'tensor_sensitivity_best_by_stride.csv');
matFile = fullfile(resultsDir, 'tensor_sensitivity_results.mat');

writetable(trial_table, trialCsv);
writetable(summary_table, summaryCsv);
writetable(best_condition_table, bestConditionCsv);
writetable(best_stride_table, bestStrideCsv);
save(matFile, 'config', 'trial_table', 'summary_table', 'best_condition_table', 'best_stride_table');

plotFiles = plot_tensor_sensitivity(summary_table, best_stride_table, resultsDir, config);

fprintf('\nBest settings by signal/SNR:\n');
disp(best_condition_table);
fprintf('\nBest settings by stride:\n');
disp(best_stride_table);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    trialCsv, summaryCsv, bestConditionCsv, bestStrideCsv, matFile, plotFiles{1}, plotFiles{2});

%% --- Local helpers ---
function row = init_trial_row()
    row = struct( ...
        'method', "TensorCPD", ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'seed', nan, ...
        'stride_label', "", ...
        'stride_order', nan, ...
        'num_strides', nan, ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'tensor_rows', nan, ...
        'tensor_cols', nan, ...
        'tensor_slices', nan, ...
        'tensor_norm', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'fit', nan, ...
        'rel_error', nan, ...
        'lambda_abs_sum', nan, ...
        'lambda_abs_max', nan, ...
        'iterations', nan);
end

function row = run_single_trial(x, signal_type, snr_db, seed, base_window_length, cpd_rank, stride_label, strides, config)
    row = init_trial_row();
    row.signal_type = string(signal_type);
    row.snr_db = snr_db;
    row.seed = seed;
    row.stride_label = stride_label;
    row.stride_order = stride_order_from_label(stride_label);
    row.num_strides = numel(strides);
    row.base_window_length = base_window_length;
    row.cpd_rank = cpd_rank;

    params = tensor_merge_params(struct( ...
        'strides', strides, ...
        'base_window_length', base_window_length, ...
        'cpd_rank', cpd_rank, ...
        'cpd_method', config.cpd_method, ...
        'cpd_options', struct()));

    try
        tensor_out = tensor_build_decimated_hankel(x, params);
        row.tensor_rows = tensor_out.tensor_size(1);
        row.tensor_cols = tensor_out.tensor_size(2);
        row.tensor_slices = tensor_out.tensor_size(3);
        row.tensor_norm = tensor_out.tensor_norm;

        cpd_out = tensor_run_cpd(tensor_out, params);
        row.success = cpd_out.success;
        row.failure_reason = string(cpd_out.failure_reason);
        row.fit = cpd_out.fit;
        row.rel_error = cpd_out.rel_error;

        if ~isempty(cpd_out.lambda)
            row.lambda_abs_sum = sum(abs(cpd_out.lambda));
            row.lambda_abs_max = max(abs(cpd_out.lambda));
        end

        if isfield(cpd_out, 'diagnostics') && isfield(cpd_out.diagnostics, 'iterations')
            row.iterations = cpd_out.diagnostics.iterations;
        end
    catch ME
        row.success = false;
        row.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
    end
end

function summary_table = summarize_by_condition(trial_table)
    keys = unique(trial_table(:, ...
        {'signal_type', 'snr_db', 'stride_label', 'stride_order', 'num_strides', 'base_window_length', 'cpd_rank'}), ...
        'rows');
    rows = repmat(init_summary_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            snr_equal(trial_table.snr_db, keys.snr_db(i)) & ...
            trial_table.stride_label == keys.stride_label(i) & ...
            trial_table.base_window_length == keys.base_window_length(i) & ...
            trial_table.cpd_rank == keys.cpd_rank(i);
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
        'stride_label', "", ...
        'stride_order', nan, ...
        'num_strides', nan, ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_fit', nan, ...
        'median_fit', nan, ...
        'mean_rel_error', nan, ...
        'median_rel_error', nan, ...
        'mean_tensor_rows', nan, ...
        'mean_tensor_cols', nan, ...
        'mean_tensor_slices', nan, ...
        'mean_lambda_abs_sum', nan, ...
        'mean_lambda_abs_max', nan, ...
        'mean_iterations', nan);
end

function row = pack_summary_row(key_row, subset)
    row = init_summary_row();
    row.signal_type = key_row.signal_type;
    row.snr_db = key_row.snr_db;
    row.stride_label = key_row.stride_label;
    row.stride_order = key_row.stride_order;
    row.num_strides = key_row.num_strides;
    row.base_window_length = key_row.base_window_length;
    row.cpd_rank = key_row.cpd_rank;
    row.num_trials = height(subset);
    row.failure_rate = mean(~subset.success);
    row.mean_fit = mean(subset.fit, 'omitnan');
    row.median_fit = median(subset.fit, 'omitnan');
    row.mean_rel_error = mean(subset.rel_error, 'omitnan');
    row.median_rel_error = median(subset.rel_error, 'omitnan');
    row.mean_tensor_rows = mean(subset.tensor_rows, 'omitnan');
    row.mean_tensor_cols = mean(subset.tensor_cols, 'omitnan');
    row.mean_tensor_slices = mean(subset.tensor_slices, 'omitnan');
    row.mean_lambda_abs_sum = mean(subset.lambda_abs_sum, 'omitnan');
    row.mean_lambda_abs_max = mean(subset.lambda_abs_max, 'omitnan');
    row.mean_iterations = mean(subset.iterations, 'omitnan');
end

function best_table = build_best_condition_table(summary_table)
    keys = unique(summary_table(:, {'signal_type', 'snr_db'}), 'rows');
    rows = repmat(init_best_row(), height(keys), 1);
    for i = 1:height(keys)
        subset = summary_table(summary_table.signal_type == keys.signal_type(i) & ...
            snr_equal(summary_table.snr_db, keys.snr_db(i)), :);
        rows(i) = choose_best_row(subset, keys.signal_type(i), keys.snr_db(i), "");
    end
    best_table = struct2table(rows);
    best_table = sortrows(best_table, {'signal_type', 'snr_db'}, {'ascend', 'descend'});
end

function best_table = build_best_stride_table(summary_table)
    keys = unique(summary_table(:, {'signal_type', 'snr_db', 'stride_label', 'stride_order'}), 'rows');
    rows = repmat(init_best_row(), height(keys), 1);
    for i = 1:height(keys)
        subset = summary_table(summary_table.signal_type == keys.signal_type(i) & ...
            snr_equal(summary_table.snr_db, keys.snr_db(i)) & ...
            summary_table.stride_label == keys.stride_label(i), :);
        rows(i) = choose_best_row(subset, keys.signal_type(i), keys.snr_db(i), keys.stride_label(i));
    end
    best_table = struct2table(rows);
    best_table = sortrows(best_table, {'signal_type', 'snr_db', 'stride_order'}, {'ascend', 'descend', 'ascend'});
end

function row = init_best_row()
    row = struct( ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'stride_label', "", ...
        'stride_order', nan, ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'failure_rate', nan, ...
        'mean_fit', nan, ...
        'median_fit', nan, ...
        'mean_rel_error', nan, ...
        'mean_tensor_rows', nan, ...
        'mean_tensor_cols', nan, ...
        'mean_tensor_slices', nan);
end

function row = choose_best_row(subset, signal_type, snr_db, stride_label)
    row = init_best_row();
    row.signal_type = signal_type;
    row.snr_db = snr_db;
    row.stride_label = stride_label;

    if isempty(subset)
        return;
    end

    if any(~isnan(subset.mean_fit))
        sort_subset = sortrows(subset, {'mean_fit', 'failure_rate', 'mean_rel_error'}, ...
            {'descend', 'ascend', 'ascend'});
        best = sort_subset(1, :);
        row.stride_label = best.stride_label;
        row.stride_order = best.stride_order;
        row.base_window_length = best.base_window_length;
        row.cpd_rank = best.cpd_rank;
        row.failure_rate = best.failure_rate;
        row.mean_fit = best.mean_fit;
        row.median_fit = best.median_fit;
        row.mean_rel_error = best.mean_rel_error;
        row.mean_tensor_rows = best.mean_tensor_rows;
        row.mean_tensor_cols = best.mean_tensor_cols;
        row.mean_tensor_slices = best.mean_tensor_slices;
    end
end

function plotFiles = plot_tensor_sensitivity(summary_table, best_stride_table, resultsDir, config)
    signal_types = string(config.signal_types);
    snr_vals = sorted_snr(config.snr_list);
    rank_vals = config.rank_list;
    window_vals = config.window_list;

    heatmapFig = figure('Name', 'Tensor sensitivity: fit heatmaps', 'NumberTitle', 'off');
    tiledlayout(numel(signal_types), numel(snr_vals), 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signal_types)
        for ni = 1:numel(snr_vals)
            ax = nexttile;
            subset = summary_table(summary_table.signal_type == signal_types(si) & ...
                snr_equal(summary_table.snr_db, snr_vals(ni)) & ...
                summary_table.stride_label == config.primary_stride_label, :);
            M = nan(numel(window_vals), numel(rank_vals));
            for wi = 1:numel(window_vals)
                for ri = 1:numel(rank_vals)
                    mask = subset.base_window_length == window_vals(wi) & subset.cpd_rank == rank_vals(ri);
                    if any(mask)
                        M(wi, ri) = subset.mean_fit(find(mask, 1));
                    end
                end
            end
            imagesc(ax, M, [0, 1]);
            xticks(ax, 1:numel(rank_vals));
            xticklabels(ax, arrayfun(@num2str, rank_vals, 'UniformOutput', false));
            yticks(ax, 1:numel(window_vals));
            yticklabels(ax, arrayfun(@num2str, window_vals, 'UniformOutput', false));
            xlabel(ax, 'CP rank');
            ylabel(ax, 'Base window length');
            title(ax, sprintf('%s, SNR=%s', signal_types(si), snr_label_scalar(snr_vals(ni))));
            colorbar(ax);
        end
    end
    sgtitle(sprintf('Mean CP fit for primary stride setting (%s)', config.primary_stride_label));
    heatmapPng = fullfile(resultsDir, 'tensor_sensitivity_fit_heatmaps.png');
    exportgraphics(heatmapFig, heatmapPng, 'Resolution', 300);

    strideFig = figure('Name', 'Tensor sensitivity: best fit by stride', 'NumberTitle', 'off');
    tiledlayout(numel(signal_types), numel(snr_vals), 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signal_types)
        for ni = 1:numel(snr_vals)
            ax = nexttile;
            subset = best_stride_table(best_stride_table.signal_type == signal_types(si) & ...
                snr_equal(best_stride_table.snr_db, snr_vals(ni)), :);
            subset = sortrows(subset, 'stride_order');
            bar(ax, subset.mean_fit);
            ylim(ax, [0, 1]);
            xticks(ax, 1:height(subset));
            xticklabels(ax, cellstr(subset.stride_label));
            xtickangle(ax, 25);
            ylabel(ax, 'Best mean fit');
            xlabel(ax, 'Stride setting');
            title(ax, sprintf('%s, SNR=%s', signal_types(si), snr_label_scalar(snr_vals(ni))));
            grid(ax, 'on');
        end
    end
    sgtitle('Best mean fit per stride setting (best rank/window chosen within each stride set)');
    stridePng = fullfile(resultsDir, 'tensor_sensitivity_best_by_stride.png');
    exportgraphics(strideFig, stridePng, 'Resolution', 300);

    plotFiles = {heatmapPng, stridePng};
end

function order = stride_order_from_label(label)
    switch string(label)
        case "S1_S2"
            order = 1;
        case "S1_S2_S4"
            order = 2;
        case "S1_S2_S4_S8"
            order = 3;
        case "S1_S4_S8"
            order = 4;
        otherwise
            error('run_tensor_sensitivity_study:UnknownStrideLabel', ...
                'Unknown stride label %s', label);
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
