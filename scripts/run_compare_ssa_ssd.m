%RUN_COMPARE_SSA_SSD Aggregate frozen SSA trials and SSD trials side-by-side.
%
%   Requires:
%     results/ssa_frozen_results.mat  (from run_ssa_frozen_experiments.m)
%     results/ssd_results.mat       (from run_ssd_experiments.m)
%
%   Writes results/ssa_ssd_comparison_by_condition.csv (wide) and
%   results/ssa_ssd_comparison_long.csv (long: one row per method and condition),
%   plus:
%     results/ssa_ssd_mean_nmse_vs_snr.png
%     results/ssa_ssd_mean_abs_rho_vs_snr.png
%     results/ssa_ssd_nmse_delta_vs_snr.png

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
resultsDir = fullfile(projRoot, 'results');

ssaMat = fullfile(resultsDir, 'ssa_frozen_results.mat');
ssdMat = fullfile(resultsDir, 'ssd_results.mat');
if ~exist(ssaMat, 'file')
    error('Missing %s. Run scripts/run_ssa_frozen_experiments.m first.', ssaMat);
end
if ~exist(ssdMat, 'file')
    error('Missing %s. Run scripts/run_ssd_experiments.m first.', ssdMat);
end

S = load(ssaMat, 'trial_table');
D = load(ssdMat, 'trial_table');
ssa_t = S.trial_table;
ssd_t = D.trial_table;

ssa_agg = aggregate_by_signal_snr(ssa_t);
ssd_agg = aggregate_by_signal_snr(ssd_t);

ssa_agg.method(:) = repmat("SSA", height(ssa_agg), 1);
ssd_agg.method(:) = repmat("SSD", height(ssd_agg), 1);

long_table = [ssa_agg; ssd_agg];
longCsv = fullfile(resultsDir, 'ssa_ssd_comparison_long.csv');
writetable(long_table, longCsv);

wide_table = build_wide_table(ssa_agg, ssd_agg);
wideCsv = fullfile(resultsDir, 'ssa_ssd_comparison_by_condition.csv');
writetable(wide_table, wideCsv);

fprintf('\n=== SSA vs SSD (aggregated over seeds; SSA also over L) ===\n');
disp(wide_table);
plotFiles = plot_comparison_figures(long_table, wide_table, resultsDir);
fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n  %s\n', ...
    wideCsv, longCsv, plotFiles{1}, plotFiles{2}, plotFiles{3});

%% --- Local helpers ---
function agg = aggregate_by_signal_snr(trial_table)
    keys = unique(trial_table(:, {'signal_type', 'snr_db'}), 'rows');
    rows = repmat(init_agg_row(), height(keys), 1);
    for i = 1:height(keys)
        mask = trial_table.signal_type == keys.signal_type(i) & ...
            trial_table.snr_db == keys.snr_db(i);
        sub = trial_table(mask, :);
        rows(i) = pack_agg_row(keys.signal_type(i), keys.snr_db(i), sub);
    end
    agg = struct2table(rows);
end

function r = init_agg_row()
    r = struct( ...
        'signal_type', "", ...
        'snr_db', nan, ...
        'num_trials', nan, ...
        'failure_rate', nan, ...
        'mean_nmse', nan, ...
        'mean_abs_rho', nan, ...
        'mean_peak_error', nan, ...
        'mean_pairs', nan);
end

function r = pack_agg_row(signal_type, snr_db, sub)
    r = init_agg_row();
    r.signal_type = signal_type;
    r.snr_db = snr_db;
    r.num_trials = height(sub);
    r.failure_rate = mean(~sub.success);
    r.mean_nmse = mean(sub.mean_nmse, 'omitnan');
    r.mean_abs_rho = mean(sub.mean_abs_rho, 'omitnan');
    r.mean_peak_error = mean(sub.mean_peak_error, 'omitnan');
    r.mean_pairs = mean(sub.num_pairs, 'omitnan');
end

function wide = build_wide_table(ssa_agg, ssd_agg)
    keys = unique([ssa_agg(:, {'signal_type', 'snr_db'}); ssd_agg(:, {'signal_type', 'snr_db'})], 'rows');
    wide = table();
    wide.signal_type = keys.signal_type;
    wide.snr_db = keys.snr_db;
    n = height(keys);
    nan_col = nan(n, 1);

    wide.ssa_num_trials = nan_col;
    wide.ssa_failure_rate = nan_col;
    wide.ssa_mean_nmse = nan_col;
    wide.ssa_mean_abs_rho = nan_col;
    wide.ssa_mean_peak_error = nan_col;
    wide.ssa_mean_pairs = nan_col;

    wide.ssd_num_trials = nan_col;
    wide.ssd_failure_rate = nan_col;
    wide.ssd_mean_nmse = nan_col;
    wide.ssd_mean_abs_rho = nan_col;
    wide.ssd_mean_peak_error = nan_col;
    wide.ssd_mean_pairs = nan_col;

    wide.nmse_delta_ssd_minus_ssa = nan_col;

    for i = 1:n
        st = keys.signal_type(i);
        sn = keys.snr_db(i);
        ia = find(ssa_agg.signal_type == st & ssa_agg.snr_db == sn, 1);
        ib = find(ssd_agg.signal_type == st & ssd_agg.snr_db == sn, 1);
        if ~isempty(ia)
            wide.ssa_num_trials(i) = ssa_agg.num_trials(ia);
            wide.ssa_failure_rate(i) = ssa_agg.failure_rate(ia);
            wide.ssa_mean_nmse(i) = ssa_agg.mean_nmse(ia);
            wide.ssa_mean_abs_rho(i) = ssa_agg.mean_abs_rho(ia);
            wide.ssa_mean_peak_error(i) = ssa_agg.mean_peak_error(ia);
            wide.ssa_mean_pairs(i) = ssa_agg.mean_pairs(ia);
        end
        if ~isempty(ib)
            wide.ssd_num_trials(i) = ssd_agg.num_trials(ib);
            wide.ssd_failure_rate(i) = ssd_agg.failure_rate(ib);
            wide.ssd_mean_nmse(i) = ssd_agg.mean_nmse(ib);
            wide.ssd_mean_abs_rho(i) = ssd_agg.mean_abs_rho(ib);
            wide.ssd_mean_peak_error(i) = ssd_agg.mean_peak_error(ib);
            wide.ssd_mean_pairs(i) = ssd_agg.mean_pairs(ib);
        end
        if ~isempty(ia) && ~isempty(ib)
            wide.nmse_delta_ssd_minus_ssa(i) = ssd_agg.mean_nmse(ib) - ssa_agg.mean_nmse(ia);
        end
    end
end

function plotFiles = plot_comparison_figures(longTable, wideTable, resultsDir)
    signalTypes = ["close", "multiscale"];
    methods = ["SSA", "SSD"];
    methodLabels = ["SSA", "SSD (oscillatory-only)"];
    colors = [0.1216, 0.4667, 0.7059; 0.8392, 0.1529, 0.1569];

    nmseFig = figure('Name', 'SSA vs oscillatory-only SSD mean NMSE vs SNR', 'NumberTitle', 'off');
    tiledlayout(1, numel(signalTypes), 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signalTypes)
        ax = nexttile;
        hold(ax, 'on');
        subset = longTable(strcmpi(string(longTable.signal_type), signalTypes(si)), :);
        snrVals = sorted_snr(unique(subset.snr_db));
        x = 1:numel(snrVals);
        for mi = 1:numel(methods)
            methodMask = strcmpi(string(subset.method), methods(mi));
            methodSubset = subset(methodMask, :);
            y = nan(size(snrVals));
            for k = 1:numel(snrVals)
                rowMask = snr_equal(methodSubset.snr_db, snrVals(k));
                if any(rowMask)
                    y(k) = methodSubset.mean_nmse(find(rowMask, 1));
                end
            end
            plot(ax, x, y, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
                'Color', colors(mi, :), 'DisplayName', methodLabels(mi));
        end
        set(ax, 'YScale', 'log');
        xticks(ax, x);
        xticklabels(ax, snr_labels(snrVals));
        xlabel(ax, 'SNR (dB)');
        ylabel(ax, 'Mean NMSE');
        title(ax, char(signalTypes(si)));
        grid(ax, 'on');
        legend(ax, 'Location', 'best');
    end
    sgtitle('SSA vs oscillatory-only SSD: mean NMSE by signal type');
    nmsePng = fullfile(resultsDir, 'ssa_ssd_mean_nmse_vs_snr.png');
    exportgraphics(nmseFig, nmsePng, 'Resolution', 300);

    rhoFig = figure('Name', 'SSA vs oscillatory-only SSD mean abs rho vs SNR', 'NumberTitle', 'off');
    tiledlayout(1, numel(signalTypes), 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signalTypes)
        ax = nexttile;
        hold(ax, 'on');
        subset = longTable(strcmpi(string(longTable.signal_type), signalTypes(si)), :);
        snrVals = sorted_snr(unique(subset.snr_db));
        x = 1:numel(snrVals);
        for mi = 1:numel(methods)
            methodMask = strcmpi(string(subset.method), methods(mi));
            methodSubset = subset(methodMask, :);
            y = nan(size(snrVals));
            for k = 1:numel(snrVals)
                rowMask = snr_equal(methodSubset.snr_db, snrVals(k));
                if any(rowMask)
                    y(k) = methodSubset.mean_abs_rho(find(rowMask, 1));
                end
            end
            plot(ax, x, y, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6, ...
                'Color', colors(mi, :), 'DisplayName', methodLabels(mi));
        end
        xticks(ax, x);
        xticklabels(ax, snr_labels(snrVals));
        ylim(ax, [0, 1.05]);
        xlabel(ax, 'SNR (dB)');
        ylabel(ax, 'Mean |Pearson \rho|');
        title(ax, char(signalTypes(si)));
        grid(ax, 'on');
        legend(ax, 'Location', 'best');
    end
    sgtitle('SSA vs oscillatory-only SSD: mean absolute correlation by signal type');
    rhoPng = fullfile(resultsDir, 'ssa_ssd_mean_abs_rho_vs_snr.png');
    exportgraphics(rhoFig, rhoPng, 'Resolution', 300);

    deltaFig = figure('Name', 'Oscillatory-only SSD minus SSA NMSE delta vs SNR', 'NumberTitle', 'off');
    tiledlayout(1, numel(signalTypes), 'Padding', 'compact', 'TileSpacing', 'compact');
    for si = 1:numel(signalTypes)
        ax = nexttile;
        subset = wideTable(strcmpi(string(wideTable.signal_type), signalTypes(si)), :);
        snrVals = sorted_snr(subset.snr_db);
        y = nan(size(snrVals));
        for k = 1:numel(snrVals)
            rowMask = snr_equal(subset.snr_db, snrVals(k));
            if any(rowMask)
                y(k) = subset.nmse_delta_ssd_minus_ssa(find(rowMask, 1));
            end
        end
        x = 1:numel(snrVals);
        bar(ax, x, y, 0.65, 'FaceColor', [0.4, 0.4, 0.4]);
        yline(ax, 0, 'k-', 'LineWidth', 1);
        xticks(ax, x);
        xticklabels(ax, snr_labels(snrVals));
        xlabel(ax, 'SNR (dB)');
        ylabel(ax, '\Delta NMSE (oscillatory-only SSD - SSA)');
        title(ax, char(signalTypes(si)));
        grid(ax, 'on');
    end
    sgtitle('SSA vs oscillatory-only SSD: NMSE difference by signal type');
    deltaPng = fullfile(resultsDir, 'ssa_ssd_nmse_delta_vs_snr.png');
    exportgraphics(deltaFig, deltaPng, 'Resolution', 300);

    plotFiles = {nmsePng, rhoPng, deltaPng};
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
        if isinf(v(i))
            labels{i} = 'Inf';
        else
            labels{i} = num2str(v(i));
        end
    end
end

function mask = snr_equal(v, target)
    if isinf(target)
        mask = isinf(v);
    else
        mask = v == target;
    end
end
