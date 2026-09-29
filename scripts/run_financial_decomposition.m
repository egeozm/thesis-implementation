%RUN_FINANCIAL_DECOMPOSITION Apply SSA, SSD, and TensorCPD to RQ3 data.
%
%   Loads the EUR/USD and S&P 500 level series used for the thesis financial
%   experiments, standardizes each asset separately, runs the three
%   decomposition families, and writes component tables plus diagnostics.

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));
addpath(fullfile(projRoot, 'ssd'));
addpath(fullfile(projRoot, 'tensor'));
if exist(fullfile(projRoot, 'tensorlab'), 'dir')
    addpath(genpath(fullfile(projRoot, 'tensorlab')));
end

resultsDir = fullfile(projRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

calibrationMat = fullfile(resultsDir, 'ssa_calibration_results.mat');
require_file(calibrationMat, 'scripts/run_ssa_calibration.m');
loaded = load(calibrationMat, 'best_params');
best_params = loaded.best_params;

config = build_config(projRoot);
diagnostic_rows = repmat(init_diagnostic_row(), 0, 1);
ssd_comparison_rows = repmat(init_ssd_comparison_row(), 0, 1);
saved_files = strings(0, 1);

fprintf('\n=== Financial decomposition (RQ3) ===\n');
fprintf('Date range : %s to %s\n', datestr(config.date_range(1)), datestr(config.date_range(2)));
fprintf('SSA L      : %s\n', char(numvec_to_str(config.ssa.L_list)));
fprintf('SSD L_list : %s\n', char(numvec_to_str(config.ssd.L_list)));
fprintf('Tensor grid: strides %s, L %s, R %s, grouping %s\n', ...
    char(numvec_to_str(config.tensor.strides)), ...
    char(numvec_to_str(config.tensor.base_window_lengths)), ...
    char(numvec_to_str(config.tensor.cpd_ranks)), ...
    char(strjoin(config.tensor.grouping_modes, ', ')));

for a = 1:numel(config.assets)
    asset = config.assets(a);
    fprintf('\n--- %s ---\n', char(asset.label));

    series = load_financial_series(asset.path, asset.date_col, asset.value_col, config.date_range);
    fprintf('Loaded %d observations from %s to %s.\n', numel(series.x_z), ...
        datestr(series.dates(1)), datestr(series.dates(end)));

    [ssa_table, ssa_diag] = run_ssa_financial(asset, series, best_params, config);
    ssa_file = fullfile(resultsDir, sprintf('financial_%s_ssa.csv', char(asset.id)));
    writetable(ssa_table, ssa_file);
    saved_files(end + 1) = string(ssa_file); %#ok<SAGROW>
    diagnostic_rows = append_rows(diagnostic_rows, ssa_diag);

    [ssd_table, ssd_diag, ssd_comparison] = run_ssd_financial(asset, series, best_params, config);
    ssd_file = fullfile(resultsDir, sprintf('financial_%s_ssd.csv', char(asset.id)));
    writetable(ssd_table, ssd_file);
    saved_files(end + 1) = string(ssd_file); %#ok<SAGROW>
    diagnostic_rows = append_rows(diagnostic_rows, ssd_diag);
    ssd_comparison_rows = append_rows(ssd_comparison_rows, ssd_comparison);

    [tensor_table, tensor_diag] = run_tensor_financial(asset, series, config);
    tensor_file = fullfile(resultsDir, sprintf('financial_%s_tensor.csv', char(asset.id)));
    writetable(tensor_table, tensor_file);
    saved_files(end + 1) = string(tensor_file); %#ok<SAGROW>
    diagnostic_rows = append_rows(diagnostic_rows, tensor_diag);
end

diagnostics_table = struct2table(diagnostic_rows);
diagnostics_file = fullfile(resultsDir, 'financial_diagnostics.csv');
writetable(diagnostics_table, diagnostics_file);
saved_files(end + 1) = string(diagnostics_file);

ssd_comparison_table = struct2table(ssd_comparison_rows);
ssd_comparison_file = fullfile(resultsDir, 'financial_ssd_variant_comparison.csv');
writetable(ssd_comparison_table, ssd_comparison_file);
saved_files(end + 1) = string(ssd_comparison_file);

fprintf('\nDiagnostics:\n');
disp(diagnostics_table);
fprintf('\nSaved:\n');
for i = 1:numel(saved_files)
    fprintf('  %s\n', char(saved_files(i)));
end

%% --- Local helpers ---
function config = build_config(projRoot)
    dataDir = fullfile(projRoot, 'financial time series datasets');

    config = struct();
    config.fs = 1;
    config.date_range = [datetime(1999, 1, 4), datetime(2026, 4, 2)];
    config.assets = [ ...
        asset_struct('eurusd', 'EUR/USD', ...
            fullfile(dataDir, 'ECB Data Portal_20260407115419.csv'), ...
            'DATE', 'US dollar/Euro ECB reference exchange rate (EXR.D.USD.EUR.SP00.A)'), ...
        asset_struct('spx', 'S&P 500', ...
            fullfile(dataDir, '_spx_d.csv'), ...
            'Date', 'Close') ...
        ];

    config.ssa = struct('L_list', [200, 400, 600, 1000]);
    config.ssd = struct( ...
        'L_list', [200, 400, 600, 1000], ...
        'extract_trend_first', true);
    config.tensor = struct();
    config.tensor.strides = [1, 2];
    config.tensor.base_window_lengths = 200;
    config.tensor.cpd_ranks = [1, 2, 3, 4, 6, 8];
    config.tensor.grouping_modes = ["auto", "none"];
    config.tensor.reconstruction_mode = 'cpd';
    config.tensor.num_restarts = 5;
    config.tensor.restart_seed_base = 7300;
    config.tensor.prefer_stable_candidate = true;
end

function asset = asset_struct(id, label, path, date_col, value_col)
    asset = struct( ...
        'id', string(id), ...
        'label', string(label), ...
        'path', string(path), ...
        'date_col', string(date_col), ...
        'value_col', string(value_col));
end

function series = load_financial_series(path, date_col, value_col, date_range)
    require_file(path, 'financial data download');

    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');
    T = readtable(path, opts);
    date_col = char(date_col);
    value_col = char(value_col);

    if ~ismember(date_col, T.Properties.VariableNames)
        error('run_financial_decomposition:MissingDateColumn', ...
            'Missing date column "%s" in %s.', date_col, path);
    end
    if ~ismember(value_col, T.Properties.VariableNames)
        error('run_financial_decomposition:MissingValueColumn', ...
            'Missing value column "%s" in %s.', value_col, path);
    end

    dates = parse_dates(T.(date_col));
    values = parse_numeric_values(T.(value_col));
    keep = dates >= date_range(1) & dates <= date_range(2) & ~isnan(values);

    dates = dates(keep);
    values = values(keep);
    [dates, order] = sort(dates);
    values = values(order);

    mu = mean(values, 'omitnan');
    sigma = std(values, 0, 'omitnan');
    if sigma <= eps
        error('run_financial_decomposition:ConstantSeries', ...
            'Cannot z-score %s because its standard deviation is zero.', path);
    end

    series = struct();
    series.dates = dates(:);
    series.x_raw = values(:);
    series.x_z = (values(:) - mu) ./ sigma;
    series.mean = mu;
    series.std = sigma;
end

function dates = parse_dates(raw_dates)
    if isdatetime(raw_dates)
        dates = raw_dates;
        return;
    end
    dates = datetime(string(raw_dates), 'InputFormat', 'yyyy-MM-dd');
end

function values = parse_numeric_values(raw_values)
    if isnumeric(raw_values)
        values = double(raw_values);
        return;
    end
    values = str2double(string(raw_values));
end

function [out_table, diag_rows] = run_ssa_financial(asset, series, best_params, config)
    x = series.x_z(:);
    out_table = base_series_table(series);
    diag_rows = repmat(init_diagnostic_row(), numel(config.ssa.L_list), 1);

    for i = 1:numel(config.ssa.L_list)
        L = config.ssa.L_list(i);
        row = init_diagnostic_row();
        row.asset = asset.id;
        row.asset_label = asset.label;
        row.method = "SSA";
        row.L = L;

        try
            features = ssa_precompute(x, L, best_params.r_max, config.fs);
            grouping = auto_group_ssa(features, best_params);
            groups = ssa_financial_groups(grouping);
            reconstructed = ssa_reconstruct_groups(features, groups);

            out_table.(sprintf('trend_L%d', L)) = reconstructed(:, 1);
            out_table.(sprintf('pair_lo_L%d', L)) = reconstructed(:, 2);
            out_table.(sprintf('pair_hi_L%d', L)) = reconstructed(:, 3);
            out_table.(sprintf('residual_L%d', L)) = reconstructed(:, 4);

            row.success = true;
            row.num_components = size(reconstructed, 2);
            row.num_oscillatory_pairs = numel(grouping.oscillatory_pairs);
            row.pair_frequencies = join_numbers(grouping.pair_frequencies);
            row.pair_scores = join_numbers(grouping.pair_scores);
            row = add_reconstruction_diagnostics(row, x, sum(reconstructed(:, 1:3), 2), reconstructed(:, 4));
        catch ME
            out_table.(sprintf('trend_L%d', L)) = nan(size(x));
            out_table.(sprintf('pair_lo_L%d', L)) = nan(size(x));
            out_table.(sprintf('pair_hi_L%d', L)) = nan(size(x));
            out_table.(sprintf('residual_L%d', L)) = nan(size(x));
            row.success = false;
            row.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
        end

        diag_rows(i) = row;
    end
end

function groups = ssa_financial_groups(grouping)
    groups = cell(1, 4);
    groups{1} = grouping.trend_idx;
    for k = 1:2
        if numel(grouping.oscillatory_pairs) >= k
            groups{k + 1} = grouping.oscillatory_pairs{k};
        else
            groups{k + 1} = [];
        end
    end
    groups{4} = grouping.residual_idx;
end

function [out_table, row, comparison_rows] = run_ssd_financial(asset, series, best_params, config)
    x = series.x_z(:);
    out_table = base_series_table(series);
    row = init_diagnostic_row();
    row.asset = asset.id;
    row.asset_label = asset.label;
    row.method = "trend-enabled SSD";
    row.L_list = join_numbers(config.ssd.L_list);
    comparison_rows = run_ssd_variant_comparison(asset, x, best_params, config);

    try
        params = ssd_merge_params(best_params, struct( ...
            'fs', config.fs, ...
            'L_list', config.ssd.L_list, ...
            'extract_trend_first', config.ssd.extract_trend_first));
        out = ssd_decompose(x, params);
        out_table.trend = out.trend(:);
        for k = 1:size(out.modes, 2)
            out_table.(sprintf('mode_%d', k)) = out.modes(:, k);
        end
        out_table.residual = out.residual(:);

        row.success = true;
        row.num_components = double(out.trend_extracted) + size(out.modes, 2) + 1;
        row.trend_extracted = out.trend_extracted;
        row.trend_L = out.trend_L;
        row.trend_idx = join_numbers(out.trend_idx);
        row.num_iterations = out.num_iterations;
        row.L_history = join_numbers(out.L_history);
        row.pair_history = join_pairs(out.pair_history);
        row.pair_scores = join_numbers(out.scores);
        row = add_reconstruction_diagnostics(row, x, ssd_total_reconstruction(out), out.residual);
    catch ME
        out_table.trend = nan(size(x));
        out_table.residual = nan(size(x));
        row.success = false;
        row.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
    end
end

function comparison_rows = run_ssd_variant_comparison(asset, x, best_params, config)
    variants = [ ...
        ssd_variant_struct("initial SSD (oscillatory-only)", false), ...
        ssd_variant_struct("trend-enabled SSD", true) ...
        ];
    comparison_rows = repmat(init_ssd_comparison_row(), numel(variants), 1);

    for i = 1:numel(variants)
        v = variants(i);
        row = init_ssd_comparison_row();
        row.asset = asset.id;
        row.asset_label = asset.label;
        row.method_variant = v.label;
        row.extract_trend_first = v.extract_trend_first;
        row.L_list = join_numbers(config.ssd.L_list);

        try
            params = ssd_merge_params(best_params, struct( ...
                'fs', config.fs, ...
                'L_list', config.ssd.L_list, ...
                'extract_trend_first', v.extract_trend_first));
            out = ssd_decompose(x, params);
            reconstruction = ssd_total_reconstruction(out);

            row.success = true;
            row.num_components = double(out.trend_extracted) + size(out.modes, 2) + 1;
            row.trend_extracted = out.trend_extracted;
            row.trend_L = out.trend_L;
            row.trend_idx = join_numbers(out.trend_idx);
            row.num_iterations = out.num_iterations;
            row.L_history = join_numbers(out.L_history);
            row.pair_history = join_pairs(out.pair_history);
            row.pair_scores = join_numbers(out.scores);
            row = add_ssd_comparison_diagnostics(row, x, reconstruction, out.residual);
        catch ME
            row.success = false;
            row.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
        end

        comparison_rows(i) = row;
    end
end

function v = ssd_variant_struct(label, extract_trend_first)
    v = struct( ...
        'label', string(label), ...
        'extract_trend_first', extract_trend_first);
end

function reconstruction = ssd_total_reconstruction(out)
    reconstruction = out.trend(:);
    if ~isempty(out.modes)
        reconstruction = reconstruction + sum(out.modes, 2);
    end
end

function [out_table, diag_rows] = run_tensor_financial(asset, series, config)
    x = series.x_z(:);
    out_table = base_series_table(series);
    diag_rows = repmat(init_diagnostic_row(), estimate_tensor_trials(config), 1);
    row_idx = 0;

    for L_base = config.tensor.base_window_lengths
        for R = config.tensor.cpd_ranks
            for gm = config.tensor.grouping_modes
                row_idx = row_idx + 1;
                prefix = sprintf('tensor_%s_L%d_R%d_%s', ...
                    stride_label_from_values(config.tensor.strides), L_base, R, char(gm));
                [setting_table, row] = run_tensor_setting(asset, x, config, L_base, R, gm, prefix);
                out_table = append_table_columns(out_table, setting_table);
                diag_rows(row_idx) = row;
            end
        end
    end

    diag_rows = diag_rows(1:row_idx);
    diag_rows = annotate_tensor_rank_recommendations(diag_rows);
end

function n = estimate_tensor_trials(config)
    n = numel(config.tensor.base_window_lengths) * ...
        numel(config.tensor.cpd_ranks) * ...
        numel(config.tensor.grouping_modes);
end

function rows = annotate_tensor_rank_recommendations(rows)
    if isempty(rows)
        return;
    end

    for i = 1:numel(rows)
        rows(i).rank_search_mode = "low_rank_grid";
    end

    grouping_modes = unique([rows.grouping_mode]);
    for gm = grouping_modes
        group_idx = find([rows.grouping_mode] == gm);
        if isempty(group_idx)
            continue;
        end

        group_rows = rows(group_idx);
        success = [group_rows.success];
        stable = success & ~[group_rows.unstable_terms];
        interpretable_stable = stable & ~[group_rows.low_interpretability];
        has_stable = any(stable);

        for k = group_idx
            rows(k).stable_rank_available = has_stable;
        end

        if any(interpretable_stable)
            eligible = group_idx(interpretable_stable);
            selected_idx = select_by_rel_error(rows, eligible);
        elseif has_stable
            eligible = group_idx(stable);
            selected_idx = select_by_rel_error(rows, eligible);
        elseif any(success)
            eligible = group_idx(success);
            selected_idx = select_by_cancellation_then_error(rows, eligible);
        else
            selected_idx = group_idx(1);
        end

        rows(selected_idx).recommended_tensor_setting = true;
    end
end

function selected_idx = select_by_rel_error(rows, eligible)
    rel_errors = safe_numeric([rows(eligible).rel_error], inf);
    [~, pos] = min(rel_errors);
    selected_idx = eligible(pos);
end

function selected_idx = select_by_cancellation_then_error(rows, eligible)
    cancellation = safe_numeric([rows(eligible).cancellation_ratio], inf);
    rel_errors = safe_numeric([rows(eligible).rel_error], inf);
    scores = [cancellation(:), rel_errors(:)];
    [~, order] = sortrows(scores, [1 2]);
    selected_idx = eligible(order(1));
end

function values = safe_numeric(values, replacement)
    values(~isfinite(values)) = replacement;
end

function [component_table, row] = run_tensor_setting(asset, x, config, L_base, R, grouping_mode, prefix)
    component_table = table();
    row = init_diagnostic_row();
    row.asset = asset.id;
    row.asset_label = asset.label;
    row.method = "TensorCPD";
    row.stride_label = stride_label_from_values(config.tensor.strides);
    row.base_window_length = L_base;
    row.cpd_rank = R;
    row.grouping_mode = string(grouping_mode);

    params = tensor_merge_params(struct( ...
        'strides', config.tensor.strides, ...
        'base_window_length', L_base, ...
        'cpd_rank', R, ...
        'grouping_mode', char(grouping_mode), ...
        'reconstruction_mode', config.tensor.reconstruction_mode, ...
        'cpd_method', 'cpd', ...
        'cpd_options', struct(), ...
        'fs', config.fs));

    try
        tensor_out = tensor_build_decimated_hankel(x, params);
        selected = select_financial_cpd_candidate(tensor_out, x, params, config.tensor);
        row.cpd_restart_count = selected.restart_count;
        row.cpd_selected_restart = selected.restart_idx;
        row.cpd_stable_restart_count = selected.stable_count;
        row.cpd_candidate_summary = selected.summary;
        if ~selected.success
            row.success = false;
            row.failure_reason = selected.failure_reason;
            return;
        end

        cpd_out = selected.cpd_out;
        recon = selected.recon;
        row.cpd_tensor_fit = cpd_out.fit;
        row.cpd_tensor_rel_error = cpd_out.rel_error;
        row.support_length = recon.support_length;

        raw_total_recon = sum(recon.est_components_full, 2);
        raw_residual = recon.residual_full(:);
        row = add_cpd_stability_diagnostics(row, x, recon.est_components_full, raw_total_recon);

        if grouping_mode == "auto"
            [components, grouping] = group_financial_auto_terms( ...
                recon.est_components_full, cpd_out.lambda, params);
            row.cluster_count = size(components, 2);
            row.auto_group_collapsed = grouping.auto_group_collapsed;
        else
            components = recon.est_components_full;
            grouping = struct( ...
                'mode', "none", ...
                'num_grouped_components', size(components, 2), ...
                'auto_group_collapsed', false);
            row.cluster_count = size(components, 2);
        end

        for k = 1:size(components, 2)
            component_table.(sprintf('%s_c%d', prefix, k)) = components(:, k);
        end
        component_table.(sprintf('%s_total_recon', prefix)) = raw_total_recon;
        component_table.(sprintf('%s_residual', prefix)) = raw_residual;

        row.success = true;
        row.num_components = size(components, 2);
        row = add_reconstruction_diagnostics(row, x, raw_total_recon, raw_residual);
        if isfield(grouping, 'num_grouped_components')
            row.cluster_count = grouping.num_grouped_components;
        end
        row = add_interpretability_diagnostics(row);
    catch ME
        row.success = false;
        row.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
    end
end

function selected = select_financial_cpd_candidate(tensor_out, x, params, tensor_config)
    restart_count = max(1, double(tensor_config.num_restarts));
    candidates = repmat(init_cpd_candidate(), restart_count, 1);

    for restart_idx = 1:restart_count
        seed = tensor_config.restart_seed_base + ...
            100000 * params.base_window_length + 1000 * params.cpd_rank + restart_idx;
        params_i = params;
        params_i.cpd_initialization = seeded_cpd_initialization( ...
            tensor_out.tensor_size, params.cpd_rank, seed);
        candidates(restart_idx) = run_financial_cpd_candidate( ...
            tensor_out, x, params_i, restart_idx, seed);
    end

    selected = choose_financial_cpd_candidate(candidates, tensor_config.prefer_stable_candidate);
    selected.restart_count = restart_count;
    selected.stable_count = sum([candidates.success] & ~[candidates.unstable_terms]);
    selected.summary = summarize_cpd_candidates(candidates);
end

function candidate = run_financial_cpd_candidate(tensor_out, x, params, restart_idx, seed)
    candidate = init_cpd_candidate();
    candidate.restart_idx = restart_idx;
    candidate.seed = seed;

    cpd_out = tensor_run_cpd(tensor_out, params);
    candidate.cpd_out = cpd_out;
    candidate.cpd_tensor_fit = cpd_out.fit;
    candidate.cpd_tensor_rel_error = cpd_out.rel_error;
    if ~cpd_out.success
        candidate.failure_reason = string(cpd_out.failure_reason);
        return;
    end

    recon = tensor_reconstruct_via_cpd(cpd_out, x, params);
    candidate.recon = recon;
    if ~recon.success
        candidate.failure_reason = string(recon.failure_reason);
        return;
    end

    raw_total_recon = sum(recon.est_components_full, 2);
    raw_residual = recon.residual_full(:);
    diag_row = init_diagnostic_row();
    diag_row.cpd_tensor_fit = cpd_out.fit;
    diag_row.cpd_tensor_rel_error = cpd_out.rel_error;
    diag_row = add_cpd_stability_diagnostics(diag_row, x, recon.est_components_full, raw_total_recon);
    diag_row = add_reconstruction_diagnostics(diag_row, x, raw_total_recon, raw_residual);

    candidate.success = true;
    candidate.failure_reason = "";
    candidate.rel_error = diag_row.rel_error;
    candidate.fit = diag_row.fit;
    candidate.unstable_terms = diag_row.unstable_terms;
    candidate.cancellation_ratio = diag_row.cancellation_ratio;
    candidate.max_term_to_signal_norm = diag_row.max_term_to_signal_norm;
    candidate.max_abs_term = diag_row.max_abs_term;
end

function selected = choose_financial_cpd_candidate(candidates, prefer_stable_candidate)
    success_mask = [candidates.success];
    if ~any(success_mask)
        selected = candidates(1);
        selected.failure_reason = summarize_candidate_failures(candidates);
        selected.restart_idx = nan;
        return;
    end

    eligible = find(success_mask);
    stable_mask = success_mask & ~[candidates.unstable_terms];
    if prefer_stable_candidate && any(stable_mask)
        eligible = find(stable_mask);
    end

    rel_errors = [candidates(eligible).rel_error];
    [~, best_pos] = min(rel_errors);
    selected = candidates(eligible(best_pos));
end

function U0 = seeded_cpd_initialization(tensor_size, rank_k, seed)
    old_rng = rng;
    cleanup = onCleanup(@() rng(old_rng)); %#ok<NASGU>
    rng(seed, 'twister');
    if exist('cpd_rnd', 'file') == 2 || exist('cpd_rnd', 'builtin') == 5
        U0 = cpd_rnd(tensor_size, rank_k);
    else
        U0 = cell(numel(tensor_size), 1);
        for m = 1:numel(tensor_size)
            U0{m} = randn(tensor_size(m), rank_k);
        end
    end
end

function candidate = init_cpd_candidate()
    candidate = struct( ...
        'success', false, ...
        'failure_reason', "", ...
        'restart_idx', nan, ...
        'seed', nan, ...
        'cpd_out', struct(), ...
        'recon', struct(), ...
        'fit', nan, ...
        'rel_error', nan, ...
        'cpd_tensor_fit', nan, ...
        'cpd_tensor_rel_error', nan, ...
        'unstable_terms', false, ...
        'cancellation_ratio', nan, ...
        'max_term_to_signal_norm', nan, ...
        'max_abs_term', nan, ...
        'restart_count', nan, ...
        'stable_count', nan, ...
        'summary', "");
end

function s = summarize_cpd_candidates(candidates)
    parts = strings(1, numel(candidates));
    for i = 1:numel(candidates)
        c = candidates(i);
        if c.success
            status = "stable";
            if c.unstable_terms
                status = "unstable";
            end
            parts(i) = sprintf('#%d:%s,rel=%.4g,cancel=%.4g', ...
                c.restart_idx, status, c.rel_error, c.cancellation_ratio);
        else
            parts(i) = sprintf('#%d:failed,%s', c.restart_idx, char(c.failure_reason));
        end
    end
    s = strjoin(parts, ' | ');
end

function s = summarize_candidate_failures(candidates)
    reasons = strings(1, numel(candidates));
    for i = 1:numel(candidates)
        reasons(i) = sprintf('#%d: %s', ...
            candidates(i).restart_idx, char(candidates(i).failure_reason));
    end
    s = "All CPD restarts failed: " + strjoin(reasons, ' | ');
end

function [components, grouping] = group_financial_auto_terms(rank1_components, lambda, params)
    [~, grouping] = group_cp_terms_auto(rank1_components, lambda, params);
    num_groups = grouping.num_grouped_components;
    components = zeros(size(rank1_components, 1), num_groups);

    for k = 1:num_groups
        members = grouping.clusters(k).members;
        components(:, k) = sum(rank1_components(:, members), 2);
    end

    grouping.financial_sum_mode = "preserve_original_cpd_signs";
    grouping.auto_group_collapsed = num_groups <= 1 && size(rank1_components, 2) > 1;
end

function T = base_series_table(series)
    T = table(series.dates(:), series.x_raw(:), series.x_z(:), ...
        'VariableNames', {'date', 'x_raw', 'x_z'});
end

function out = append_table_columns(out, extra)
    for i = 1:numel(extra.Properties.VariableNames)
        name = extra.Properties.VariableNames{i};
        out.(name) = extra.(name);
    end
end

function rows = append_rows(rows, new_rows)
    if isempty(rows)
        rows = new_rows(:);
    else
        rows = [rows(:); new_rows(:)]; %#ok<AGROW>
    end
end

function row = add_reconstruction_diagnostics(row, x, reconstruction, residual)
    x = x(:);
    reconstruction = reconstruction(:);
    residual = residual(:);
    x_norm = norm(x);
    row.x_z_norm = x_norm;
    row.reconstruction_norm = norm(reconstruction);
    row.residual_norm = norm(residual);
    if x_norm <= eps
        row.rel_error = nan;
        row.relative_error = nan;
        row.fit = nan;
    else
        row.rel_error = row.residual_norm / x_norm;
        row.relative_error = row.rel_error;
        row.fit = 1 - row.rel_error;
    end
    row.max_abs_reconstruction = max(abs(reconstruction), [], 'omitnan');
end

function row = add_cpd_stability_diagnostics(row, x, rank1_components, total_recon)
    x = x(:);
    rank1_components = rank1_components(:, :);
    total_recon = total_recon(:);

    if isempty(rank1_components)
        return;
    end

    term_norms = vecnorm(rank1_components, 2, 1);
    total_norm = norm(total_recon);
    x_norm = norm(x);
    max_abs_x = max(abs(x), [], 'omitnan');

    row.term_norm_sum = sum(term_norms, 'omitnan');
    row.max_term_norm = max(term_norms, [], 'omitnan');
    row.max_abs_term = max(abs(rank1_components(:)), [], 'omitnan');
    row.cancellation_ratio = row.term_norm_sum / max(total_norm, eps);
    row.max_term_to_signal_norm = row.max_term_norm / max(x_norm, eps);
    row.unstable_terms = row.cancellation_ratio > 10 || ...
        row.max_term_to_signal_norm > 5 || ...
        row.max_abs_term > 10 * max(max_abs_x, eps);
end

function row = add_interpretability_diagnostics(row)
    reasons = strings(0, 1);
    if row.unstable_terms
        reasons(end + 1) = "unstable cancellation-heavy CP terms"; %#ok<AGROW>
    end
    if row.auto_group_collapsed
        reasons(end + 1) = "auto grouping collapsed CP terms into one component"; %#ok<AGROW>
    end

    row.low_interpretability = ~isempty(reasons);
    row.interpretability_reason = strjoin(reasons, "; ");
end

function row = add_ssd_comparison_diagnostics(row, x, reconstruction, residual)
    x = x(:);
    reconstruction = reconstruction(:);
    residual = residual(:);
    x_norm = norm(x);
    row.x_z_norm = x_norm;
    row.reconstruction_norm = norm(reconstruction);
    row.residual_norm = norm(residual);
    if x_norm <= eps
        row.relative_error = nan;
        row.fit = nan;
    else
        row.relative_error = row.residual_norm / x_norm;
        row.fit = 1 - row.relative_error;
    end
    row.max_abs_reconstruction = max(abs(reconstruction), [], 'omitnan');
end

function row = init_diagnostic_row()
    row = struct( ...
        'asset', "", ...
        'asset_label', "", ...
        'method', "", ...
        'L', nan, ...
        'L_list', "", ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'grouping_mode', "", ...
        'num_components', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'trend_extracted', false, ...
        'trend_L', nan, ...
        'trend_idx', "", ...
        'num_oscillatory_pairs', nan, ...
        'pair_frequencies', "", ...
        'pair_scores', "", ...
        'num_iterations', nan, ...
        'L_history', "", ...
        'pair_history', "", ...
        'fit', nan, ...
        'rel_error', nan, ...
        'relative_error', nan, ...
        'cpd_tensor_fit', nan, ...
        'cpd_tensor_rel_error', nan, ...
        'cpd_restart_count', nan, ...
        'cpd_selected_restart', nan, ...
        'cpd_stable_restart_count', nan, ...
        'cpd_candidate_summary', "", ...
        'rank_search_mode', "", ...
        'stable_rank_available', false, ...
        'recommended_tensor_setting', false, ...
        'x_z_norm', nan, ...
        'reconstruction_norm', nan, ...
        'max_abs_reconstruction', nan, ...
        'support_length', nan, ...
        'cluster_count', nan, ...
        'term_norm_sum', nan, ...
        'max_term_norm', nan, ...
        'max_abs_term', nan, ...
        'cancellation_ratio', nan, ...
        'max_term_to_signal_norm', nan, ...
        'unstable_terms', false, ...
        'auto_group_collapsed', false, ...
        'low_interpretability', false, ...
        'interpretability_reason', "", ...
        'residual_norm', nan);
end

function row = init_ssd_comparison_row()
    row = struct( ...
        'asset', "", ...
        'asset_label', "", ...
        'method_variant', "", ...
        'extract_trend_first', false, ...
        'L_list', "", ...
        'success', false, ...
        'failure_reason', "", ...
        'num_components', nan, ...
        'trend_extracted', false, ...
        'trend_L', nan, ...
        'trend_idx', "", ...
        'num_iterations', nan, ...
        'L_history', "", ...
        'pair_history', "", ...
        'pair_scores', "", ...
        'x_z_norm', nan, ...
        'reconstruction_norm', nan, ...
        'residual_norm', nan, ...
        'relative_error', nan, ...
        'fit', nan, ...
        'max_abs_reconstruction', nan);
end

function require_file(path_str, producer)
    if exist(path_str, 'file') ~= 2
        error('run_financial_decomposition:MissingFile', ...
            'Missing %s. Run or provide %s first.', char(path_str), char(producer));
    end
end

function label = stride_label_from_values(v)
    label = char("S" + join(string(v), "_S"));
end

function s = join_numbers(v)
    if isempty(v)
        s = "";
        return;
    end
    s = strjoin(compose('%.12g', v(:).'), ';');
end

function s = join_pairs(pairs)
    if isempty(pairs)
        s = "";
        return;
    end
    parts = strings(1, numel(pairs));
    for i = 1:numel(pairs)
        parts(i) = "[" + join(string(pairs{i}), " ") + "]";
    end
    s = strjoin(parts, ';');
end

function s = numvec_to_str(v)
    s = strjoin(string(v), ', ');
end
