%RUN_FINANCIAL_VISUALIZATIONS Render qualitative RQ3 financial figures.
%
%   Consumes the CSV outputs from scripts/run_financial_decomposition.m and
%   exports per-asset overview plots for SSA, SSD, TensorCPD, plus a compact
%   cross-method comparison panel.

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
resultsDir = fullfile(projRoot, 'results');

diagnosticsFile = fullfile(resultsDir, 'financial_diagnostics.csv');
require_file(diagnosticsFile, 'scripts/run_financial_decomposition.m');
diagnostics = readtable(diagnosticsFile, 'VariableNamingRule', 'preserve');
diagnostics = normalize_diagnostics(diagnostics);

assets = [ ...
    asset_struct('eurusd', 'EUR/USD'), ...
    asset_struct('spx', 'S&P 500') ...
    ];

saved_files = strings(0, 1);

fprintf('\n=== Financial visualizations (RQ3) ===\n');
for a = 1:numel(assets)
    asset = assets(a);
    fprintf('\n--- %s ---\n', char(asset.label));

    ssa_table = read_required_table(resultsDir, asset.id, 'ssa');
    ssd_table = read_required_table(resultsDir, asset.id, 'ssd');
    tensor_table = read_required_table(resultsDir, asset.id, 'tensor');
    asset_diag = diagnostics(diagnostics.asset == string(asset.id), :);

    fig = make_ssa_overview(asset, ssa_table, asset_diag);
    out = fullfile(resultsDir, sprintf('financial_%s_ssa_overview.png', char(asset.id)));
    exportgraphics(fig, out, 'Resolution', 300);
    close(fig);
    saved_files(end + 1) = string(out); %#ok<SAGROW>

    fig = make_ssd_overview(asset, ssd_table, asset_diag);
    out = fullfile(resultsDir, sprintf('financial_%s_ssd_overview.png', char(asset.id)));
    exportgraphics(fig, out, 'Resolution', 300);
    close(fig);
    saved_files(end + 1) = string(out); %#ok<SAGROW>

    reps = choose_tensor_representatives(asset_diag);
    fig = make_tensor_overview(asset, tensor_table, reps, asset_diag);
    out = fullfile(resultsDir, sprintf('financial_%s_tensor_overview.png', char(asset.id)));
    exportgraphics(fig, out, 'Resolution', 300);
    close(fig);
    saved_files(end + 1) = string(out); %#ok<SAGROW>

    fig = make_method_compare(asset, ssa_table, ssd_table, tensor_table, asset_diag, reps);
    out = fullfile(resultsDir, sprintf('financial_%s_methods_compare.png', char(asset.id)));
    exportgraphics(fig, out, 'Resolution', 300);
    close(fig);
    saved_files(end + 1) = string(out); %#ok<SAGROW>

    fig = make_tensor_rank_contrast(asset, tensor_table, asset_diag);
    out = fullfile(resultsDir, sprintf('financial_%s_tensor_best_fit_vs_stable_rank.png', char(asset.id)));
    exportgraphics(fig, out, 'Resolution', 300);
    close(fig);
    saved_files(end + 1) = string(out); %#ok<SAGROW>
end

fprintf('\nSaved:\n');
for i = 1:numel(saved_files)
    fprintf('  %s\n', char(saved_files(i)));
end

%% --- Local helpers ---
function asset = asset_struct(id, label)
    asset = struct('id', string(id), 'label', string(label));
end

function T = normalize_diagnostics(T)
    string_vars = {'asset', 'asset_label', 'method', 'L_list', 'stride_label', ...
        'grouping_mode', 'failure_reason', 'pair_frequencies', 'pair_scores', ...
        'L_history', 'pair_history', 'interpretability_reason', ...
        'cpd_candidate_summary', 'rank_search_mode'};
    for i = 1:numel(string_vars)
        name = string_vars{i};
        if ismember(name, T.Properties.VariableNames)
            T.(name) = string(T.(name));
        end
    end
    logical_vars = {'success', 'unstable_terms', 'auto_group_collapsed', ...
        'low_interpretability', 'stable_rank_available', 'recommended_tensor_setting'};
    for i = 1:numel(logical_vars)
        name = logical_vars{i};
        if ismember(name, T.Properties.VariableNames) && ~islogical(T.(name))
            if isnumeric(T.(name))
                T.(name) = T.(name) ~= 0;
            else
                T.(name) = strcmpi(string(T.(name)), "true") | string(T.(name)) == "1";
            end
        end
    end
end

function T = read_required_table(resultsDir, asset_id, suffix)
    path = fullfile(resultsDir, sprintf('financial_%s_%s.csv', char(asset_id), suffix));
    require_file(path, 'scripts/run_financial_decomposition.m');
    T = readtable(path, 'VariableNamingRule', 'preserve');
    if ~isdatetime(T.date)
        T.date = datetime(T.date);
    end
end

function fig = make_ssa_overview(asset, T, diag)
    L_values = unique(diag.L(diag.method == "SSA" & ~isnan(diag.L))).';
    if isempty(L_values)
        L_values = [200, 400, 600, 1000];
    end

    fig = figure('Name', sprintf('Financial SSA overview: %s', char(asset.id)), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1500 950]);
    tiledlayout(numel(L_values), 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    for i = 1:numel(L_values)
        L = L_values(i);
        ax = nexttile;
        plot(ax, T.date, T.x_z, 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8, 'DisplayName', 'z-scored level');
        hold(ax, 'on');
        plot_if_present(ax, T, sprintf('trend_L%d', L), 'trend');
        plot_if_present(ax, T, sprintf('pair_lo_L%d', L), 'pair lo');
        plot_if_present(ax, T, sprintf('pair_hi_L%d', L), 'pair hi');
        plot_if_present(ax, T, sprintf('residual_L%d', L), 'residual');
        title(ax, sprintf('%s SSA components, L = %d', char(asset.label), L));
        xlabel(ax, 'date');
        ylabel(ax, 'z-score');
        grid(ax, 'on');
        legend(ax, 'Location', 'best');
    end
    sgtitle(sprintf('Financial SSA overview: %s', char(asset.label)));
end

function fig = make_ssd_overview(asset, T, diag)
    mode_names = names_starting_with(T, 'mode_');
    ssd_diag = diag(is_ssd_method(diag.method), :);

    fig = figure('Name', sprintf('Financial trend-enabled SSD overview: %s', char(asset.id)), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1500 900]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile;
    plot(ax, T.date, T.x_z, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9, 'DisplayName', 'z-scored level');
    hold(ax, 'on');
    plot_if_present(ax, T, 'trend', 'trend');
    plot_if_present(ax, T, 'residual', 'residual');
    title(ax, sprintf('%s trend-enabled SSD trend, input, and final residual', char(asset.label)));
    xlabel(ax, 'date');
    ylabel(ax, 'z-score');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');

    ax = nexttile;
    if isempty(mode_names)
        text(ax, 0.5, 0.5, 'No trend-enabled SSD modes saved', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
    else
        plot_offset_series(ax, T.date, table_matrix(T, mode_names), mode_names);
        title(ax, 'Trend-enabled SSD extracted modes');
        xlabel(ax, 'date');
        ylabel(ax, 'offset modes');
        grid(ax, 'on');
    end

    ax = nexttile;
    L_hist = parse_number_list(first_string(ssd_diag.L_history));
    if isempty(L_hist)
        text(ax, 0.5, 0.5, 'No L history available', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
    else
        bar(ax, 1:numel(L_hist), L_hist);
        title(ax, 'Trend-enabled SSD selected window length by extraction');
        xlabel(ax, 'extraction');
        ylabel(ax, 'L');
        grid(ax, 'on');
    end
    sgtitle(sprintf('Financial trend-enabled SSD overview: %s', char(asset.label)));
end

function fig = make_tensor_overview(asset, T, reps, diag)
    fig = figure('Name', sprintf('Financial TensorCPD overview: %s', char(asset.id)), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1500 1200]);
    tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile;
    plot(ax, T.date, T.x_z, 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8, 'DisplayName', 'z-scored level');
    hold(ax, 'on');
    for i = 1:numel(reps)
        total_col = reps(i).prefix + "_total_recon";
        plot_if_present(ax, T, total_col, char(reps(i).grouping_mode + " total"));
    end
    title(ax, sprintf('%s TensorCPD total reconstructions', char(asset.label)));
    xlabel(ax, 'date');
    ylabel(ax, 'z-score');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');

    for i = 1:2
        ax = nexttile;
        if numel(reps) < i || reps(i).prefix == ""
            text(ax, 0.5, 0.5, 'No successful TensorCPD setting available', 'HorizontalAlignment', 'center');
            axis(ax, 'off');
            continue;
        end
        component_names = component_columns_for_prefix(T, reps(i).prefix);
        residual_col = reps(i).prefix + "_residual";
        if ismember(char(residual_col), T.Properties.VariableNames)
            component_names{end + 1} = char(residual_col); %#ok<AGROW>
        end
        plot_offset_series(ax, T.date, table_matrix(T, component_names), component_names);
        title(ax, sprintf('TensorCPD %s grouping: %s', char(reps(i).grouping_mode), char(reps(i).label)));
        xlabel(ax, 'date');
        ylabel(ax, 'offset components');
        grid(ax, 'on');
    end

    ax = nexttile;
    render_tensor_stability_table(ax, diag);
    sgtitle(sprintf('Financial TensorCPD overview: %s', char(asset.label)));
end

function fig = make_method_compare(asset, ssaT, ssdT, tensorT, diag, reps)
    ssa_L = choose_ssa_L(diag);
    fig = figure('Name', sprintf('Financial methods compare: %s', char(asset.id)), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1550 950]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile;
    plot(ax, ssaT.date, ssaT.x_z, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'z-scored level');
    hold(ax, 'on');
    plot_if_present(ax, ssaT, sprintf('trend_L%d', ssa_L), sprintf('SSA trend L=%d', ssa_L));
    plot_if_present(ax, ssaT, sprintf('pair_lo_L%d', ssa_L), 'SSA pair lo');
    plot_if_present(ax, ssaT, sprintf('pair_hi_L%d', ssa_L), 'SSA pair hi');
    title(ax, sprintf('%s SSA representative components', char(asset.label)));
    xlabel(ax, 'date');
    ylabel(ax, 'z-score');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');

    ax = nexttile;
    plot(ax, ssdT.date, ssdT.x_z, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'z-scored level');
    hold(ax, 'on');
    plot_if_present(ax, ssdT, 'trend', 'trend-enabled SSD trend');
    mode_names = names_starting_with(ssdT, 'mode_');
    for i = 1:numel(mode_names)
        plot(ax, ssdT.date, ssdT.(mode_names{i}), 'LineWidth', 0.9, 'DisplayName', mode_names{i});
    end
    plot_if_present(ax, ssdT, 'residual', 'trend-enabled SSD residual');
    title(ax, 'Trend-enabled SSD modes and residual');
    xlabel(ax, 'date');
    ylabel(ax, 'z-score');
    legend(ax, 'Location', 'best');
    grid(ax, 'on');

    ax = nexttile;
    plot(ax, tensorT.date, tensorT.x_z, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'z-scored level');
    hold(ax, 'on');
    if ~isempty(reps) && reps(1).prefix ~= ""
        component_names = component_columns_for_prefix(tensorT, reps(1).prefix);
        for i = 1:min(numel(component_names), 4)
            plot(ax, tensorT.date, tensorT.(component_names{i}), 'LineWidth', 0.9, 'DisplayName', component_names{i});
        end
        plot_if_present(ax, tensorT, reps(1).prefix + "_residual", char(reps(1).grouping_mode + " residual"));
    end
    title(ax, 'TensorCPD representative components');
    xlabel(ax, 'date');
    ylabel(ax, 'z-score');
    legend(ax, 'Location', 'best', 'Interpreter', 'none');
    grid(ax, 'on');

    sgtitle(sprintf('Financial method comparison: %s', char(asset.label)));
end

function fig = make_tensor_rank_contrast(asset, T, diag)
    [best_rep, stable_rep] = choose_tensor_rank_contrast_reps(diag);

    fig = figure('Name', sprintf('Financial TensorCPD rank contrast: %s', char(asset.id)), ...
        'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1550 1050]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile;
    plot(ax, T.date, T.x_z, 'Color', [0.70 0.70 0.70], 'LineWidth', 0.8, ...
        'DisplayName', 'z-scored level');
    hold(ax, 'on');
    if best_rep.prefix ~= ""
        plot_if_present(ax, T, best_rep.prefix + "_total_recon", ...
            sprintf('best fit total: %s', char(best_rep.label)));
    end
    if stable_rep.prefix ~= "" && stable_rep.prefix ~= best_rep.prefix
        plot_if_present(ax, T, stable_rep.prefix + "_total_recon", ...
            sprintf('stable total: %s', char(stable_rep.label)));
    end
    title(ax, sprintf('%s TensorCPD total reconstruction: best fit vs stable rank', char(asset.label)));
    xlabel(ax, 'date');
    ylabel(ax, 'z-score');
    legend(ax, 'Location', 'best', 'Interpreter', 'none');
    grid(ax, 'on');

    ax = nexttile;
    plot_tensor_rep_components(ax, T, stable_rep, 'Stable/recommended raw CP terms');

    ax = nexttile;
    plot_tensor_rep_components(ax, T, best_rep, 'Best-fit raw CP terms');

    sgtitle(sprintf('TensorCPD rank contrast: %s', char(asset.label)));
end

function plot_tensor_rep_components(ax, T, rep, title_prefix)
    if rep.prefix == ""
        text(ax, 0.5, 0.5, 'No TensorCPD setting available', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        return;
    end

    component_names = component_columns_for_prefix(T, rep.prefix);
    residual_col = rep.prefix + "_residual";
    if ismember(char(residual_col), T.Properties.VariableNames)
        component_names{end + 1} = char(residual_col); %#ok<AGROW>
    end
    plot_offset_series(ax, T.date, table_matrix(T, component_names), component_names);
    title(ax, sprintf('%s: %s', title_prefix, char(rep.label)), 'Interpreter', 'none');
    xlabel(ax, 'date');
    ylabel(ax, 'offset components');
    grid(ax, 'on');
end

function render_tensor_stability_table(ax, diag)
    tensor_diag = diag(diag.method == "TensorCPD" & diag.success, :);
    if isempty(tensor_diag)
        text(ax, 0.5, 0.5, 'No successful TensorCPD diagnostics available', ...
            'HorizontalAlignment', 'center');
        axis(ax, 'off');
        return;
    end

    tensor_diag = ensure_tensor_interpretability_columns(tensor_diag);
    tensor_diag = sortrows(tensor_diag, {'grouping_mode', 'cpd_rank'}, {'ascend', 'ascend'});

    lines = strings(0, 1);
    lines(end + 1) = "TensorCPD stability and recommendation table";
    lines(end + 1) = "Rec  Group  R   RelErr    Stable  LowInterp  Collapse  Cancel";
    lines(end + 1) = "---- ------ --- -------- -------- ---------- -------- --------";

    for i = 1:height(tensor_diag)
        row = tensor_diag(i, :);
        rec = marker(row.recommended_tensor_setting);
        stable = marker(~row.unstable_terms);
        low_interp = marker(row.low_interpretability);
        collapsed = marker(row.auto_group_collapsed);
        cancel = numeric_or_na(row.cancellation_ratio);
        lines(end + 1) = sprintf('%-4s %-6s %3d %8.4g %-8s %-10s %-8s %8s', ...
            char(rec), char(row.grouping_mode), row.cpd_rank, row.rel_error, ...
            char(stable), char(low_interp), char(collapsed), char(cancel)); %#ok<AGROW>
    end

    text(ax, 0.01, 0.98, strjoin(lines, newline), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Menlo', ...
        'FontSize', 9, ...
        'Interpreter', 'none');
    axis(ax, 'off');
end

function s = marker(tf)
    if tf
        s = "yes";
    else
        s = "no";
    end
end

function s = numeric_or_na(x)
    if isnan(x)
        s = "n/a";
    else
        s = string(sprintf('%.4g', x));
    end
end

function reps = choose_tensor_representatives(diag)
    modes = ["auto", "none"];
    reps = repmat(init_rep(), 1, numel(modes));
    for i = 1:numel(modes)
        subset = diag(diag.method == "TensorCPD" & diag.grouping_mode == modes(i) & diag.success, :);
        if isempty(subset)
            continue;
        end
        subset = ensure_tensor_interpretability_columns(subset);
        recommended_subset = subset(subset.recommended_tensor_setting, :);
        if ~isempty(recommended_subset)
            subset = recommended_subset;
        end

        stable_subset = subset(~subset.low_interpretability & ...
            ~subset.unstable_terms & ~subset.auto_group_collapsed, :);
        if ~isempty(stable_subset)
            subset = stable_subset;
        end

        multi_component_subset = subset(subset.num_components > 1, :);
        if ~isempty(multi_component_subset)
            subset = multi_component_subset;
        end

        subset = sortrows(subset, {'rel_error', 'fit'}, {'ascend', 'descend'});
        row = subset(1, :);
        reps(i).grouping_mode = modes(i);
        reps(i).stride_label = string(row.stride_label(1));
        reps(i).base_window_length = row.base_window_length;
        reps(i).cpd_rank = row.cpd_rank;
        reps(i).prefix = sprintf('tensor_%s_L%d_R%d_%s', ...
            char(reps(i).stride_label), reps(i).base_window_length, reps(i).cpd_rank, char(modes(i)));
        reps(i).label = tensor_rep_label(row, reps(i));
    end
end

function [best_rep, stable_rep] = choose_tensor_rank_contrast_reps(diag)
    subset = diag(diag.method == "TensorCPD" & diag.grouping_mode == "none" & diag.success, :);
    best_rep = init_rep();
    stable_rep = init_rep();
    if isempty(subset)
        return;
    end

    subset = ensure_tensor_interpretability_columns(subset);

    best_subset = sortrows(subset, {'rel_error', 'fit'}, {'ascend', 'descend'});
    best_rep = tensor_diag_row_to_rep(best_subset(1, :), "none");

    stable_subset = subset(subset.recommended_tensor_setting, :);
    if isempty(stable_subset)
        stable_subset = subset(~subset.low_interpretability & ~subset.unstable_terms, :);
    end
    if isempty(stable_subset)
        stable_subset = subset(~subset.unstable_terms, :);
    end
    if isempty(stable_subset)
        stable_subset = best_subset(1, :);
    else
        stable_subset = sortrows(stable_subset, {'rel_error', 'fit'}, {'ascend', 'descend'});
    end
    stable_rep = tensor_diag_row_to_rep(stable_subset(1, :), "none");
end

function rep = tensor_diag_row_to_rep(row, grouping_mode)
    rep = init_rep();
    rep.grouping_mode = grouping_mode;
    rep.stride_label = string(row.stride_label(1));
    rep.base_window_length = row.base_window_length;
    rep.cpd_rank = row.cpd_rank;
    rep.prefix = sprintf('tensor_%s_L%d_R%d_%s', ...
        char(rep.stride_label), rep.base_window_length, rep.cpd_rank, char(grouping_mode));
    rep.label = tensor_rep_label(row, rep);
end

function T = ensure_tensor_interpretability_columns(T)
    if ~ismember('unstable_terms', T.Properties.VariableNames)
        T.unstable_terms = false(height(T), 1);
    end
    if ~ismember('auto_group_collapsed', T.Properties.VariableNames)
        T.auto_group_collapsed = false(height(T), 1);
    end
    if ~ismember('low_interpretability', T.Properties.VariableNames)
        T.low_interpretability = T.unstable_terms | T.auto_group_collapsed;
    end
    if ~ismember('interpretability_reason', T.Properties.VariableNames)
        T.interpretability_reason = strings(height(T), 1);
    end
    if ~ismember('stable_rank_available', T.Properties.VariableNames)
        T.stable_rank_available = ~T.unstable_terms;
    end
    if ~ismember('recommended_tensor_setting', T.Properties.VariableNames)
        T.recommended_tensor_setting = false(height(T), 1);
    end
end

function label = tensor_rep_label(row, rep)
    label = sprintf('%s, L=%d, R=%d, relerr=%.4g', ...
        char(rep.stride_label), rep.base_window_length, rep.cpd_rank, row.rel_error);
    reason = string(row.interpretability_reason(1));
    if row.low_interpretability(1)
        if strlength(reason) == 0 || ismissing(reason)
            reason = "low interpretability";
        end
        label = sprintf('%s (%s)', label, char(reason));
    end
end

function rep = init_rep()
    rep = struct( ...
        'grouping_mode', "", ...
        'stride_label', "", ...
        'base_window_length', nan, ...
        'cpd_rank', nan, ...
        'prefix', "", ...
        'label', "");
end

function mask = is_ssd_method(method_values)
    method_values = string(method_values);
    mask = method_values == "SSD" | method_values == "trend-enabled SSD";
end

function L = choose_ssa_L(diag)
    subset = diag(diag.method == "SSA" & diag.success, :);
    if any(subset.L == 1000)
        L = 1000;
    elseif ~isempty(subset)
        L = subset.L(1);
    else
        L = 1000;
    end
end

function plot_if_present(ax, T, var_name, label)
    var_name = char(var_name);
    if ismember(var_name, T.Properties.VariableNames)
        plot(ax, T.date, T.(var_name), 'LineWidth', 1.0, 'DisplayName', label);
    end
end

function names = names_starting_with(T, prefix)
    all_names = T.Properties.VariableNames;
    names = all_names(startsWith(all_names, prefix));
end

function names = component_columns_for_prefix(T, prefix)
    all_names = T.Properties.VariableNames;
    prefix = char(prefix);
    mask = startsWith(all_names, [prefix '_c']);
    names = all_names(mask);
    names = sort_component_names(names);
end

function names = sort_component_names(names)
    if isempty(names)
        return;
    end
    idx = zeros(size(names));
    for i = 1:numel(names)
        tok = regexp(names{i}, '_c(\d+)$', 'tokens', 'once');
        if isempty(tok)
            idx(i) = inf;
        else
            idx(i) = str2double(tok{1});
        end
    end
    [~, order] = sort(idx);
    names = names(order);
end

function X = table_matrix(T, names)
    X = zeros(height(T), numel(names));
    for i = 1:numel(names)
        X(:, i) = T.(names{i});
    end
end

function plot_offset_series(ax, t, Y, names)
    Y = Y(:, :);
    if isempty(Y)
        text(ax, 0.5, 0.5, 'No components available', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        return;
    end

    colors = lines(size(Y, 2));
    amp = max(abs(Y(:)), [], 'omitnan');
    if isnan(amp) || amp <= eps
        amp = 1;
    end
    spacing = 2.5 * amp;
    hold(ax, 'on');
    for k = 1:size(Y, 2)
        offset = (size(Y, 2) - k) * spacing;
        plot(ax, t, Y(:, k) + offset, 'LineWidth', 1.0, 'Color', colors(k, :));
        text(ax, t(1), offset, ['  ' names{k}], 'Color', colors(k, :), ...
            'VerticalAlignment', 'bottom', 'Interpreter', 'none');
    end
    yticks(ax, []);
end

function nums = parse_number_list(s)
    s = string(s);
    if strlength(s) == 0 || ismissing(s)
        nums = [];
        return;
    end
    parts = split(s, ';');
    nums = str2double(parts).';
    nums = nums(~isnan(nums));
end

function s = first_string(values)
    if isempty(values)
        s = "";
        return;
    end
    if iscell(values)
        s = string(values{1});
    else
        s = string(values(1));
    end
end

function require_file(path_str, producer)
    if exist(path_str, 'file') ~= 2
        error('run_financial_visualizations:MissingFile', ...
            'Missing %s. Run %s first.', char(path_str), char(producer));
    end
end
