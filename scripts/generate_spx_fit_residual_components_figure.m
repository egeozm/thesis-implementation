%GENERATE_SPX_FIT_RESIDUAL_COMPONENTS_FIGURE Generate clearer RQ3 figure.
%
%   Exports:
%     results/financial_spx_fit_residual_components.png

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
resultsDir = fullfile(projRoot, 'results');

ssa = readtable(fullfile(resultsDir, 'financial_spx_ssa.csv'));
ssd = readtable(fullfile(resultsDir, 'financial_spx_ssd.csv'));
tensor = readtable(fullfile(resultsDir, 'financial_spx_tensor.csv'));

if ~isdatetime(ssa.date)
    ssa.date = datetime(ssa.date);
end
if ~isdatetime(ssd.date)
    ssd.date = datetime(ssd.date);
end
if ~isdatetime(tensor.date)
    tensor.date = datetime(tensor.date);
end

ssa_recon = ssa.trend_L400 + ssa.pair_lo_L400 + ssa.pair_hi_L400;
ssa_residual = ssa.x_z - ssa_recon;

mode_names = ssd.Properties.VariableNames(startsWith(ssd.Properties.VariableNames, 'mode_'));
ssd_recon = ssd.trend;
for i = 1:numel(mode_names)
    ssd_recon = ssd_recon + ssd.(mode_names{i});
end
ssd_residual = ssd.x_z - ssd_recon;

tensor_prefix = 'tensor_S1_S2_L200_R6_none';
tensor_recon = tensor.([tensor_prefix '_total_recon']);
tensor_residual = tensor.([tensor_prefix '_residual']);
tensor_component_names = component_columns_for_prefix(tensor, tensor_prefix);
tensor_components = table_matrix(tensor, tensor_component_names);

fig = figure('Name', 'S&P 500 fit residual components', ...
    'NumberTitle', 'off', 'Color', 'w', 'Position', [80 80 1450 1250]);
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'loose');

ax = nexttile;
plot(ax, ssa.date, ssa.x_z, 'Color', [0.72 0.72 0.72], ...
    'LineWidth', 1.0, 'DisplayName', 'observed');
hold(ax, 'on');
plot(ax, ssa.date, ssa_recon, 'LineWidth', 1.3, ...
    'DisplayName', 'SSA');
plot(ax, ssd.date, ssd_recon, 'LineWidth', 1.3, ...
    'DisplayName', 'SSD');
plot(ax, tensor.date, tensor_recon, 'LineWidth', 1.3, ...
    'DisplayName', 'TensorCPD R=6');
title(ax, 'Total reconstruction: all methods capture the long-term market movement');
format_axis(ax, 'z-score');
format_legend(legend(ax, 'Location', 'northwest', 'FontSize', 15, 'NumColumns', 2));

ax = nexttile;
plot(ax, ssa.date, ssa_residual, 'LineWidth', 1.0, ...
    'DisplayName', 'SSA residual');
hold(ax, 'on');
plot(ax, ssd.date, ssd_residual, 'LineWidth', 1.0, ...
    'DisplayName', 'SSD residual');
plot(ax, tensor.date, tensor_residual, 'LineWidth', 1.0, ...
    'DisplayName', 'TensorCPD residual (R=6)');
yline(ax, 0, 'k:', 'LineWidth', 0.8, 'DisplayName', 'zero');
title(ax, 'Residual comparison: TensorCPD has the lowest selected relative error');
format_axis(ax, 'residual');
format_legend(legend(ax, 'Location', 'northwest', 'FontSize', 15, 'NumColumns', 2));

ax = nexttile;
plot_offset_components(ax, tensor.date, tensor_components, tensor_component_names);
title(ax, 'TensorCPD R=6 components shown with offsets for readability');
xlabel(ax, 'date');
ylabel(ax, 'offset components');
grid(ax, 'on');
set(ax, 'FontSize', 16, 'LineWidth', 1.1, ...
    'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.75 0.75 0.75], 'GridAlpha', 0.45);
ax.Title.FontSize = 17;
ax.Title.FontWeight = 'bold';
ax.Title.Color = 'k';
ax.XLabel.FontSize = 16;
ax.XLabel.Color = 'k';
ax.YLabel.FontSize = 16;
ax.YLabel.Color = 'k';

fig_title = sgtitle('S&P 500 financial decomposition: fit, residual, and TensorCPD components', ...
    'FontSize', 20, 'FontWeight', 'bold');
fig_title.Color = 'k';

axes_handles = findall(fig, 'Type', 'axes');
for i = 1:numel(axes_handles)
    if isprop(axes_handles(i), 'Toolbar') && ~isempty(axes_handles(i).Toolbar)
        axes_handles(i).Toolbar.Visible = 'off';
    end
end

outFile = fullfile(resultsDir, 'financial_spx_fit_residual_components.png');
exportgraphics(fig, outFile, 'Resolution', 300);
close(fig);
fprintf('Saved %s\n', outFile);

function format_axis(ax, y_label)
    xlabel(ax, 'date');
    ylabel(ax, y_label);
    grid(ax, 'on');
    set(ax, 'FontSize', 16, 'LineWidth', 1.1, ...
        'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
        'GridColor', [0.75 0.75 0.75], 'GridAlpha', 0.45);
    ax.Title.FontSize = 17;
    ax.Title.FontWeight = 'bold';
    ax.Title.Color = 'k';
    ax.XLabel.FontSize = 16;
    ax.XLabel.Color = 'k';
    ax.YLabel.FontSize = 16;
    ax.YLabel.Color = 'k';
end

function format_legend(lgd)
    lgd.Color = 'w';
    lgd.TextColor = 'k';
    lgd.EdgeColor = [0.8 0.8 0.8];
end

function names = component_columns_for_prefix(T, prefix)
    all_names = T.Properties.VariableNames;
    prefix = char(prefix);
    mask = startsWith(all_names, [prefix '_c']);
    names = all_names(mask);
    names = sort_component_names(names);
end

function names = sort_component_names(names)
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

function plot_offset_components(ax, dates, components, names)
    hold(ax, 'on');
    if isempty(components)
        text(ax, 0.5, 0.5, 'No TensorCPD components available');
        return;
    end

    amp = max(abs(components(:)), [], 'omitnan');
    if amp <= eps || isnan(amp)
        amp = 1;
    end
    spacing = 1.4 * amp;
    colors = lines(size(components, 2));
    for i = 1:size(components, 2)
        offset = (size(components, 2) - i) * spacing;
        plot(ax, dates, components(:, i) + offset, ...
            'LineWidth', 1.1, 'Color', colors(i, :));
        short_name = regexp(names{i}, '_c(\d+)$', 'tokens', 'once');
        if isempty(short_name)
            label = names{i};
        else
            label = sprintf('CP term %s', short_name{1});
        end
        text(ax, dates(1), offset, ['  ' label], ...
            'Color', colors(i, :), 'FontWeight', 'bold', ...
            'FontSize', 15, 'VerticalAlignment', 'bottom', ...
            'Interpreter', 'none');
    end
    yticks(ax, []);
end
