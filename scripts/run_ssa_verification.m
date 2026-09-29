%RUN_SSA_VERIFICATION End-to-end SSA verification (thesis Sec. 3.1–3.2, 3.6).
%
%   From project root:
%     run('scripts/run_ssa_verification.m')
%
%   Or from scripts folder with project root on path (adds ../ssa automatically).

thisDir = fileparts(mfilename('fullpath'));
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'ssa'));

%% Parameters (thesis)
N = 2000;
seed = 42;
L_ref = 1000;
L_sweep = [200, 400, 600, 1000];
snr_list = [Inf, 20];
fs = 1; % normalized time → FFT peaks in cycles/sample
signal_names = {'trend', 'slow', 'fast'};
all_group_names = {'trend', 'slow', 'fast', 'residual'};

%% --- Noiseless reference: scree + reconstructions (L = 1000) ---
[x0, comp0, t] = generate_signal('multiscale', N, Inf, seed);
true_mat0 = [comp0.trend, comp0.slow, comp0.fast];
groups0 = infer_multiscale_groups(x0, L_ref, true_mat0);

[rec0, sv0, ~, ~, ~] = ssa_decompose(x0, L_ref, groups0);
m0 = compute_metrics(true_mat0, rec0(:, 1:3), (1:3).', fs);

figure('Name', 'SSA singular values (noiseless, L=1000)', 'NumberTitle', 'off');
num_sv_plot = min(20, numel(sv0));
semilogy(1:num_sv_plot, sv0(1:num_sv_plot), 'o-', 'LineWidth', 1);
xlabel('Index r');
ylabel('Singular value \sigma_r');
title('Scree plot (first 20) — multiscale synthetic, SNR = \infty');
grid on;

figure('Name', 'SSA reconstructions vs ground truth', 'NumberTitle', 'off');
sgtitle(sprintf('Noiseless multiscale signal, L = %d (truth-guided grouping)', L_ref));
for g = 1:4
    subplot(2, 2, g);
    if g <= 3
        plot(t, true_mat0(:, g), 'k-', 'LineWidth', 1); hold on;
        plot(t, rec0(:, g), 'r--', 'LineWidth', 1);
        legend('true', 'SSA group', 'Location', 'best');
    else
        plot(t, rec0(:, g), 'b-', 'LineWidth', 1);
        yline(0, 'k:');
        legend('residual', 'zero line', 'Location', 'best');
    end
    xlabel('n');
    ylabel('amplitude');
    title(sprintf('%s group', all_group_names{g}));
    grid on;
end
report_groups(groups0, 'Noiseless, L=1000');
print_metrics_table(m0, signal_names, 'Noiseless signal components, L=1000');

fprintf('\nSum of reconstructions vs signal: rel L2 err = %.3e\n', ...
    norm(x0 - sum(rec0, 2)) / norm(x0));

%% --- Window length and SNR sensitivity ---
num_signal_components = numel(signal_names);
nmse_L = nan(numel(L_sweep), numel(snr_list), num_signal_components);
corr_L = nan(numel(L_sweep), numel(snr_list), num_signal_components);

for si = 1:numel(snr_list)
    snr_db = snr_list(si);
    [x_s, comp_s, ~] = generate_signal('multiscale', N, snr_db, seed);
    true_m = [comp_s.trend, comp_s.slow, comp_s.fast];

    for Li = 1:numel(L_sweep)
        L = L_sweep(Li);
        groups = infer_multiscale_groups(x_s, L, true_m);
        rec = ssa_decompose(x_s, L, groups);
        m = compute_metrics(true_m, rec(:, 1:3), (1:3).', fs);
        for k = 1:num_signal_components
            nmse_L(Li, si, k) = m.nmse(k);
            corr_L(Li, si, k) = m.pearson(k);
        end
    end
end

figure('Name', 'NMSE vs L (by component)', 'NumberTitle', 'off');
for k = 1:num_signal_components
    subplot(1, 3, k);
    hold on;
    for si = 1:numel(snr_list)
        lbl = ternary(isinf(snr_list(si)), 'SNR=\infty', sprintf('SNR=%g dB', snr_list(si)));
        plot(L_sweep, nmse_L(:, si, k), 'o-', 'LineWidth', 1.2, 'DisplayName', lbl);
    end
    set(gca, 'YScale', 'log');
    xlabel('L');
    ylabel('NMSE');
    title(signal_names{k});
    legend('Location', 'best');
    grid on;
end
sgtitle('Normalized MSE vs window length (signal components only)');

figure('Name', 'Pearson correlation vs L', 'NumberTitle', 'off');
for k = 1:num_signal_components
    subplot(1, 3, k);
    hold on;
    for si = 1:numel(snr_list)
        lbl = ternary(isinf(snr_list(si)), 'SNR=\infty', sprintf('SNR=%g dB', snr_list(si)));
        plot(L_sweep, corr_L(:, si, k), 'o-', 'LineWidth', 1.2, 'DisplayName', lbl);
    end
    xlabel('L');
    ylabel('Pearson \rho');
    ylim([-0.1, 1.05]);
    title(signal_names{k});
    legend('Location', 'best');
    grid on;
end
sgtitle('Pearson correlation vs window length (signal components only)');

fprintf('\n=== Sensitivity summary (mean NMSE over 3 signal components) ===\n');
for si = 1:numel(snr_list)
    for Li = 1:numel(L_sweep)
        fprintf('  SNR %s, L=%4d: mean NMSE = %.4e, mean |rho| = %.4f\n', ...
            ternary(isinf(snr_list(si)), 'Inf', sprintf('%g', snr_list(si))), ...
            L_sweep(Li), ...
            mean(squeeze(nmse_L(Li, si, :)), 'omitnan'), ...
            mean(abs(squeeze(corr_L(Li, si, :))), 'omitnan'));
    end
end

fprintf('\nrun_ssa_verification finished.\n');

%% --- Local helpers ---
function groups = infer_multiscale_groups(x, L, true_components)
    % For the synthetic verification signal, use the known ground truth to guide
    % grouping of the leading eigentriples. This avoids a brittle hard-coded split.
    x = x(:);
    N_local = numel(x);
    H = embed_hankel(x, L);
    [U, S, V] = svd(H, 'econ');
    R = size(S, 1);
    inspect_r = min(12, R);

    elementary = zeros(N_local, inspect_r);
    corr_to_truth = zeros(3, inspect_r);
    for r = 1:inspect_r
        Hr = S(r, r) * (U(:, r) * V(:, r).');
        elementary(:, r) = diagonal_average(Hr, N_local);
        for c = 1:3
            corr_to_truth(c, r) = abs(pearson_corr(elementary(:, r), true_components(:, c)));
        end
    end

    used = false(1, inspect_r);
    slow_idx = select_oscillatory_pair(corr_to_truth(2, :), used);
    used(slow_idx) = true;

    fast_idx = select_oscillatory_pair(corr_to_truth(3, :), used);
    used(fast_idx) = true;

    trend_idx = select_trend_group(corr_to_truth(1, :), used);
    used(trend_idx) = true;

    signal_idx = unique([trend_idx(:); slow_idx(:); fast_idx(:)]).';
    residual_idx = setdiff(1:R, signal_idx);

    groups = {sort(trend_idx), sort(slow_idx), sort(fast_idx), residual_idx};
end

function idx = select_oscillatory_pair(scores, used)
    available = find(~used);
    if isempty(available)
        idx = [];
        return;
    end

    [~, best_pos] = max(scores(available));
    seed = available(best_pos);
    remaining = available(available ~= seed);
    neighbors = remaining(abs(remaining - seed) == 1);

    if ~isempty(neighbors)
        [~, partner_pos] = max(scores(neighbors));
        partner = neighbors(partner_pos);
    elseif ~isempty(remaining)
        [~, partner_pos] = max(scores(remaining));
        partner = remaining(partner_pos);
    else
        partner = [];
    end

    idx = unique(sort([seed, partner]));
end

function idx = select_trend_group(scores, used)
    available = find(~used);
    if isempty(available)
        idx = [];
        return;
    end

    available_scores = scores(available);
    [sorted_scores, order] = sort(available_scores, 'descend');
    if isempty(sorted_scores)
        idx = [];
        return;
    end

    threshold = max(0.15 * sorted_scores(1), 0.05);
    keep_mask = sorted_scores >= threshold;
    keep = available(order(keep_mask));
    if isempty(keep)
        keep = available(order(1));
    end

    idx = sort(keep(1:min(4, numel(keep))));
end

function report_groups(groups, label_str)
    fprintf('\n--- %s grouping ---\n', label_str);
    fprintf('  trend   : %s\n', mat2str(groups{1}));
    fprintf('  slow    : %s\n', mat2str(groups{2}));
    fprintf('  fast    : %s\n', mat2str(groups{3}));
    if numel(groups{4}) <= 20
        residual_repr = mat2str(groups{4});
    else
        residual_repr = sprintf('[%d:%d] (%d components)', ...
            groups{4}(1), groups{4}(end), numel(groups{4}));
    end
    fprintf('  residual: %s\n', residual_repr);
end

function print_metrics_table(m, names, title_str)
    fprintf('\n--- %s ---\n', title_str);
    fprintf('%-8s %10s %10s %14s %10s\n', 'comp', 'NMSE', 'Pearson', 'peak|err|', 'leak(est)');
    Nt = numel(names);
    Ne = numel(m.leakage);
    for k = 1:Nt
        j = m.est_idx(k);
        lk = nan;
        if ~isnan(j) && j >= 1 && j <= Ne
            lk = m.leakage(j);
        end
        fprintf('%-8s %10.4e %10.4f %14.6f %10.4f\n', ...
            names{k}, m.nmse(k), m.pearson(k), m.peak_freq_error(k), lk);
    end
end

function c = pearson_corr(a, b)
    a = a(:) - mean(a);
    b = b(:) - mean(b);
    den = norm(a) * norm(b);
    if den == 0
        c = 0;
    else
        c = (a.' * b) / den;
    end
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
