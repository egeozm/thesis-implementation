%RUN_TENSOR_DIAGNOSTIC Deep inspection of the pure-CPD tensor reconstruction.
%
%   Runs one noiseless multiscale case, prints per-stride CPD reconstruction
%   diagnostics, and plots the reconstructed components against the ground truth.

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
config.signal_type = 'multiscale';
config.N = 2000;
config.snr_db = Inf;
config.seed = 11;
config.fs = 1;
config.params = tensor_merge_params(struct( ...
    'strides', [1, 2], ...
    'base_window_length', 200, ...
    'cpd_rank', 4, ...
    'grouping_mode', 'oracle', ...
    'reconstruction_mode', 'cpd', ...
    'cpd_method', 'cpd', ...
    'cpd_options', struct()));

[x, components, t] = generate_signal(config.signal_type, config.N, config.snr_db, config.seed);
outcome = evaluate_tensor_cpd(config.signal_type, x, components, config.params, config.fs);
if ~outcome.success
    error('run_tensor_diagnostic:Failed', 'Tensor diagnostic failed: %s', outcome.failure_reason);
end

true_mat = truth_matrix(config.signal_type, components);
signal_names = outcome.component_names;
recon = outcome.reconstruction;
stride_table = build_stride_table(recon.per_stride, true_mat, outcome.est_components(:, 1:size(true_mat, 2)));

fprintf('\n=== Tensor diagnostic ===\n');
fprintf('Signal              : %s\n', config.signal_type);
fprintf('N                   : %d\n', config.N);
fprintf('SNR                 : %s\n', snr_label(config.snr_db));
fprintf('Seed                : %d\n', config.seed);
fprintf('Grouping            : %s\n', config.params.grouping_mode);
fprintf('Reconstruction mode : %s\n', config.params.reconstruction_mode);
fprintf('CP rank             : %d\n', config.params.cpd_rank);
fprintf('CP fit              : %.6f\n', outcome.fit);
fprintf('CP rel. error       : %.6f\n', outcome.rel_error);
fprintf('Mean NMSE           : %.6e\n', outcome.mean_nmse);
fprintf('Mean |rho|          : %.6f\n', outcome.mean_abs_rho);
fprintf('Mean peak error     : %.6f\n', outcome.mean_peak_error);
fprintf('Successful strides  : %d / %d\n', ...
    outcome.reconstruction.grouping.num_successful_strides, numel(recon.per_stride));

disp('Per-stride diagnostic summary:');
disp(stride_table);
print_metrics_table(outcome.metrics, signal_names);

figure('Name', 'Tensor diagnostic: reconstructed components', 'NumberTitle', 'off');
sgtitle('Tensor pure-CPD reconstruction vs truth (noiseless multiscale)');
num_components = size(true_mat, 2);
for g = 1:num_components
    subplot(num_components, 1, g);
    plot(t, true_mat(:, g), 'k-', 'LineWidth', 1); hold on;
    plot(t, outcome.est_components(:, g), 'r--', 'LineWidth', 1);
    title(sprintf('%s component', signal_names{g}));
    legend('truth', 'tensor', 'Location', 'best');
    xlabel('n');
    ylabel('amplitude');
    grid on;
end

figure('Name', 'Tensor diagnostic: per-stride total reconstructions', 'NumberTitle', 'off');
num_stride_plots = numel(recon.per_stride);
for i = 1:num_stride_plots
    subplot(num_stride_plots, 1, i);
    if ~recon.per_stride(i).success
        title(sprintf('stride %d failed: %s', recon.per_stride(i).stride, recon.per_stride(i).failure_reason));
        axis off;
        continue;
    end
    t_dec = recon.per_stride(i).sample_indices;
    plot(t_dec, x(t_dec), 'k-', 'LineWidth', 1); hold on;
    plot(t_dec, recon.per_stride(i).total_decimated_reconstruction, 'r--', 'LineWidth', 1);
    title(sprintf('stride %d total reconstruction', recon.per_stride(i).stride));
    legend('truth samples', 'CPD recon', 'Location', 'best');
    xlabel('n');
    ylabel('amplitude');
    grid on;
end

resultFile = fullfile(resultsDir, 'tensor_diagnostic_results.mat');
save(resultFile, 'config', 'x', 'components', 'outcome', 'stride_table');
fprintf('Saved               : %s\n', resultFile);

function true_mat = truth_matrix(signal_type, components)
    switch lower(string(signal_type))
        case "multiscale"
            true_mat = [components.trend, components.slow, components.fast];
        case "close"
            true_mat = [components.comp1, components.comp2];
        otherwise
            error('run_tensor_diagnostic:UnknownSignal', 'Unknown signal type %s.', signal_type);
    end
end

function tbl = build_stride_table(stride_results, true_mat, final_est)
    rows = repmat(init_stride_row(), numel(stride_results), 1);
    for i = 1:numel(stride_results)
        rows(i).stride = stride_results(i).stride;
        rows(i).success = stride_results(i).success;
        rows(i).failure_reason = string(stride_results(i).failure_reason);
        rows(i).num_blocks = stride_results(i).num_blocks;
        rows(i).signal_rel_error = stride_results(i).signal_rel_error;
        rows(i).combine_weight = stride_results(i).combine_weight;
        rows(i).coverage_fraction = stride_results(i).coverage_fraction;
        if stride_results(i).success
            truth_dec = true_mat(stride_results(i).sample_indices, :);
            est_dec = final_est(stride_results(i).sample_indices, :);
            num_components = min(size(truth_dec, 2), size(est_dec, 2));
            rho = nan(1, 3);
            for g = 1:num_components
                rho(g) = pearson_raw(truth_dec(:, g), est_dec(:, g));
            end
            rows(i).comp1_rho = rho(1);
            rows(i).comp2_rho = rho(2);
            rows(i).comp3_rho = rho(3);
        end
    end
    tbl = struct2table(rows);
end

function row = init_stride_row()
    row = struct( ...
        'stride', nan, ...
        'success', false, ...
        'failure_reason', "", ...
        'num_blocks', nan, ...
        'signal_rel_error', nan, ...
        'combine_weight', nan, ...
        'coverage_fraction', nan, ...
        'comp1_rho', nan, ...
        'comp2_rho', nan, ...
        'comp3_rho', nan);
end

function print_metrics_table(m, names)
    fprintf('\n--- Final component metrics ---\n');
    fprintf('%-8s %12s %12s %12s\n', 'comp', 'NMSE', 'Pearson', 'peak|err|');
    for k = 1:numel(names)
        fprintf('%-8s %12.4e %12.4f %12.6f\n', ...
            names{k}, m.nmse(k), m.pearson(k), m.peak_freq_error(k));
    end
end

function label = snr_label(v)
    if isinf(v)
        label = 'Inf';
    else
        label = num2str(v);
    end
end

function c = pearson_raw(a, b)
    a = a(:) - mean(a);
    b = b(:) - mean(b);
    den = norm(a) * norm(b);
    if den == 0
        c = 0;
    else
        c = (a.' * b) / den;
    end
end
