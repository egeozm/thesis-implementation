%RUN_TENSOR_VERIFICATION Smoke test for decimated Hankel tensor + CPD.

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
config.snr_db = 20;
config.seed = 11;

params = tensor_merge_params(struct( ...
    'strides', [1, 2, 4, 8], ...
    'base_window_length', 400, ...
    'cpd_rank', 3, ...
    'cpd_method', 'cpd', ...
    'cpd_options', struct()));

[x, components] = generate_signal(config.signal_type, config.N, config.snr_db, config.seed);
tensor_out = tensor_build_decimated_hankel(x, params);
cpd_out = tensor_run_cpd(tensor_out, params);
if cpd_out.success
    recon_out = tensor_reconstruct_components(cpd_out, x);
else
    recon_out = struct('success', false, 'failure_reason', cpd_out.failure_reason, 'support_length', nan);
end

fprintf('\n=== Tensor verification ===\n');
fprintf('Signal        : %s\n', config.signal_type);
fprintf('Tensor size   : [%s]\n', numvec_to_str(tensor_out.tensor_size));
fprintf('Valid strides : %s\n', numvec_to_str(tensor_out.valid_strides));
fprintf('Num blocks    : %d\n', tensor_out.num_blocks);
fprintf('Coverage      : %.2f%%\n', 100 * tensor_out.coverage.coverage_fraction);
fprintf('Tensor norm   : %.6g\n', tensor_out.tensor_norm);

if cpd_out.success
    fprintf('CPD method    : %s\n', cpd_out.method);
    fprintf('CP rank       : %d\n', cpd_out.rank);
    fprintf('Fit           : %.6f\n', cpd_out.fit);
    fprintf('Rel. error    : %.6f\n', cpd_out.rel_error);
else
    fprintf('CPD skipped   : %s\n', cpd_out.failure_reason);
end
if recon_out.success
    fprintf('Recon length  : %d\n', recon_out.support_length);
else
    fprintf('Recon skipped : %s\n', recon_out.failure_reason);
end

disp('Stride summary:');
disp(stride_table(tensor_out.stride_info));

resultFile = fullfile(resultsDir, 'tensor_verification_results.mat');
save(resultFile, 'config', 'params', 'x', 'components', 'tensor_out', 'cpd_out', 'recon_out');
fprintf('Saved         : %s\n', resultFile);

function s = numvec_to_str(v)
    s = strjoin(arrayfun(@num2str, v, 'UniformOutput', false), ', ');
end

function tbl = stride_table(stride_info)
    num_slices = numel(stride_info);
    rows = repmat(struct( ...
        'stride', nan, ...
        'decimated_length', nan, ...
        'raw_rows', nan, ...
        'raw_cols', nan, ...
        'block_rows', nan, ...
        'block_cols', nan, ...
        'num_row_blocks', nan, ...
        'num_col_blocks', nan, ...
        'num_blocks', nan), num_slices, 1);

    for i = 1:num_slices
        rows(i).stride = stride_info(i).stride;
        rows(i).decimated_length = stride_info(i).decimated_length;
        rows(i).raw_rows = stride_info(i).raw_rows;
        rows(i).raw_cols = stride_info(i).raw_cols;
        rows(i).block_rows = stride_info(i).block_num_rows;
        rows(i).block_cols = stride_info(i).block_num_cols;
        rows(i).num_row_blocks = stride_info(i).num_row_blocks;
        rows(i).num_col_blocks = stride_info(i).num_col_blocks;
        rows(i).num_blocks = stride_info(i).num_blocks;
    end

    tbl = struct2table(rows);
end
