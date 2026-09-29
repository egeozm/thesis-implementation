function [best_params, best_row, eligible_rows] = select_ssa_params(summary_table, margin)
%SELECT_SSA_PARAMS Select a robust fixed parameter set from calibration summary.
%
%   [best_params, best_row, eligible_rows] = select_ssa_params(summary_table, margin)

    if nargin < 2 || isempty(margin)
        margin = 0.05;
    end
    if isempty(summary_table)
        error('select_ssa_params:EmptyTable', 'summary_table must not be empty.');
    end

    min_score = min(summary_table.score);
    eligible_rows = summary_table(summary_table.score <= (1 + margin) * min_score, :);
    eligible_rows = sortrows(eligible_rows, ...
        {'failure_rate', 'mean_nmse', 'r_max', 'freq_tol', 'min_osc_freq', 'sv_ratio_thresh'}, ...
        {'ascend', 'ascend', 'ascend', 'ascend', 'ascend', 'descend'});

    best_row = eligible_rows(1, :);
    best_params = struct(...
        'r_max', best_row.r_max, ...
        'sv_ratio_thresh', best_row.sv_ratio_thresh, ...
        'freq_tol', best_row.freq_tol, ...
        'min_osc_freq', best_row.min_osc_freq);
end
