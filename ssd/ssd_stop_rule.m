function stop = ssd_stop_rule(mode_signal, residual_norm_before, iter, score, first_score, prev_score, params)
%SSD_STOP_RULE True if SSD should stop after extracting mode_signal.
%
%   stop = ssd_stop_rule(mode_signal, residual_norm_before, iter, score, ...
%       first_score, prev_score, params)

    if iter > params.max_modes
        stop = true;
        return;
    end

    nr = max(residual_norm_before, eps);
    if norm(mode_signal) / nr < params.min_mode_energy_frac
        stop = true;
        return;
    end

    if iter > params.min_modes_before_score_stop
        if ~isnan(first_score) && score < params.min_score_ratio * first_score
            stop = true;
            return;
        end
        if ~isnan(prev_score) && score < params.min_score_drop_ratio * prev_score
            stop = true;
            return;
        end
    end

    stop = false;
end
