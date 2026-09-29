function outcome = evaluate_ssd(signal_type, x, components, params, fs)
%EVALUATE_SSD Evaluate SSD decomposition against synthetic ground truth.
%
%   outcome = evaluate_ssd(signal_type, x, components, params, fs)
%
%   params: ssd_merge_params struct (frozen grouping hyperparameters + SSD fields).
%   Builds est_components = [modes, residual] for matching/metrics, or
%   [trend, modes, residual] when params.extract_trend_first is enabled.
%
%   Interpretation note:
%   - In the multiscale case, the final residual is treated as a trend-dominated
%     remainder for evaluation, not as a guaranteed pure trend estimate. It may still
%     contain residual noise or unmodeled low-frequency content.
%   - The success flag is only a minimal implementation check (enough oscillatory
%     modes extracted); actual separation quality is determined by the reported
%     component metrics.

    if nargin < 5 || isempty(fs)
        fs = 1;
    end

    signal_type = string(lower(string(signal_type)));
    params = ssd_merge_params(params);
    params.fs = fs;

    dec = ssd_decompose(x, params);
    [true_components, component_names] = truth_components(signal_type, components);
    [success, failure_reason] = check_success(signal_type, dec);

    outcome = struct();
    outcome.signal_type = signal_type;
    outcome.decomposition = dec;
    outcome.component_names = component_names;
    outcome.success = success;
    outcome.failure_reason = string(failure_reason);
    outcome.num_ssd_modes = dec.num_iterations;

    if ~success
        outcome.metrics = empty_metrics(size(true_components, 2));
        outcome.est_components = [];
        outcome.mean_nmse = nan;
        outcome.mean_abs_rho = nan;
        outcome.mean_peak_error = nan;
        outcome.score = 1e3;
        return;
    end

    if isfield(dec, 'trend_extracted') && dec.trend_extracted
        est_components = [dec.trend, dec.modes, dec.residual];
    else
        est_components = [dec.modes, dec.residual];
    end
    [assignment, ~] = match_components(true_components, est_components);
    metrics = compute_metrics(true_components, est_components, assignment, fs);

    outcome.metrics = metrics;
    outcome.est_components = est_components;
    outcome.mean_nmse = mean(metrics.nmse, 'omitnan');
    outcome.mean_abs_rho = mean(abs(metrics.pearson), 'omitnan');
    outcome.mean_peak_error = mean(metrics.peak_freq_error, 'omitnan');
    outcome.score = outcome.mean_nmse;
end

function [true_components, component_names] = truth_components(signal_type, components)
    switch signal_type
        case "multiscale"
            true_components = [components.trend, components.slow, components.fast];
            component_names = {"trend", "slow", "fast"};
        case "close"
            true_components = [components.comp1, components.comp2];
            component_names = {"comp1", "comp2"};
        otherwise
            error('evaluate_ssd:UnknownSignal', 'Unknown signal type %s', signal_type);
    end
end

function [success, failure_reason] = check_success(signal_type, dec)
    M = dec.num_iterations;
    switch signal_type
        case "multiscale"
            if M < 2
                success = false;
                failure_reason = "fewer than two SSD oscillatory modes; inspect component metrics separately";
            else
                success = true;
                failure_reason = "";
            end
        case "close"
            if M < 2
                success = false;
                failure_reason = "fewer than two SSD oscillatory modes; inspect component metrics separately";
            else
                success = true;
                failure_reason = "";
            end
        otherwise
            error('evaluate_ssd:UnknownSignal', 'Unknown signal type %s', signal_type);
    end
end

function metrics = empty_metrics(num_components)
    metrics = struct();
    metrics.nmse = nan(num_components, 1);
    metrics.pearson = nan(num_components, 1);
    metrics.peak_freq_error = nan(num_components, 1);
    metrics.leakage = nan(num_components, 1);
    metrics.dominant_freq_true = nan(num_components, 1);
    metrics.dominant_freq_est = nan(num_components, 1);
    metrics.true_idx = (1:num_components).';
    metrics.est_idx = nan(num_components, 1);
end
