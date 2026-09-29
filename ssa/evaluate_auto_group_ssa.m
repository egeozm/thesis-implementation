function outcome = evaluate_auto_group_ssa(signal_type, features, components, params, fs)
%EVALUATE_AUTO_GROUP_SSA Evaluate automatic SSA grouping against synthetic truth.
%
%   outcome = evaluate_auto_group_ssa(signal_type, features, components, params, fs)

    if nargin < 5 || isempty(fs)
        fs = 1;
    end

    signal_type = string(lower(string(signal_type)));
    grouping = auto_group_ssa(features, params);
    [true_components, component_names] = truth_components(signal_type, components);
    [component_groups, success, failure_reason] = groups_for_signal(signal_type, grouping);

    outcome = struct();
    outcome.signal_type = signal_type;
    outcome.grouping = grouping;
    outcome.component_names = component_names;
    outcome.success = success;
    outcome.failure_reason = string(failure_reason);
    outcome.num_pairs = numel(grouping.oscillatory_pairs);
    outcome.trend_group_size = numel(grouping.trend_idx);

    if ~success
        outcome.metrics = empty_metrics(size(true_components, 2));
        outcome.est_components = [];
        outcome.mean_nmse = nan;
        outcome.mean_abs_rho = nan;
        outcome.mean_peak_error = nan;
        outcome.score = 1e3;
        return;
    end

    est_components = ssa_reconstruct_groups(features, component_groups);
    assignment = (1:size(true_components, 2)).';
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
            error('evaluate_auto_group_ssa:UnknownSignal', 'Unknown signal type %s', signal_type);
    end
end

function [component_groups, success, failure_reason] = groups_for_signal(signal_type, grouping)
    pairs = grouping.oscillatory_pairs;
    switch signal_type
        case "multiscale"
            if numel(pairs) < 2
                component_groups = {};
                success = false;
                failure_reason = "fewer than two oscillatory pairs";
                return;
            end
            if isempty(grouping.trend_idx)
                component_groups = {};
                success = false;
                failure_reason = "empty trend group";
                return;
            end
            component_groups = {grouping.trend_idx, pairs{1}, pairs{2}};
            success = true;
            failure_reason = "";
        case "close"
            if numel(pairs) < 2
                component_groups = {};
                success = false;
                failure_reason = "fewer than two oscillatory pairs";
                return;
            end
            component_groups = {pairs{1}, pairs{2}};
            success = true;
            failure_reason = "";
        otherwise
            error('evaluate_auto_group_ssa:UnknownSignal', 'Unknown signal type %s', signal_type);
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
