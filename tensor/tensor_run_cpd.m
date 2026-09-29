function result = tensor_run_cpd(tensor_input, params)
%TENSOR_RUN_CPD Run Tensorlab CPD and package diagnostics.
%
%   result = tensor_run_cpd(T, params)
%   result = tensor_run_cpd(tensor_build_decimated_hankel(...), params)

    params = tensor_merge_params(params);
    [T, build_metadata] = parse_input(tensor_input);

    result = struct();
    result.success = false;
    result.failure_reason = "";
    result.method = string(params.cpd_method);
    result.rank = params.cpd_rank;
    result.tensor_size = size(T);
    result.tensor_norm = norm(T(:));
    result.lambda = [];
    result.factors = {};
    result.fit = nan;
    result.rel_error = nan;
    result.diagnostics = struct();
    result.raw_output = struct();
    result.build_metadata = build_metadata;

    if ndims(T) < 3
        result.failure_reason = "Tensor input must be at least third-order.";
        return;
    end

    method_name = char(params.cpd_method);
    if exist(method_name, 'file') ~= 2 && exist(method_name, 'builtin') ~= 5
        result.failure_reason = "Requested Tensorlab CPD method is not on the MATLAB path.";
        return;
    end

    try
        switch lower(method_name)
            case 'cpd'
                if isempty(params.cpd_initialization)
                    [U, output] = cpd(T, params.cpd_rank, params.cpd_options);
                else
                    [U, output] = cpd(T, params.cpd_initialization, params.cpd_options);
                end
            case 'cpd_nls'
                if isempty(params.cpd_initialization)
                    U0 = initial_guess(size(T), params.cpd_rank);
                else
                    U0 = params.cpd_initialization;
                end
                [U, output] = cpd_nls(T, U0, params.cpd_options);
            otherwise
                error('tensor_run_cpd:UnknownMethod', ...
                    'Unsupported CPD method "%s".', method_name);
        end
    catch ME
        result.failure_reason = string(sprintf('%s: %s', ME.identifier, ME.message));
        return;
    end

    [U, lambda] = normalize_factors(U);
    [fit, rel_error] = approximation_quality(T, U, lambda);

    result.success = true;
    result.lambda = lambda;
    result.factors = U;
    result.fit = fit;
    result.rel_error = rel_error;
    result.raw_output = output;
    result.diagnostics = build_diagnostics(U, lambda, output, fit, rel_error);
end

function [T, build_metadata] = parse_input(tensor_input)
    build_metadata = struct();
    if isnumeric(tensor_input)
        T = tensor_input;
        return;
    end
    if isstruct(tensor_input) && isfield(tensor_input, 'tensor')
        T = tensor_input.tensor;
        build_metadata = tensor_input;
        return;
    end
    error('tensor_run_cpd:BadInput', ...
        'Input must be a numeric tensor or tensor_build_decimated_hankel output.');
end

function U0 = initial_guess(tensor_size, rank_k)
    if exist('cpd_rnd', 'file') == 2 || exist('cpd_rnd', 'builtin') == 5
        U0 = cpd_rnd(tensor_size, rank_k);
    else
        U0 = cell(numel(tensor_size), 1);
        for m = 1:numel(tensor_size)
            U0{m} = randn(tensor_size(m), rank_k);
        end
    end
end

function [U_norm, lambda] = normalize_factors(U)
    num_modes = numel(U);
    rank_k = size(U{1}, 2);
    U_norm = U;
    lambda = ones(rank_k, 1);

    for r = 1:rank_k
        for m = 1:num_modes
            col = U_norm{m}(:, r);
            scale = norm(col);
            if scale > 0
                U_norm{m}(:, r) = col / scale;
                lambda(r) = lambda(r) * scale;
            end
        end

        [~, idx] = max(abs(U_norm{1}(:, r)));
        if ~isempty(idx) && U_norm{1}(idx, r) < 0
            U_norm{1}(:, r) = -U_norm{1}(:, r);
            lambda(r) = -lambda(r);
        end
    end

    [~, order] = sort(abs(lambda), 'descend');
    lambda = lambda(order);
    for m = 1:num_modes
        U_norm{m} = U_norm{m}(:, order);
    end
end

function [fit, rel_error] = approximation_quality(T, U, lambda)
    fit = nan;
    rel_error = nan;

    if ~(exist('cpdgen', 'file') == 2 || exist('cpdgen', 'builtin') == 5)
        return;
    end

    try
        U_weighted = U;
        U_weighted{1} = U_weighted{1} * diag(lambda);
        T_hat = cpdgen(U_weighted);
        rel_error = norm(T(:) - T_hat(:)) / max(norm(T(:)), eps);
        fit = 1 - rel_error;
    catch
        fit = nan;
        rel_error = nan;
    end
end

function diagnostics = build_diagnostics(U, lambda, output, fit, rel_error)
    diagnostics = struct();
    diagnostics.fit = fit;
    diagnostics.rel_error = rel_error;
    diagnostics.output_fields = string(fieldnames(output));
    diagnostics.component_weights = lambda;
    diagnostics.component_weight_abs = abs(lambda);
    diagnostics.mode_column_norms = cellfun(@(A) vecnorm(A, 2, 1).', U, 'UniformOutput', false);
    diagnostics.num_modes = numel(U);
    diagnostics.rank = numel(lambda);

    [iterations, has_iterations] = try_scalar_field(output, ...
        {'iterations', 'Iterations', 'iter', 'Iter'});
    if has_iterations
        diagnostics.iterations = iterations;
    else
        diagnostics.iterations = nan;
    end

    [obj, has_obj] = try_vector_field(output, ...
        {'fval', 'Fval', 'relerr', 'RelErr', 'objective', 'Objective'});
    if has_obj
        diagnostics.objective_trace = obj;
    else
        diagnostics.objective_trace = [];
    end
end

function [value, ok] = try_scalar_field(s, names)
    ok = false;
    value = nan;
    for i = 1:numel(names)
        name = names{i};
        if isfield(s, name)
            candidate = s.(name);
            if isnumeric(candidate) && isscalar(candidate)
                value = candidate;
                ok = true;
                return;
            end
        end
    end
end

function [value, ok] = try_vector_field(s, names)
    ok = false;
    value = [];
    for i = 1:numel(names)
        name = names{i};
        if isfield(s, name)
            candidate = s.(name);
            if isnumeric(candidate)
                value = candidate(:);
                ok = true;
                return;
            end
        end
    end
end
