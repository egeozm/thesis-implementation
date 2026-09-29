function [grouped_est, grouping, component_names] = group_cp_terms_benchmark(est_rank1, lambda, signal_type, params)
%GROUP_CP_TERMS_BENCHMARK Benchmark-oriented grouping for thesis synthetics.
%
%   [grouped_est, grouping, component_names] = group_cp_terms_benchmark(est_rank1, lambda, signal_type, params)
%
%   This helper is intentionally tied to the synthetic benchmark families
%   used throughout the repository. It maps reconstructed CP rank-1 terms onto
%   the known component families so that tensor outputs can be compared
%   directly against SSA/SSD benchmark tables on the same synthetic tasks.

    if nargin < 4
        params = struct();
    end
    params = tensor_merge_params(params);

    est_rank1 = est_rank1(:,:);
    lambda = lambda(:);
    [N, rank_k] = size(est_rank1);
    if numel(lambda) ~= rank_k
        error('group_cp_terms_benchmark:SizeMismatch', ...
            'lambda must have one entry per reconstructed rank-1 term.');
    end

    specs = benchmark_specs(signal_type, N, params);
    num_groups = numel(specs);
    component_names = {specs.name};
    grouped_est = zeros(N, num_groups);
    assignments = zeros(rank_k, 1);
    assignment_scores = zeros(rank_k, num_groups);
    dominant_freqs = nan(rank_k, 1);
    sign_flips = ones(rank_k, 1);
    members = cell(num_groups, 1);

    for r = 1:rank_k
        term = est_rank1(:, r);
        dominant_freqs(r) = estimate_dominant_frequency(term, params.fs);
        for g = 1:num_groups
            assignment_scores(r, g) = score_term(term, dominant_freqs(r), specs(g), params);
        end
        [~, assignments(r)] = max(assignment_scores(r, :));
    end

    x_hat = sum(est_rank1, 2);
    [grouped_est, coeff_blocks] = project_total_reconstruction(x_hat, specs);

    for g = 1:num_groups
        members{g} = find(assignments == g).';
    end

    grouping = struct();
    grouping.mode = "benchmark";
    grouping.signal_type = string(signal_type);
    grouping.params = params;
    grouping.component_names = component_names;
    grouping.specs = specs;
    grouping.lambda = lambda;
    grouping.dominant_freqs = dominant_freqs;
    grouping.assignment_scores = assignment_scores;
    grouping.assignments = assignments;
    grouping.sign_flips = sign_flips;
    grouping.component_members = members;
    grouping.total_reconstruction = x_hat;
    grouping.component_coefficients = coeff_blocks;
    grouping.num_grouped_components = num_groups;
end

function specs = benchmark_specs(signal_type, N, params)
    n = (1:N).';
    n_norm = n / N;

    switch lower(string(signal_type))
        case "multiscale"
            specs = [ ...
                make_spec("trend", nan, [ones(N, 1), n_norm, n_norm.^2], true), ...
                make_spec("slow", 0.01, sinusoid_basis(n, 0.01, params.fs), false), ...
                make_spec("fast", 0.08, sinusoid_basis(n, 0.08, params.fs), false) ...
                ];
        case "close"
            specs = [ ...
                make_spec("comp1", 0.05, sinusoid_basis(n, 0.05, params.fs), false), ...
                make_spec("comp2", 0.055, sinusoid_basis(n, 0.055, params.fs), false) ...
                ];
        otherwise
            error('group_cp_terms_benchmark:UnknownSignal', ...
                'Unknown benchmark signal type %s.', string(signal_type));
    end
end

function spec = make_spec(name, target_freq, basis, is_trend)
    spec = struct();
    spec.name = char(name);
    spec.target_freq = target_freq;
    spec.basis = orthonormal_basis(basis);
    spec.is_trend = is_trend;
end

function basis = sinusoid_basis(n, freq, fs)
    omega = 2 * pi * freq / max(fs, eps);
    basis = [sin(omega * n), cos(omega * n)];
end

function score = score_term(term, dominant_freq, spec, params)
    projected = project_onto_basis(term, spec.basis);
    proj_ratio = sum(projected.^2) / max(sum(term.^2), eps);

    if spec.is_trend
        if dominant_freq < params.benchmark_trend_freq_max
            freq_weight = 1;
        else
            freq_weight = 0.1;
        end
    else
        if isnan(dominant_freq) || dominant_freq <= 0
            freq_weight = 0.1;
        else
            freq_diff = abs(dominant_freq - spec.target_freq);
            freq_weight = 1 / (1 + freq_diff / max(params.benchmark_freq_tol, eps));
        end
    end

    score = proj_ratio * freq_weight;
end

function Q = orthonormal_basis(B)
    if isempty(B)
        Q = zeros(size(B, 1), 0);
        return;
    end
    [Q, ~] = qr(B, 0);
end

function projected = project_onto_basis(x, basis)
    if isempty(basis)
        projected = zeros(size(x));
        return;
    end
    projected = basis * (basis' * x(:));
end

function [grouped_est, coeff_blocks] = project_total_reconstruction(x_hat, specs)
    num_groups = numel(specs);
    basis_sizes = zeros(num_groups, 1);
    full_basis = [];
    for g = 1:num_groups
        Bg = specs(g).basis;
        basis_sizes(g) = size(Bg, 2);
        full_basis = [full_basis, Bg]; %#ok<AGROW>
    end

    if isempty(full_basis)
        grouped_est = zeros(numel(x_hat), num_groups);
        coeff_blocks = cell(num_groups, 1);
        return;
    end

    coeffs = full_basis \ x_hat(:);
    grouped_est = zeros(numel(x_hat), num_groups);
    coeff_blocks = cell(num_groups, 1);
    cursor = 0;
    for g = 1:num_groups
        block = coeffs((cursor + 1):(cursor + basis_sizes(g)));
        coeff_blocks{g} = block;
        grouped_est(:, g) = specs(g).basis * block;
        cursor = cursor + basis_sizes(g);
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
