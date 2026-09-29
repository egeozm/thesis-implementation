function reconstructed = ssa_reconstruct_groups(features, groups)
%SSA_RECONSTRUCT_GROUPS Reconstruct grouped SSA components from precomputed features.
%
%   reconstructed = ssa_reconstruct_groups(features, groups)

    num_groups = numel(groups);
    reconstructed = zeros(features.N, num_groups);

    for g = 1:num_groups
        idx = groups{g}(:).';
        idx = idx(idx >= 1 & idx <= features.rank);
        if isempty(idx)
            continue;
        end

        HG = zeros(features.L, features.K);
        for r = idx
            HG = HG + features.S(r, r) * (features.U(:, r) * features.V(:, r).');
        end
        reconstructed(:, g) = diagonal_average(HG, features.N);
    end
end
