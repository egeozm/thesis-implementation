function y = ssd_extract_mode_at_features(features, pair)
%SSD_EXTRACT_MODE_AT_FEATURES Reconstruct one SSD mode from an eigentriple pair.
%
%   y = ssd_extract_mode_at_features(features, pair)
%
%   pair: 1x2 indices into leading singular triples (same convention as auto_group_ssa).

    if isempty(pair)
        y = zeros(features.N, 1);
        return;
    end
    y = ssa_reconstruct_groups(features, {pair});
    y = y(:);
end
