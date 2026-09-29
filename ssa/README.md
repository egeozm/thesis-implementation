# SSA baseline

This folder contains the repository's classical Singular Spectrum Analysis (SSA) baseline. In this project, SSA is the fixed single-matrix comparator against which the later adaptive SSD baseline and the tensor-based methods are interpreted.

## Signal embedding

Given a univariate signal `x` of length `N` and a window length `L`, the implementation builds the standard trajectory/Hankel matrix

- `H in R^(L x K)`
- `K = N - L + 1`
- `H(i,j) = x(i + j - 1)`

using `embed_hankel.m`.

This is the only matrix construction used by the SSA baseline: one signal, one window length, one Hankel matrix.

## Chosen window lengths

The main SSA experiments in this repository use

- `L = [200, 400, 600, 1000]` for calibration, frozen experiments, and sensitivity studies
- `L_ref = 1000` as the main reference setting in `scripts/run_ssa_verification.m`

All experiments use signals of length `N = 2000`, so the tested windows cover a range from relatively local to half-length embeddings.

The purpose of keeping several fixed `L` values is to show explicitly how much the quality of classical SSA depends on the chosen embedding scale.

## Decomposition step

For a chosen `L`, the baseline performs the standard SSA factorization

`H = U S V^T`

using MATLAB's economy SVD:

- `svd(H, 'econ')`
- singular values are stored in descending order
- each eigentriple `(u_r, sigma_r, v_r)` defines one elementary rank-1 matrix

`ssa_precompute.m` reconstructs the leading elementary components and estimates a dominant frequency for each one. Those features are then used by the fixed automatic grouping rule.

## Grouping and reconstruction rule

The main automatic grouping logic is implemented in `auto_group_ssa.m`.

For the first `r_max` eigentriples, the code inspects adjacent pairs `(r, r+1)` and accepts them as oscillatory candidates only if all of the following hold:

- the pair has non-negligible spectral energy
- the singular-value balance satisfies `sv_ratio >= sv_ratio_thresh`
- the dominant-frequency gap satisfies `|f_r - f_(r+1)| <= freq_tol`
- the mean dominant frequency is at least `min_osc_freq`

Accepted adjacent pairs are scored by

`score(r, r+1) = (sigma_r + sigma_(r+1)) * min(sigma_r, sigma_(r+1)) / max(sigma_r, sigma_(r+1)) / (1 + |f_r - f_(r+1)| / max(freq_tol, eps))`

The implementation then:

1. keeps the best two non-overlapping oscillatory pairs
2. forms the trend group from remaining low-frequency unpaired leading components
3. treats everything else as residual

Grouped reconstruction is done in `ssa_reconstruct_groups.m`:

1. sum the selected elementary matrices for a group
2. apply diagonal averaging (Hankelization)
3. return one reconstructed time-domain series per group

For the thesis synthetic tasks, the intended grouped outputs are:

- `multiscale`: `{trend, slow pair, fast pair}`
- `close`: `{pair 1, pair 2}`

## Fixed implementation choices

Several choices are intentionally frozen in this repository's SSA baseline:

- the method is always single-resolution SSA on one Hankel matrix
- only adjacent eigentriple pairs are considered as oscillatory pairs
- at most two oscillatory pairs are retained by the automatic grouping rule
- the trend group is built from low-frequency unpaired leading components rather than from a separate optimization step
- evaluation uses the same synthetic tasks for all tested `L` values: `multiscale` and `close`

The grouping hyperparameters are calibrated once in `scripts/run_ssa_calibration.m` over:

- `r_max = [8, 10, 12]`
- `sv_ratio_thresh = [0.5, 0.6, 0.7, 0.8]`
- `freq_tol = [0.005, 0.01, 0.02]`
- `min_osc_freq = [0.005, 0.01, 0.015]`

Then `select_ssa_params.m` freezes one parameter set by choosing from configurations within 5% of the best score and tie-breaking toward simpler settings. The frozen parameters are reused in later SSA batch experiments and also passed into SSD so that the matrix baselines stay comparable.

## Files

- `embed_hankel.m`: build the trajectory matrix
- `ssa_decompose.m`: direct SSA pipeline for user-supplied groups
- `ssa_precompute.m`: SVD plus leading elementary reconstructions and frequency features
- `auto_group_ssa.m`: fixed automatic grouping rule
- `ssa_reconstruct_groups.m`: grouped reconstruction by diagonal averaging
- `evaluate_auto_group_ssa.m`: evaluate automatic grouping against synthetic truth
- `select_ssa_params.m`: choose the frozen grouping parameter set
