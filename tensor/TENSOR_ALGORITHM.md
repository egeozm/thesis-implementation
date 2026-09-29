# Decimated Hankel Tensor + CPD

This module implements the finalized tensor CPD workflow for the thesis:

1. build a third-order tensor from full-coverage decimated Hankel blocks of a 1D signal,
2. fit a CP decomposition with Tensorlab,
3. reconstruct CP rank-1 terms back to full-length 1D estimates,
4. group CP rank-1 terms into higher-level 1D component estimates,
5. return factors, diagnostics, and component-evaluation hooks.

## Tensorization Contract

### Inputs

- `x`: column vector of length `N`.
- `params`: merged via `tensor_merge_params`, with key fields:
  - `strides`: positive integer decimation factors, e.g. `[1 2 4 8]`
  - `base_window_length`: target window span on the original sample grid
  - `start_index`: first original-signal sample used in each decimation
  - `block_num_rows`, `block_num_cols`: shared block size used inside each Hankel matrix
  - `block_row_step`, `block_col_step`: block extraction strides inside each Hankel matrix
  - `min_common_rows`, `min_common_cols`: minimum admissible shared block size
  - `min_num_slices`: minimum number of valid decimated slices
  - `min_blocks_per_stride`: minimum number of extracted blocks required per stride
  - `recon_weighting`: overlap-add weighting policy for inverse mapping
  - `reconstruction_mode`: one of `cpd`, `svd`, `blockwise`
  - `cpd_rank`, `cpd_method`, `cpd_options`: CPD settings
  - `grouping_mode`: one of `benchmark`, `auto`, `oracle`, `none`
  - `group_freq_tol`, `group_rho_merge_thresh`, `group_min_osc_freq`: automatic grouping controls
  - `group_oracle_min_rho`: truth-guided grouping threshold
  - `benchmark_freq_tol`, `benchmark_trend_freq_max`: synthetic-benchmark grouping controls

### Slice Construction

For each stride `s`:

1. Form decimated indices `idx_s = start_index:s:N`.
2. Extract the decimated signal `x_s = x(idx_s)`.
3. Convert the base window length from original-grid samples to decimated samples:

   `L_s = floor((base_window_length - 1) / s) + 1`

   so that the retained rows correspond to roughly the same original-time span.
4. Build the Hankel slice `H_s = embed_hankel(x_s, L_s)`.

## Common-Shape Policy

All stride-specific Hankel matrices are converted into shared-size local blocks:

- compute each valid full-slice size `(L_s, K_s)`,
- define
  - `L_block <= min_s L_s`
  - `K_block <= min_s K_s`
- extract many valid sub-blocks `H_s(r:r+L_block-1, c:c+K_block-1)`,
- stack those blocks along the third mode.

This module does **not** pad missing values. Every tensor block contains only observed Hankel entries, while coverage of the full signal is achieved by collecting many overlapping blocks rather than one shared top-left crop.

## Metadata Contract

The tensor builder returns:

- `tensor`: numeric array of size `L_block x K_block x B`
- `tensor_size`: same dimensions as a row vector
- `valid_strides`: strides that produced admissible slices
- `stride_info`: per-stride metadata including
  - raw Hankel sizes,
  - row/column block starts,
  - number of blocks per stride
- `block_info`: per-block metadata including
  - `stride`
  - `row_start`, `row_end`
  - `col_start`, `col_end`
  - retained original-grid indices covered by the block
- `skipped_slices`: stride-specific validation failures
- `common_shape`: struct with `rows`, `cols`, `num_blocks`, `num_strides`
- `coverage`: counts and uncovered indices on the original time grid

This metadata is intended to support later inverse mapping from CP factors back to the original 1D signal.

## CPD Contract

`tensor_run_cpd` accepts either:

- the numeric tensor itself, or
- the full tensor-builder output struct.

The CPD wrapper:

- checks whether Tensorlab is available,
- runs the requested CPD method (`cpd` by default),
- normalizes factor columns into a separate weight vector `lambda`,
- reports fit diagnostics and raw Tensorlab output.

## Reconstruction Modes

`evaluate_tensor_cpd` supports three inverse-map variants:

- `cpd` (default): reconstruct each CP rank-1 term into a full stride-level
  Hankel estimate by reassembling overlapping CP-generated blocks in matrix
  space, diagonal-average that stride-level term Hankel to the decimated grid,
  then fuse the stride-specific term estimates back onto the original signal
  grid. This is the repository's fully-CPD path because it does not apply any
  post-CPD matrix factorization.
- `svd`: first rebuild a stride-level Hankel approximation from the CP factors,
  then run an SSA-like SVD/grouping/diagonal-averaging stage on each stride.
  This is a hybrid CPD+SSA diagnostic path.
- `blockwise`: the original debugging baseline that diagonal-averages each
  reconstructed block directly onto the original grid, mixing all strides in one
  overlap-add pass.

Under `cpd`, each rank-1 waveform is formed by:

1. reconstructing every retained block contribution for that CP term,
2. reassembling those contributions into a full Hankel matrix within each
   stride by averaging overlapping block entries,
3. diagonal-averaging the full stride-level term Hankel to the decimated grid,
4. combining the stride-specific term estimates back on the original grid.

## CP Term Grouping

`evaluate_tensor_cpd` supports four grouping modes:

- `none`: use the raw rank-1 waveforms directly as candidate components. This
  is kept as a debugging baseline.
- `benchmark`: use the known synthetic benchmark families to map the **full CPD
  reconstruction** into the same final component layout used by SSA/SSD
  evaluation. For the thesis synthetic tasks this is the intended reporting
  mode.
- `oracle`: assign each rank-1 waveform to the truth component with the
  strongest absolute correlation, optionally flipping its sign before summing.
  This is only for synthetic-signal diagnostics and provides an upper bound on
  what the current reconstruction can recover.
- `auto`: cluster rank-1 waveforms by dominant-frequency proximity and waveform
  correlation, align their signs within each cluster, then sum each cluster
  into a higher-level component estimate.

For the RQ3 financial experiments, no ground-truth components are available, so
`benchmark` and `oracle` grouping are not used. Instead,
`scripts/run_financial_decomposition.m` reports both `auto` grouping and `none`
grouping for the same fixed financial grid. The financial driver also bypasses
`evaluate_tensor_cpd`, which is synthetic-evaluation oriented, and directly
chains `tensor_build_decimated_hankel` -> `tensor_run_cpd` ->
`tensor_reconstruct_via_cpd` before applying the selected financial grouping
mode.

For the reported benchmark mode:

1. CPD reconstructs the full signal estimate from the tensor factors.
2. The reconstructed signal is projected onto benchmark subspaces:
   - `multiscale`: quadratic trend basis plus slow/fast sinusoidal bases,
   - `close`: the two known sinusoidal bases.
3. The evaluation wrapper appends a residual column
   `x - sum(grouped_components, 2)` and scores the resulting full-length
   columns against synthetic truth using `match_components` and
   `compute_metrics`.

## Deferred Work

This tensor pipeline still does not yet:

- provide a signal-agnostic automatic grouping rule that matches the
  benchmark-oriented grouping quality on these synthetic tasks.
