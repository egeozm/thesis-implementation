# Singular Spectrum Decomposition (SSD) — implementation spec

This implementation follows the thesis description of SSD (Bonizzi et al., *Singular Spectrum Decomposition*, Signal Processing, 2014): an **iterative**, **single-resolution** procedure that repeatedly applies an SSA-like embedding and SVD on the **current residual**, extracts a **dominant oscillatory mode** (as one reconstructed component), subtracts it, and continues until a stopping rule fires.

## Inputs

- `x`: column vector, length `N`.
- `params` (struct), merged via `ssd_merge_params`:
  - **From frozen SSA calibration** (same semantics as `auto_group_ssa`): `r_max`, `sv_ratio_thresh`, `freq_tol`, `min_osc_freq`.
  - **SSD-specific**: `L_list` (candidate window lengths), `max_modes` (max extractions), `min_mode_energy_frac` (ignore extractions smaller than this fraction of the residual norm before subtraction), `min_modes_before_score_stop`, `min_score_ratio`, `min_score_drop_ratio`, `fs` (sampling rate for FFT features, default `1`), `extract_trend_first` (default `false`), and optional `trend_L`.

## Per-iteration steps

If `extract_trend_first = true`, the implementation first computes an SSA
decomposition of the original signal using `trend_L` when provided, otherwise
the candidate `L_list` in descending order. It reconstructs the leading
low-frequency eigentriple or eigentriples selected by the same frozen SSA trend
rule, subtracts that estimate from the signal, and then starts the oscillatory
SSD loop below on the trend-removed residual. With the default
`extract_trend_first = false`, SSD remains oscillatory-only.

1. **Window selection** (`ssd_select_window`): For each `L` in `L_list` with `2 <= L < N`, form the Hankel trajectory matrix of the residual, run SVD, and compute leading elementary reconstructions and dominant frequencies (via `ssa_precompute` + the same adjacent-pair feasibility rules as `auto_group_ssa`). Among all feasible **single** adjacent oscillatory pairs across all `L`, choose the pair with **largest** internal `pair_score` (same construction as in `auto_group_ssa`: balances singular-value mass, singular-value balance, and frequency tightness). Tie-break: larger `L`, then lexicographic pair index.

## Pair score

For an adjacent candidate pair `(r, r+1)` with singular values `sigma_r`, `sigma_{r+1}` and dominant frequencies `f_r`, `f_{r+1}`, define

- `sigma_mass = sigma_r + sigma_{r+1}`
- `sv_ratio = min(sigma_r, sigma_{r+1}) / max(sigma_r, sigma_{r+1})`
- `freq_diff = |f_r - f_{r+1}|`
- `freq_penalty = 1 + freq_diff / max(freq_tol, eps)`

The score used in both the implementation and the thesis experiments is the multiplicative heuristic

`score(r, r+1) = sigma_mass * sv_ratio / freq_penalty`

There are no additional learned or tuned weights in this score: the only tuning entering the frequency term is `freq_tol`, which also appears in the feasibility test.

## Clarifying what “best” means across `L`

The implementation considers **all feasible adjacent pairs across all candidate `L` values**, but computes that efficiently:

1. For each `L`, rank all feasible adjacent pairs by the score above.
2. Keep only the top-scoring pair for that `L`.
3. Choose the single best pair across `L` by comparing those per-`L` maxima.

This is equivalent to taking the global maximum over all feasible adjacent pairs across all `L`, because a non-maximal pair within one fixed `L` can never beat that same `L`'s maximal pair.

2. **Mode extraction** (`ssd_extract_mode_at_features`): Reconstruct the chosen pair as **one** time-domain mode using grouped Hankelization (`ssa_reconstruct_groups`).

3. **Residual update**: `residual := residual - mode`.

4. **Stopping** (`ssd_stop_rule`):
   - No feasible pair at any `L` in `L_list`.
   - `norm(mode) / max(norm(residual_before), eps) < min_mode_energy_frac`.
   - After the first `min_modes_before_score_stop` accepted modes, stop if the candidate pair score is too small relative to the **first** accepted score (`min_score_ratio`) or falls too sharply relative to the **previous** accepted score (`min_score_drop_ratio`).
   - Iteration count reaches `max_modes`.

## Outputs (`ssd_decompose`)

- `modes`: `N × M` matrix, columns = extracted oscillatory modes in order.
- `trend`: `N × 1` trend estimate when `extract_trend_first` succeeds, otherwise zeros.
- `trend_extracted`, `trend_L`, `trend_idx`: diagnostics for the optional trend extraction.
- `residual`: final residual (trend + noise + non-oscillatory remainder).
- `L_history`, `pair_history`, `scores`: diagnostics per iteration.

## Evaluation vs ground truth

- Build `est_components = [modes, residual]` so there are `M+1` columns for `match_components` / `compute_metrics`. When `extract_trend_first` is enabled and succeeds, `evaluate_ssd` instead uses `[trend, modes, residual]`.
- In the `multiscale` synthetic case, the final residual is interpreted as a **trend-dominated remainder**, not as a guaranteed pure trend estimate. It may still contain residual noise or low-frequency content not absorbed by the extracted oscillatory modes.
- **Success flags** are only a minimal implementation check for whether SSD extracted enough oscillatory modes to attempt the comparison: `multiscale` requires `M >= 2` (slow + fast candidates before the final residual), and `close` requires `M >= 2`.
- The primary evaluation remains the actual component metrics (`NMSE`, Pearson correlation, peak-frequency error, leakage). Extracting two modes is necessary for the current implementation to proceed, but it is **not** sufficient to claim correct separation.

## Relation to classical SSA in this repo

- Reuses `ssa/embed_hankel`, `ssa/ssa_precompute`, `ssa/ssa_reconstruct_groups`, `ssa/diagonal_average` (indirectly), `ssa/match_components`, `ssa/compute_metrics`.
- **Does not** reuse frozen automatic **grouping** as a one-shot triple (trend + two pairs); instead it repeatedly extracts **one** oscillatory pair per iteration with **adaptive** `L`.

## Limitations

- This is a **practical MATLAB realization** aligned with the thesis narrative; full parity with every detail of Bonizzi et al. (e.g. alternative subspace or window criteria) is not guaranteed without their exact procedural specification.
