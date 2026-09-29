# SSD baseline

This folder contains the repository's Singular Spectrum Decomposition (SSD) baseline. In this thesis codebase, SSD is used as a controlled adaptive matrix comparator: it still works with one Hankel matrix at a time, but it updates the residual and can change the window length between extraction steps.

## How SSD is implemented here

The implementation follows an iterative SSA-like pattern:

1. start from the current residual, initially the observed signal
2. try several candidate window lengths
3. build a Hankel matrix for each candidate window
4. run SVD and identify the best feasible adjacent oscillatory pair
5. reconstruct that pair as one mode
6. subtract the mode from the residual
7. repeat until a stopping rule fires

The main loop is in `ssd_decompose.m`.

Optionally, `ssd_decompose.m` can extract a low-frequency trend before the
oscillatory loop by setting `params.extract_trend_first = true`. In that mode,
the implementation reconstructs the leading low-frequency eigentriple or
eigentriples from the original signal, subtracts this trend estimate, and then
continues with the same adaptive oscillatory-pair extraction on the residual.
The default remains oscillatory-only extraction for the synthetic SSD
experiments.

## Adaptive window selection

At each iteration, `ssd_select_window.m` evaluates every `L` in `L_list`. The default candidate set in this repository is

- `L_list = [200, 400, 600, 1000]`

For each candidate `L`, the code:

1. embeds the current residual into a trajectory/Hankel matrix
2. runs `ssa_precompute` to obtain singular values, elementary reconstructions, and dominant frequencies
3. ranks feasible adjacent oscillatory pairs using the same feasibility tests used by the frozen SSA grouping rule
4. keeps only the top-scoring pair for that `L`

Across window lengths, the implementation then chooses the single best pair overall. Tie-breaks prefer:

1. larger `L`
2. lexicographically smaller pair indices

So the SSD baseline is adaptive in `L`, but every extraction step is still based on one ordinary Hankel matrix.

## Pair selection criterion

The same adjacent-pair feasibility thresholds as SSA are reused:

- `r_max`
- `sv_ratio_thresh`
- `freq_tol`
- `min_osc_freq`

For a feasible adjacent pair `(r, r+1)`, the score is

`score(r, r+1) = (sigma_r + sigma_(r+1)) * min(sigma_r, sigma_(r+1)) / max(sigma_r, sigma_(r+1)) / (1 + |f_r - f_(r+1)| / max(freq_tol, eps))`

This favors:

- large singular-value mass
- balanced singular values inside the pair
- small dominant-frequency mismatch

There are no extra learned weights beyond this formula.

## Residual update

Once the best pair and window length are chosen, `ssd_extract_mode_at_features.m` reconstructs the pair as one time-domain mode using SSA grouped reconstruction. The residual is then updated as

`residual := residual - mode`

The next iteration is performed on that updated residual, not on the original signal.

This is the main behavioral difference from classical SSA in this repository: SSA performs one decomposition at a fixed `L`, while SSD repeatedly decomposes the evolving residual.

## Stopping rule

The stopping logic is implemented in `ssd_stop_rule.m`. Extraction stops when any of the following happens:

- no feasible adjacent pair exists for any `L` in `L_list`
- the candidate mode energy is too small relative to the current residual norm
- after enough accepted modes, the pair score becomes too weak relative to earlier accepted scores
- the iteration limit is reached

The default SSD-specific thresholds are:

- `max_modes = 4`
- `min_mode_energy_frac = 1e-2`
- `min_modes_before_score_stop = 2`
- `min_score_ratio = 0.1`
- `min_score_drop_ratio = 0.5`

## Fixed implementation choices

This repository deliberately fixes several SSD details to keep the baseline controlled and comparable:

- SSD reuses the frozen SSA calibration parameters instead of introducing a separate SSD calibration stage
- only adjacent eigentriple pairs are eligible for extraction
- only one pair is extracted per iteration
- reconstruction is done through the same SSA regrouping machinery already used elsewhere in the repo
- the final residual is kept as a trend-dominated remainder, especially for the `multiscale` signal, rather than forcing a separate pure-trend extraction stage

For batch experiments in `scripts/run_ssd_experiments.m`, trial rows store `L = NaN` because the effective window length is adaptive over iterations. The number of extracted oscillatory modes is stored separately.

## Financial usage

For the RQ3 financial experiments, `scripts/run_financial_decomposition.m`
applies the same SSD implementation to each standardized level series with
`L_list = [200, 400, 600, 1000]` and `extract_trend_first = true`. The financial
path reuses the frozen SSA feasibility thresholds and does not introduce a
separate asset-specific calibration stage. Because EUR/USD and S&P 500 do not
provide ground-truth latent components, the trend estimate, extracted
oscillatory modes, and final residual are interpreted qualitatively.

## Relation to the SSA baseline

Relative to the SSA baseline in `ssa/`, this SSD implementation changes four things:

- window length is selected adaptively from `L_list`
- decomposition is repeated on the residual
- only one oscillatory pair is taken per iteration
- stopping is controlled by residual-energy and score-based criteria

What stays fixed is just as important:

- each extraction still uses an ordinary single-resolution Hankel matrix
- oscillatory-pair feasibility is still defined by the same frozen thresholds
- the baseline remains matrix-based, not tensor-based

## Files

- `SSD_ALGORITHM.md`: detailed implementation spec
- `ssd_merge_params.m`: merge frozen SSA parameters with SSD defaults
- `ssd_decompose.m`: main iterative SSD loop
- `ssd_select_window.m`: adaptive window and pair selection
- `ssd_list_oscillatory_pairs.m`: rank feasible adjacent pairs
- `ssd_stop_rule.m`: stopping conditions
- `evaluate_ssd.m`: compare extracted modes and final residual against synthetic truth
