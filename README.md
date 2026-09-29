# Thesis implementation: SSA and Hankel tensor methods

MATLAB code for studying how different Hankelization strategies affect component separation in univariate time series. The current repository contains:

- **Classical Singular Spectrum Analysis (SSA)** as the main single-matrix baseline
- **Singular Spectrum Decomposition (SSD)** as a controlled adaptive matrix-based comparator
- a finalized **tensor / CPD** benchmark pipeline based on decimated Hankel slices

The main thesis focus is the comparison between **single-resolution matrix Hankelization** and later **multi-resolution / tensor Hankelization**. In that context, SSD is included to show what an adaptive matrix-based method can achieve before moving to the tensor construction; it is **not** the central methodological contribution of the repository.

## Requirements

- MATLAB R2025b (or recent release with `svd(..., 'econ')`, `hankel`, standard toolboxes)

## Setup

1. Add the project root to the MATLAB path (or run scripts from the project root):

   ```matlab
   cd('/path/to/thesis-implementation');
   addpath(fullfile(pwd, 'ssa'));
   addpath(fullfile(pwd, 'ssd'));
   addpath(fullfile(pwd, 'tensor'));
   ```

2. **Tensorlab** (required for the CPD / tensor experiments): this archive includes
   Tensorlab files under `tensorlab/`. Add them with:

   ```matlab
   addpath(genpath(fullfile(pwd, 'tensorlab')));
   ```

   See [tensorlab/README.md](tensorlab/README.md) and
   [tensorlab/license_tensorlab.txt](tensorlab/license_tensorlab.txt) for layout
   and license details.

## Generated artifacts and run order

The repository includes the main CSV result tables used for reporting. The
MATLAB `.mat` files and figure exports are generated artifacts: they are written
when the corresponding scripts are run, but they are not required to inspect the
included CSV tables directly.

Some downstream aggregation and plotting scripts load `.mat` workspace files
instead of reconstructing all state from CSVs. To reproduce the full pipeline
from a fresh checkout, run the generating scripts in order, starting with
`scripts/run_ssa_calibration.m`, then the frozen SSA / SSD / tensor benchmark
scripts, and finally the comparison and visualization scripts.

## Synthetic signals and reproducibility

All SSA and SSD experiments in this repository use the generator in `ssa/generate_signal.m`. There are two synthetic signal families:

- `multiscale`: a smooth quadratic trend plus a slow sinusoid (`f = 0.01`), a faster sinusoid (`f = 0.08`) with random phase, and additive Gaussian noise
- `close`: two nearby sinusoids (`f1 = 0.05`, `f2 = 0.055`), one with random phase, plus additive Gaussian noise

For each call,

```matlab
[x, components, t] = generate_signal(signal_type, N, snr_db, seed);
```

the generator:

1. sets the random number generator with `rng(seed)` when a seed is provided
2. builds the clean synthetic signal from the deterministic formulas above
3. computes the noise level from the requested `snr_db`
4. draws the random phase and Gaussian noise
5. returns the observed signal `x`, the ground-truth components in `components`, and the sample index vector `t`

This means the synthetic signals are **not** stored explicitly for every trial. Instead, they are regenerated on demand from the metadata recorded in the experiment tables, especially:

- `signal_type`
- `N`
- `snr_db`
- `seed`

For SSA experiments, `L` is also recorded because it affects the decomposition, but not the signal generation itself.

Because the same seed is reused, the same random phase and the same noise realization are reproduced exactly. In practice, this lets you recreate any experiment signal later without saving every waveform, as long as `ssa/generate_signal.m` and the experiment configuration remain unchanged.

## Run SSA verification

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_ssa_verification.m');
```

This generates the thesis multi-scale synthetic signal, runs the full SSA pipeline, plots singular values and reconstructions, and reports NMSE, correlation, spectral peak error, and cross-component leakage for several window lengths `L`.

## Automatic SSA calibration

Calibrate a fixed automatic grouping rule on the two thesis synthetic signals (`multiscale`, `close`) using a small parameter grid:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_ssa_calibration.m');
```

The calibration script sweeps:

- `r_max = [8, 10, 12]`
- `sv_ratio_thresh = [0.5, 0.6, 0.7, 0.8]`
- `freq_tol = [0.005, 0.01, 0.02]`
- `min_osc_freq = [0.005, 0.01, 0.015]`

and saves:

- `results/ssa_calibration_trials.csv`
- `results/ssa_calibration_summary.csv`
- `results/ssa_calibration_results.mat`

The selected fixed parameter set is chosen from configurations within 5% of the best score, then tie-broken toward simpler settings.
The automatic grouping rule now:

- only accepts adjacent candidate pairs whose mean dominant frequency exceeds `min_osc_freq`
- keeps only the best two non-overlapping oscillatory pairs
- defines the trend group from remaining low-frequency unpaired leading components

## Frozen automatic SSA experiments

After calibration, run the frozen automatic SSA experiments:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_ssa_frozen_experiments.m');
```

This loads the calibrated `best_params`, evaluates the same two thesis synthetic signals on seeds `11:50`, and saves:

- `results/ssa_frozen_trials.csv`
- `results/ssa_frozen_summary.csv`
- `results/ssa_frozen_overall.csv`
- `results/ssa_frozen_results.mat`

## SSD experiments (adaptive window, controlled comparator)

**Singular Spectrum Decomposition (SSD)** here is implemented as a **controlled adaptive matrix comparator** for the thesis, not as a claim of a fully paper-complete standalone SSD benchmark. The current implementation repeatedly embeds the **current residual** in a Hankel matrix, runs SVD, picks the **best single** adjacent oscillatory eigentriple pair, reconstructs it, subtracts it, and repeats. The embedding length `L` is **adaptive** each iteration: the implementation scores all candidates in `L_list` and picks the best.

For comparability with frozen SSA, the current SSD implementation reuses the same calibrated oscillatory-pair feasibility thresholds (`r_max`, `sv_ratio_thresh`, `freq_tol`, `min_osc_freq`) rather than introducing a separate full SSD-specific calibration stage.

The adjacent-pair score used inside SSD is defined explicitly as

`score(r, r+1) = (sigma_r + sigma_{r+1}) * min(sigma_r, sigma_{r+1}) / max(sigma_r, sigma_{r+1}) / (1 + |f_r - f_{r+1}| / max(freq_tol, eps))`

where `sigma_r`, `sigma_{r+1}` are the two singular values and `f_r`, `f_{r+1}` are the dominant frequencies of the two elementary reconstructions. There are no extra weights beyond this formula.

Across `L`, the implementation is equivalent to selecting the global best feasible adjacent pair across all candidate window lengths: for each `L`, it keeps the top-scoring feasible pair, then compares those per-`L` maxima and picks the overall best one.

In other words, SSD is used here to test whether **adaptive iterative decomposition within a matrix Hankel framework** helps before moving on to the main tensor-based Hankelization experiments. Algorithm details, assumptions, and limitations are documented in [ssd/SSD_ALGORITHM.md](ssd/SSD_ALGORITHM.md).

Requires calibrated grouping hyperparameters:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_ssd_experiments.m');
```

This loads `best_params` from `results/ssa_calibration_results.mat`, evaluates the same thesis synthetic signals and SNRs on seeds `11:50`, and saves:

- `results/ssd_trials.csv`
- `results/ssd_summary.csv`
- `results/ssd_overall.csv`
- `results/ssd_results.mat`

Trial rows use `L = NaN` (adaptive `L`); `num_pairs` stores the number of SSD oscillatory modes extracted.

Important interpretation note:

- The current SSD implementation should be read as a **controlled adaptive comparator** alongside SSA, not as a definitive claim about the best possible SSD variant.
- In the `multiscale` synthetic case, the final SSD residual is interpreted as a **trend-dominated remainder**, not as a guaranteed pure trend estimate.
- The SSD `success` flag is only a minimal implementation check (`M >= 2` extracted oscillatory modes for the current synthetic tasks). Actual separation quality must be judged from the reported component metrics (`NMSE`, correlation, peak-frequency error, leakage).
- If SSD underperforms SSA on some signals, that does **not** by itself answer the thesis negatively; the main thesis comparison is still between **single-matrix Hankelization** and later **decimation-based tensor Hankelization**.

### SSA vs SSD comparison

After both frozen SSA and SSD batches exist, aggregate metrics side-by-side (SSA averaged over seeds **and** over `L`; SSD over seeds):

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_compare_ssa_ssd.m');
```

Writes:

- `results/ssa_ssd_comparison_by_condition.csv` — wide table (per `signal_type`, `snr_db`)
- `results/ssa_ssd_comparison_long.csv` — long table with a `method` column (`SSA` / `SSD`)

This comparison is intended to answer a controlled matrix-level question: does an adaptive iterative matrix method help relative to classical SSA on the same synthetic benchmarks before introducing the tensor model?

## Matrix-method sensitivity study

To support the thesis claim that `SSA` depends strongly on window length, while `SSD` reduces some of that dependence but remains a single-resolution matrix method, run:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_matrix_sensitivity_study.m');
```

This script performs a controlled study on the same two synthetic signals and seeds `11:50`:

- `SSA` is run at fixed `L = [200, 400, 600, 1000]`
- `SSD` is run with fixed singleton settings `L_list = [200]`, `[400]`, `[600]`, `[1000]`
- `SSD` is also run with the adaptive setting `L_list = [200, 400, 600, 1000]`

and saves:

- `results/matrix_sensitivity_trials.csv`
- `results/matrix_sensitivity_summary.csv`
- `results/matrix_sensitivity_stability.csv`
- `results/matrix_sensitivity_results.mat`
- `results/matrix_sensitivity_nmse_by_setting.png`
- `results/matrix_sensitivity_rho_by_setting.png`
- `results/matrix_sensitivity_nmse_span.png`

Interpretation:

- `matrix_sensitivity_summary.csv` compares mean metrics for each method under each fixed/adaptive setting
- `matrix_sensitivity_stability.csv` summarizes the span across **fixed** settings only, which makes the dependence on resolution choice easier to compare between `SSA` and `SSD`
- the adaptive SSD setting is included to show that even when `SSD` reduces sensitivity by choosing among candidate `L` values, it still remains a matrix-based single-resolution method at each extraction step

## Close-frequency diagnostic

To diagnose whether the noiseless close-frequency failure comes from automatic grouping or from SSA separability itself, run:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_ssa_close_truth_guided_diagnostic.m');
```

This compares:

- frozen automatic grouping
- truth-guided grouping over adjacent eigentriple pairs

for `L = [200, 400, 600, 1000]` on the noiseless close-frequency signal, and saves:

- `results/ssa_close_truth_guided_diagnostic.csv`
- `results/ssa_close_truth_guided_diagnostic.mat`

## Tensor benchmark experiments

To evaluate the finalized tensor CPD pipeline on the same synthetic benchmarks
in a form that is directly comparable to SSA and SSD, run:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_tensor_benchmark_experiments.m');
```

This uses the tensor pipeline with:

- `reconstruction_mode = 'cpd'` as the main reported inverse path,
- `grouping_mode = 'benchmark'` as the reported benchmark mapping,
- the same synthetic signal families, SNR values, and seeds as the SSA/SSD
  batches,

and saves:

- `results/tensor_benchmark_trials.csv`
- `results/tensor_benchmark_summary.csv`
- `results/tensor_benchmark_best.csv`
- `results/tensor_benchmark_results.mat`

The current tensor benchmark workflow chooses the best setting within the
candidate tensor grid separately for each synthetic signal family and SNR.

To generate the optional RQ2 component-level TensorCPD table from the selected
best settings, run:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_tensor_benchmark_component_metrics.m');
```

This writes:

- `results/tensor_benchmark_component_summary.csv`
- `results/tensor_benchmark_component_trials.csv`
- `results/tensor_benchmark_component_metrics.mat`

The summary reports NMSE and correlation separately for `trend` / `slow` /
`fast` in `multiscale` and `comp1` / `comp2` in `close`.

Interpretation note:

- the reported tensor method is a **CPD-based benchmark method**, not the
  earlier CPD+SVD hybrid diagnostic path;
- CPD first reconstructs the full signal estimate from decimated Hankel tensor
  factors;
- for the synthetic benchmark tables, the final reported components are obtained
  by projecting the CPD reconstruction onto the known benchmark subspaces
  (trend/slow/fast for `multiscale`, two sinusoidal bases for `close`);
- `oracle` and `svd` modes remain in the repository for diagnostics, but the
  benchmark CSVs are produced from the finalized `cpd` + `benchmark` path.

## SSA vs SSD vs Tensor comparison

After running the frozen SSA experiments, SSD experiments, and tensor benchmark
experiments, aggregate them side-by-side with:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_compare_ssa_ssd_tensor.m');
```

Writes:

- `results/ssa_ssd_tensor_comparison_by_condition.csv`
- `results/ssa_ssd_tensor_comparison_long.csv`

## Representative reconstruction plots

To generate a small set of representative signal-space figures for `SSA`,
`SSD`, and `Tensor CPD`, run:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_representative_reconstructions.m');
```

This script regenerates a few representative benchmark cases using the frozen /
best settings already selected by the repository results and exports, for each
case:

- a component-by-component reconstruction comparison across methods
- a signal-flow view showing input, method-specific intermediate signal objects,
  and outputs

The exported figures are written to `results/` with filenames starting with
`representative_`.

## Slow-component recovery results

To generate the focused `multiscale` slow-component comparison across `SSA`,
`SSD`, and `TensorCPD`, run:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_slow_component_recovery.m');
```

This writes:

- `results/slow_component_recovery_trials.csv`
- `results/slow_component_recovery_summary.csv`
- `results/slow_component_recovery_plot_data.csv`
- `results/slow_component_recovery_multiscale_<snr>db_seed<seed>.png`

The summary reports slow-component correlation and RMSE by SNR. The plot overlays
the true slow sinusoid with each method's matched extracted slow component.

## Financial decomposition (RQ3)

The RQ3 financial experiments use the daily EUR/USD reference exchange rate
from the ECB Data Portal and daily S&P 500 closes from Stooq over
`1999-01-04` to `2026-04-02`. Each level series is filtered to this date range,
rows with missing values are removed, and the remaining values are z-score
standardized separately per asset before decomposition.

Run the financial decomposition pipeline with:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_financial_decomposition.m');
```

The script applies:

- `SSA` with `L = [200, 400, 600, 1000]`
- `trend-enabled SSD` with `L_list = [200, 400, 600, 1000]` and trend-first
  extraction enabled before oscillatory-pair extraction
- `TensorCPD` with strides `[1, 2]`, base window `200`,
  ranks `[1, 2, 3, 4, 6, 8]`, and both `auto` and `none` grouping

It writes:

- `results/financial_<asset>_ssa.csv`
- `results/financial_<asset>_ssd.csv`
- `results/financial_<asset>_tensor.csv`
- `results/financial_diagnostics.csv`
- `results/financial_ssd_variant_comparison.csv`

The SSD comparison table explicitly reports the initial oscillatory-only SSD
beside the `trend-enabled SSD` for each financial series. It includes
`norm(x_z)`, `norm(reconstruction)`, `norm(residual)`, `relative_error`, `fit`,
and `max(abs(reconstruction))` sanity-check diagnostics.

After the CSVs exist, generate the qualitative RQ3 figures with:

```matlab
cd('/path/to/thesis-implementation');
run('scripts/run_financial_visualizations.m');
```

For each asset, this exports:

- `financial_<asset>_ssa_overview.png`
- `financial_<asset>_ssd_overview.png`
- `financial_<asset>_tensor_overview.png`
- `financial_<asset>_methods_compare.png`

## Layout

| Path | Role |
|------|------|
| `ssa/generate_signal.m` | Synthetic signals (multi-scale, close-frequency) |
| `ssa/embed_hankel.m` | Trajectory (Hankel) matrix |
| `ssa/diagonal_average.m` | Hankelization / diagonal averaging |
| `ssa/ssa_decompose.m` | Embed → SVD → group → reconstruct |
| `ssa/match_components.m` | Greedy matching by \|correlation\| |
| `ssa/compute_metrics.m` | NMSE, Pearson, peak frequency error, leakage |
| `ssa/ssa_precompute.m` | Precompute SSA factors and leading elementary components |
| `ssa/ssa_reconstruct_groups.m` | Reconstruct grouped components from precomputed SSA factors |
| `ssa/estimate_dominant_frequency.m` | Dominant FFT peak helper for automatic grouping |
| `ssa/auto_group_ssa.m` | Fixed automatic grouping rule over leading eigentriples |
| `ssa/evaluate_auto_group_ssa.m` | Compare automatic grouping against synthetic ground truth |
| `ssa/select_ssa_params.m` | Choose a robust fixed parameter set from calibration |
| `scripts/run_ssa_verification.m` | End-to-end demo and sensitivity study |
| `scripts/run_ssa_calibration.m` | Parameter-grid calibration for automatic SSA grouping |
| `scripts/run_ssa_frozen_experiments.m` | Frozen automatic SSA experiments after calibration |
| `scripts/run_ssa_close_truth_guided_diagnostic.m` | Diagnose close-frequency noiseless failure with truth-guided grouping |
| `ssd/SSD_ALGORITHM.md` | SSD algorithm spec (iterative SSA on residual, adaptive `L`) |
| `ssd/ssd_merge_params.m` | Merge `best_params` with SSD defaults (`L_list`, `max_modes`, …) |
| `ssd/ssd_decompose.m` | SSD main loop |
| `ssd/ssd_select_window.m` | Choose `L` and best oscillatory pair for current residual |
| `ssd/ssd_list_oscillatory_pairs.m` | Rank feasible adjacent pairs (same tests as `auto_group_ssa`) |
| `ssd/evaluate_ssd.m` | Metrics vs synthetic ground truth (`[modes, residual]` vs truth columns) |
| `scripts/run_ssd_experiments.m` | Monte Carlo SSD experiments after calibration |
| `scripts/run_compare_ssa_ssd.m` | Aggregate SSA vs SSD summaries |
| `tensor/evaluate_tensor_cpd.m` | Tensor CPD evaluation with `cpd` / `svd` / `blockwise` reconstruction modes |
| `tensor/tensor_reconstruct_via_cpd.m` | Finalized pure-CPD inverse mapping used for benchmark experiments |
| `tensor/group_cp_terms_benchmark.m` | Benchmark-oriented grouping for the thesis synthetic signals |
| `tensor/tensor_reconstruct_via_svd.m` | Hybrid CPD+SVD diagnostic reconstruction path |
| `scripts/run_tensor_benchmark_experiments.m` | Tensor Monte Carlo benchmark experiments |
| `scripts/run_compare_ssa_ssd_tensor.m` | Aggregate SSA vs SSD vs Tensor summaries |
| `scripts/run_representative_reconstructions.m` | Representative signal-space visualizations for SSA, SSD, and Tensor |
| `scripts/run_slow_component_recovery.m` | Slow-component correlation, RMSE, and true-vs-extracted recovery plot |
| `scripts/run_financial_decomposition.m` | RQ3 financial decomposition driver for EUR/USD and S&P 500 |
| `scripts/run_financial_visualizations.m` | Qualitative RQ3 financial overview and method-comparison figures |
| `scripts/run_matrix_sensitivity_study.m` | Controlled sensitivity study for SSA and SSD matrix methods |

## References

Thesis: *Exploring Hankel Tensorization Strategies for Time Series Decomposition* (Omer Ege Ozmen, Maastricht University).
