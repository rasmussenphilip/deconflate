# deconflate (development version)

Work towards 0.2.0.

## Culling hazard ratios
* `cm_impacts()` accepts `scale = "hazard_ratio"`. Hazard ratios are adjusted under a multiplicative (Cox-type) model:
  * `method = "global"` solves the model exactly over the distribution of disease combinations, so that the adjusted hazard ratios reproduce the raw ones;
  * `method = "simultaneous"` uses the first-order (log-linear) version;
  * `method = "published"` reproduces Rasmussen et al. (2024) (`HR - 1` with eq. 16).
* New `attributable_risk()`: the culling (or mortality) attributable to disease, without counting an animal with several diseases more than once, allocated to diseases by Shapley values. `summary()`, `contribution_table()` and `compare_methods()` use it to value hazard-ratio outcomes.
* `example_global_dairy()` and `sampler_global_dairy()` gain `culling_scale = c("excess_hr", "hazard_ratio")`.
* New vignette: "Culling and hazard ratios".

## Reading your own data
* New `cm_read_inputs()`: builds a model (and a Monte Carlo sampler) from CSV files or data frames, after checking every table and reporting all problems at once. `cm_check_inputs()` returns the problems without stopping.
* New `cm_template()` (example CSV files) and `cm_dist_table()` (distributions from an uncertainty table).
* Example input files installed with the package: `system.file("extdata", "global_dairy_2024", package = "deconflate")` (the 2024 analysis inputs, with distributions) and `"example_with_errors"` (deliberate mistakes).
* New vignette: "Using your own data".

## Comparing methods
* **Breaking:** `compare_methods()` now returns a `cm_comparison` object with side-by-side tables of adjusted impacts, relative changes, sign changes, totals and (with `economics`) gaps and values per method, and records methods that fail. The adjusted impacts are in `$impacts`.
* `cm_monte_carlo()` accepts several methods (`method = c("published", "simultaneous")`) and applies them to the same draws; `compare_methods()` compares their summaries.

## Monte Carlo stability
* `summary.cm_mc()` gains a `method` column, a trimmed mean and stability diagnostics (`rel_mcse`, `tail_share`, `stability`), and prints suggestions when estimates look unstable. It distinguishes means that do not exist (the published approximation dividing by a quantity that changes sign) from heavy tails and imprecision.
* New `cm_diagnose()` (flagged estimates with suggestions) and `cm_suggest_proposal()` (a defensive-mixture importance-sampling proposal).
* `cm_monte_carlo()` gains `sampling = "lhs"` (Latin hypercube sampling) and `proposal` (importance sampling). Importance weights carry through `cm_scenario()` and `cm_reweight()`.

## Other changes
* `attribute_burden()` and `productivity_gap()` skip hazard-ratio outcomes (use `attributable_risk()`).
* `plot.cm_mc()` gains a `method` argument and plots hazard ratios on their own scale.
* Added `inst/validation/reference_v02.py`.

# deconflate 0.1.0

## Reproduction review

Checked against the published 2024 analysis code (1st revision) and the 2022 corrigendum.

* **Breaking:** `example_global_dairy()` and `sampler_global_dairy()` gain `inputs = c("analysis", "tables")` and `culling = TRUE`. The default, `"analysis"`, uses the inputs of the published analysis code, which reproduce Table 5 of Rasmussen et al. (2024):
  * fixed disease probabilities (`1 - exp(-incidence)` at the unrounded global means; subclinical mastitis entered unconverted, as in the analysis);
  * unrounded impact parameters;
  * a culling outcome entered as hazard ratio minus 1 (with the analysis's rescaled standard deviations in the sampler).

  `inputs = "tables", culling = FALSE` gives the previous behaviour (Tables 2-4 as printed).
* `adjusted_hr()` gains `method = "excess_hr"` (adjusted HR - 1, plus 1), for which `conversion` is not needed.
* `hr_to_risk()` and `hr_conversion()` gain `method = "overall_odds"`, the excess culling risk used for the 2024 losses: `HR * r / (HR * r + 1 - r) - r`.
* The reproduction vignette now has a comparison with Table 9 of Rasmussen et al. (2022), notes on the yield and fertility-allocation differences and on the corrigendum (which changes Table 7 only), and an updated 2024 section with long-run Monte Carlo means for yield and culling.
* Added `inst/validation/reference_2024_analysis.py`.

## 0.0.0.9001

Builds the remaining v0.1 features.

### Reproduction of the published analyses
* Added the UK dairy model of Rasmussen et al. (2022): `example_uk_dairy_2022()` and `uk_dairy_2022_economics()`. It reproduces the fertility and culling results of Tables 8-10 and documents why the yield results differ (the printed Table 4 differs from the inputs implied by Tables 8-10).
* Added `sampler_global_dairy()` with the input distributions of Rasmussen et al. (2024). Monte Carlo means of the adjusted yield impacts reproduce Table 5. The central value of the PERT distributions in Tables 2-4 is the mode.

### Culling and valuation
* Added `hr_conversion()`, `as_impacts()` and `adjusted_hr()`, which take culling hazard ratios through to adjusted hazard ratios. They support proportional-hazards and published (odds-ratio) conversions.
* `productivity_gap()` now supports impacts on the absolute scale.
* Added `value_losses()` to value gaps in money, and `combine_impacts()` to stack impact tables.

### Uncertainty
* Added input distributions: `dist_fixed()`, `dist_normal()` (optionally truncated), `dist_lognormal()`, `dist_lognormal_ci()`, `dist_beta()`, `dist_pert()`, `dist_pert_mean()`, `dist_uniform()` and `dist_mixture()`.
* Added `cm_sampler()`, which draws inputs on their entry scale. It supports optional correlation of impacts across outcomes (Gaussian copula).
* `cm_monte_carlo()` now:
  * accepts `economics` (observed means and unit values, fixed or uncertain) and records gaps and losses for each draw;
  * counts sampler failures as rejected draws;
  * takes `seed` and `progress` arguments.
* Added `cm_scenario()` for scenario reweighting from distributions. `summary.cm_mc()` gains `what = "loss"`, `"total"` and `"rejections"`, and reports Monte Carlo standard errors.
* Parameter keys are now `prob:<id>`, `assoc:<d1>:<d2>`, `impact:<outcome>:<disease>` and `inter:<outcome>:<d1>:<d2>`.

### Feasibility, attribution and sensitivity
* Added `check_feasibility()`. It runs a triple screen without extra packages and an exact linear-programming check when `lpSolve` is installed.
* Added `shapley_by_cell()` with `loss_additive()` and `loss_multiplicative()`, for Shapley attribution of general loss functions.
* Added sensitivity tools: `screen_associations()`, `screen_interactions()`, `sensitivity_oat()` and `compare_scenarios()`. Also added `set_association()` and `set_interaction()` for editing models.

### Reporting
* Added `contribution_table()`, `summary()` for adjustment results, and base-graphics plots: `plot()` for results, Monte Carlo runs, screens and one-at-a-time sensitivity, plus `plot_burden()`.
* Added four vignettes and a pkgdown configuration.

### Other changes
* `pair_tables()` is vectorised, which makes Monte Carlo runs faster.
* Hazard ratios of exactly 1 now give exactly zero excess risk.

## 0.0.0.9000

* Initial skeleton: input constructors, the three adjustment methods, the joint distribution, productivity gaps, Shapley attribution for pairwise interactions, hazard-ratio conversion and Monte Carlo with reweighting.
