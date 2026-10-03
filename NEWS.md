# deconflate (development version)

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
