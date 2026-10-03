# Changelog

## deconflate (development version)

### Reproduction review

Checked against the published 2024 analysis code (1st revision) and the
2022 corrigendum.

- **Breaking:**
  [`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)
  and
  [`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
  gain `inputs = c("analysis", "tables")` and `culling = TRUE`. The
  default, `"analysis"`, uses the inputs of the published analysis code,
  which reproduce Table 5 of Rasmussen et al. (2024):
  - fixed disease probabilities (`1 - exp(-incidence)` at the unrounded
    global means; subclinical mastitis entered unconverted, as in the
    analysis);
  - unrounded impact parameters;
  - a culling outcome entered as hazard ratio minus 1 (with the
    analysis’s rescaled standard deviations in the sampler).

  `inputs = "tables", culling = FALSE` gives the previous behaviour
  (Tables 2-4 as printed).
- [`adjusted_hr()`](https://rasmussenphilip.github.io/deconflate/reference/adjusted_hr.md)
  gains `method = "excess_hr"` (adjusted HR - 1, plus 1), for which
  `conversion` is not needed.
- [`hr_to_risk()`](https://rasmussenphilip.github.io/deconflate/reference/hr_to_risk.md)
  and
  [`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md)
  gain `method = "overall_odds"`, the excess culling risk used for the
  2024 losses: `HR * r / (HR * r + 1 - r) - r`.
- The reproduction vignette now has a comparison with Table 9 of
  Rasmussen et al. (2022), notes on the yield and fertility-allocation
  differences and on the corrigendum (which changes Table 7 only), and
  an updated 2024 section with long-run Monte Carlo means for yield and
  culling.
- Added `inst/validation/reference_2024_analysis.py`.

### 0.0.0.9001

Builds the remaining v0.1 features.

#### Reproduction of the published analyses

- Added the UK dairy model of Rasmussen et al. (2022):
  [`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md)
  and
  [`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md).
  It reproduces the fertility and culling results of Tables 8-10 and
  documents why the yield results differ (the printed Table 4 differs
  from the inputs implied by Tables 8-10).
- Added
  [`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
  with the input distributions of Rasmussen et al. (2024). Monte Carlo
  means of the adjusted yield impacts reproduce Table 5. The central
  value of the PERT distributions in Tables 2-4 is the mode.

#### Culling and valuation

- Added
  [`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md),
  [`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md)
  and
  [`adjusted_hr()`](https://rasmussenphilip.github.io/deconflate/reference/adjusted_hr.md),
  which take culling hazard ratios through to adjusted hazard ratios.
  They support proportional-hazards and published (odds-ratio)
  conversions.
- [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
  now supports impacts on the absolute scale.
- Added
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)
  to value gaps in money, and
  [`combine_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/combine_impacts.md)
  to stack impact tables.

#### Uncertainty

- Added input distributions:
  [`dist_fixed()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
  [`dist_normal()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  (optionally truncated),
  [`dist_lognormal()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
  [`dist_lognormal_ci()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
  [`dist_beta()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
  [`dist_pert()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
  [`dist_pert_mean()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
  [`dist_uniform()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  and
  [`dist_mixture()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md).
- Added
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md),
  which draws inputs on their entry scale. It supports optional
  correlation of impacts across outcomes (Gaussian copula).
- [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  now:
  - accepts `economics` (observed means and unit values, fixed or
    uncertain) and records gaps and losses for each draw;
  - counts sampler failures as rejected draws;
  - takes `seed` and `progress` arguments.
- Added
  [`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
  for scenario reweighting from distributions.
  [`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md)
  gains `what = "loss"`, `"total"` and `"rejections"`, and reports Monte
  Carlo standard errors.
- Parameter keys are now `prob:<id>`, `assoc:<d1>:<d2>`,
  `impact:<outcome>:<disease>` and `inter:<outcome>:<d1>:<d2>`.

#### Feasibility, attribution and sensitivity

- Added
  [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md).
  It runs a triple screen without extra packages and an exact
  linear-programming check when `lpSolve` is installed.
- Added
  [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
  with
  [`loss_additive()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md)
  and
  [`loss_multiplicative()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md),
  for Shapley attribution of general loss functions.
- Added sensitivity tools:
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
  [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md),
  [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md)
  and
  [`compare_scenarios()`](https://rasmussenphilip.github.io/deconflate/reference/compare_scenarios.md).
  Also added
  [`set_association()`](https://rasmussenphilip.github.io/deconflate/reference/set_association.md)
  and
  [`set_interaction()`](https://rasmussenphilip.github.io/deconflate/reference/set_interaction.md)
  for editing models.

#### Reporting

- Added
  [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md),
  [`summary()`](https://rdrr.io/r/base/summary.html) for adjustment
  results, and base-graphics plots:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for results,
  Monte Carlo runs, screens and one-at-a-time sensitivity, plus
  [`plot_burden()`](https://rasmussenphilip.github.io/deconflate/reference/plots.md).
- Added four vignettes and a pkgdown configuration.

#### Other changes

- [`pair_tables()`](https://rasmussenphilip.github.io/deconflate/reference/pair_tables.md)
  is vectorised, which makes Monte Carlo runs faster.
- Hazard ratios of exactly 1 now give exactly zero excess risk.

### 0.0.0.9000

- Initial skeleton: input constructors, the three adjustment methods,
  the joint distribution, productivity gaps, Shapley attribution for
  pairwise interactions, hazard-ratio conversion and Monte Carlo with
  reweighting.
