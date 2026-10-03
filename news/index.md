# Changelog

## deconflate (development version)

## deconflate 0.2.0

Work towards 0.2.0. This version implements an external review of 0.1.0
and the planned hazard-ratio, import and Monte Carlo work. It changes
the interface: models written for 0.1.0 need updating (see “Breaking
changes”).

### Breaking changes

- **One impact vector per analysis, in any units.**
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  no longer has `outcome`, `scale` or `direction`. Each analysis adjusts
  one set of additive impacts in the units supplied (kg, percent, days,
  euros, welfare scores); results come back in the same units, with
  optional `label` and `units` metadata. The engine does no unit
  conversion.
- **Population and analyses.** New
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  holds the diseases, associations and three-way terms.
  `cm_model(population, impacts, interactions)` adds one impact vector;
  `cm_analyses(population, yield = ..., fertility = ...)` bundles
  several analyses on one population, and
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  adjusts each of them (one joint fit is shared by the global method).
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  still accepts a `cm_diseases` object with `associations =` given by
  name.
- **Estimands.** Each impact declares its estimand: `"crude"` (default)
  or `"adjusted_linear"` (the coefficient of an additive regression
  adjusted for the diseases in `adjusted_for`, or `"all"`). Adjusted
  coefficients are mapped exactly by the population projection of the
  omitted diseases, with an identifiability check. `adjusted_for`
  without `"adjusted_linear"` is an error: the estimand is never
  inferred from it. The published method is limited to crude estimates.
- **Hazard ratios have their own adapter.**
  [`cm_hazard_ratios()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hazard_ratios.md),
  [`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md)
  and
  [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
  (methods `"snapshot"`, a snapshot hazard-multiplier model solved
  exactly over the joint distribution, `"first_order"` and
  `"published"`), with stratified `adjusted_for` and `"all"`.
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  takes their results.
- **Removed:** `combine_impacts()`, `hr_to_risk()`, `excess_to_hr()`,
  `hr_conversion()`, `as_impacts()` and `adjusted_hr()`. The conversions
  used in the papers are kept only to reproduce them: hazard ratios as
  odds ratios and eq. 23 (internal, used by
  [`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md)
  and
  [`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)),
  and HR - 1 adjusted additively (the `culling_hr_minus_1` analysis of
  `example_global_dairy(culling = TRUE)`, used by
  [`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)).
- **Gaps and values are helpers outside the engine.**
  `productivity_gap(result, observed, direction, effect)` and
  `value_losses(gap, unit_value, additional)`; the `valuation` argument
  of [`summary()`](https://rdrr.io/r/base/summary.html),
  [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)
  and
  [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
  replaces `economics` and must state `direction` and `effect`. Lump-sum
  costs (e.g. veterinary expenditure) are added with `additional`, never
  inside the adjustment.
- **Contributions.**
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  returns `contributions`: the closed-form Shapley value of each disease
  (its own term plus half of each interaction term it is involved in).
  They add up to the aggregate; shares are `NA` when the aggregate is
  zero. Gaps are attributed without dividing by the aggregate.
- **Covariate-adjusted associations** are rejected unless
  `adjusted_associations = "use_as_marginal"`.
- **Contingency tables:** empty tables are rejected; the zero-cell
  correction is explicit (`zero_cell = c("haldane", "error")`) and
  recorded in `corrected`.
- **Disease ids** must not contain `|`, `;` or `:`; `"all"` is reserved.
- **Monte Carlo keys** in
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
  are `impact:<disease>` and `inter:<d1>:<d2>` (input files use
  `impact:<analysis>:<disease>` and `inter:<analysis>:<d1>:<d2>`); the
  sampler no longer correlates outcomes (`outcome_correlation` removed).
- [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  reads one impact table per analysis (see below).

### New features

- **Reading your own data:** `impacts_<analysis>.csv` files (or a named
  list of data frames), `interactions_<analysis>.csv`, `three_way.csv`,
  `hazard_ratios.csv`, and `uncertainty.csv` /
  `uncertainty_<analysis>.csv`. The result has `population`, `analyses`,
  `model` (one analysis), `hr_model` and a sampler (a batch sampler for
  several analyses). Every problem is reported with its table, row and
  column;
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  stops with class `deconflate_input_problems` when there are errors.
  The example folders in `inst/extdata` (`global_dairy_2024`,
  `example_with_errors`) use the new layout.
- [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
  works on models, analyses, hazard-ratio models and Monte Carlo runs.
  Methods that cannot be run, and valuations that cannot be evaluated
  for a method, are listed in `failed` with the reason.
- **Three-way scenarios:**
  [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md)
  sets the ratio of conditional odds ratios for a triple (a log-linear
  three-way term in the global fit, with all pairwise tables still
  matched);
  [`set_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/set_three_way.md)
  and
  [`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md).
- **Batch Monte Carlo:**
  [`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md)
  runs several analyses on shared draws of the disease probabilities and
  associations.
- **Reproduction:**
  [`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
  (Tables 8-10) and
  [`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
  (Table 5), with the published values beside the package’s.
- [`cm_mc_gap()`](https://rasmussenphilip.github.io/deconflate/reference/cm_mc_gap.md):
  productivity gaps and values per Monte Carlo draw.
- [`example_global_dairy_hr()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy_hr.md):
  the 2024 culling hazard ratios as a hazard-ratio model.

### Corrections and checks

- **Feasibility:** the pairwise methods screen every triple of diseases
  by default (`feasibility = c("screen", "lp", "none")`); an invertible
  conflation matrix does not mean the pairs are jointly feasible.
- **Joint distribution:** the IPF rescales every marginal in each sweep
  and recomputes the final residual; a supplied `joint` is validated
  against the model (diseases, probabilities, pairs and three-way
  terms).
- **Importance sampling:** proposals are checked against the support of
  each input (bounds and point masses), and the Monte Carlo standard
  error uses the self-normalised estimate `sqrt(sum(w^2 (x - mean)^2))`.
- **Latin hypercube sampling** runs in independent replicate blocks
  (`lhs_replicates`), and the standard error comes from the spread of
  the block means.
- **Stability statuses:** `no_mean` is replaced by `possible_pole`,
  flagged only when the published denominator `m + c` changes sign
  within the sampled inputs and differs in sign from `m`; the removable
  case (`c = 0`, where the formula reduces to `m`) is not flagged.
- A draw with a non-finite result is rejected for every method.
- **Conditions have classes:** `deconflate_infeasible`,
  `deconflate_singular`, `deconflate_nonconvergence`,
  `deconflate_unsupported` and `deconflate_nonfinite`.
- [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
  reports the probability and loss of any skipped combinations, and
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  the unallocated part; the allocated and unallocated parts add up to
  the total.
- Sensitivity screens report failed scenarios in a `failed` column
  instead of dropping them, and reject unknown pairs or triples.

### Validation

- New Python references: `inst/validation/reference_v020.py` (estimands,
  interactions, feasibility, three-way terms) and
  `reference_v020_tests.py`.
- Tests cover the reviewer’s list: known-truth recovery for crude,
  adjusted and fully adjusted estimands, with interactions, in arbitrary
  units and after rescaling; the projection counterexample; zero,
  cancelling, singular and non-finite cases; infeasible pairs with a
  well-conditioned matrix; a stale joint; removable poles;
  importance-sampling support and standard errors; repeated-run
  calibration of the standard errors; separate impact files on one
  population; allocation completeness; and the three-way term of the
  fitted joint.

### Known limitations

- The adjustment is exact for additive impacts. Adjusted estimands other
  than coefficients of additive regressions (e.g. matched or
  propensity-score estimates, or coefficients of non-linear models) are
  not supported.
- The probabilities and associations must describe the population the
  impact estimates come from.
- Pairwise associations do not identify three-way structure; the global
  method assumes maximum entropy unless three-way scenarios are given.
- The global method, the snapshot hazard model and
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  enumerate all 2^n disease combinations (about 20 diseases at most).
- The triple screen is a necessary condition for joint feasibility only;
  the exact check needs `lpSolve`.

## deconflate 0.1.0

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
- `adjusted_hr()` gains `method = "excess_hr"` (adjusted HR - 1, plus
  1), for which `conversion` is not needed.
- `hr_to_risk()` and `hr_conversion()` gain `method = "overall_odds"`,
  the excess culling risk used for the 2024 losses:
  `HR * r / (HR * r + 1 - r) - r`.
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

- Added `hr_conversion()`, `as_impacts()` and `adjusted_hr()`, which
  take culling hazard ratios through to adjusted hazard ratios. They
  support proportional-hazards and published (odds-ratio) conversions.
- [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
  now supports impacts on the absolute scale.
- Added
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)
  to value gaps in money, and `combine_impacts()` to stack impact
  tables.

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
