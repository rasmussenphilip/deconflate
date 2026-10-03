# deconflate (development version)

This version responds to a second external review of 0.2.0 and adds
threshold searches, a sampled backend for the global model beyond about 20
diseases, uncertainty given in the input tables themselves, and a
five-disease example that uses every feature.

## Breaking changes

* **No separate uncertainty files.** An uncertain value now gets its
  distribution in its own row, in the columns `dist` and `p1`-`p4` of
  `diseases`, `associations`, `three_way`, `impacts_<analysis>` and
  `interactions_<analysis>` tables; rows with an empty `dist` are point
  values, so point values and uncertain values can be mixed in one table.
  `cm_read_inputs()` and `cm_check_inputs()` lose their `uncertainty`
  argument, and files named `uncertainty.csv` or
  `uncertainty_<analysis>.csv` are reported as an error explaining the
  change. The checks report a distribution with missing or invalid
  parameters, parameters without a `dist`, a distribution on a non-numeric
  association measure (`table`, `independent`, `unknown`) and parameters that
  are not numbers (once); a point value outside its distribution's range is
  a note. Distributions given for hazard ratios are ignored with a note.
  The shipped `global_dairy_2024` and `example_with_errors` folders and
  `cm_template()` use the new layout. Every table can also have a free-text
  `note` column.
* **Hazard-ratio estimands are named for the snapshot model and must be
  stated.** `cm_hazard_ratios()` has no default estimand; use
  `estimand = "snapshot_crude"` or `"snapshot_stratified"` (the old `"crude"`
  and `"adjusted"` give an error explaining the change). The documentation
  explains how these differ from Cox hazard ratios estimated over follow-up.
  `hazard_ratios.csv` needs an `estimand` column.
* **Historical conversions are confined to the reproduction functions.**
  `example_uk_dairy_2022()` has the yield and fertility analyses only (the
  2022 culling analysis, hazard ratios treated as odds ratios, is built
  inside `reproduce_rasmussen_2022()`); `example_global_dairy()` and
  `sampler_global_dairy()` lose `culling` (the 2024 HR - 1 analysis is built
  inside `reproduce_rasmussen_2024()`).
* `uk_dairy_2022_economics()` returns the yield and fertility valuations
  (and `additional`); the paper's culling valuation is used inside
  `reproduce_rasmussen_2022()` only.
* Non-positive odds ratios, risk ratios and three-way ratios now raise
  `deconflate_infeasible` (still a `deconflate_error`).
* **`cm_template()` writes one analysis by default** (`impacts.csv`,
  `interactions.csv` and the population files, with distributions in the
  rows of four values),
  which reads into a model and a single sampler with Latin hypercube and
  importance sampling. `cm_template(dir, type = "analyses")` writes the
  multi-analysis example. Hazard ratios are no longer in the template.

## New features

* **A five-disease example that uses every feature:**
  `system.file("extdata", "five_diseases", package = "deconflate")` has all
  probability types and association measures, a three-way term, crude and
  adjusted impacts, interactions, snapshot hazard ratios, and every
  distribution type mixed with point values, in three analyses (yield,
  calving interval, welfare). Its `run_all_features.R` runs every part of
  the package on it; `README.md` describes the files. Reference values are
  in `inst/validation/reference_five_diseases.py`.
* **Three-way ratios can be uncertain:** `cm_sampler()` and
  `cm_batch_sampler()` take `three_way` distributions (keys
  `three:<d1>:<d2>:<d3>`, diseases in any order); in batch runs they are
  shared population inputs.

* **Threshold searches:** `cm_threshold()` varies one input (an association,
  including unknown or unlisted pairs, an interaction, a disease probability,
  a raw impact or a three-way ratio) over a range and finds where a
  conclusion changes: two diseases swap rank (by contribution or adjusted
  impact), an adjusted impact changes sign (a model-implied sign change), the
  aggregate crosses a target, or the aggregate departs from its baseline by
  a given relative amount. The range is scanned before crossings are refined
  by bisection, every crossing is reported, unusable parts of the range
  (infeasible, singular, undefined, unresolved) are listed separately, and a
  jump across a pole or a nearly singular system is reported as a
  discontinuity, never as a threshold. Results just below and above each
  threshold are returned; `plot()` shows the scan. For the global method,
  inputs that do not change the population (impacts and interactions) reuse
  one fitted joint distribution; `inter:` and `three:` inputs need the
  global method. No random numbers are drawn.
* `screen_interactions()` passes `...` to `fit_joint()` (e.g.
  `backend = "sampled"`), and `compare_methods()` for hazard-ratio models
  passes `fit_joint()` arguments through `...`.
* **Sampled backend for the joint distribution:**
  `fit_joint(backend = "sampled")` fits the maximum-entropy model without
  enumerating 2^n combinations: Monte Carlo moment matching with parallel
  Gibbs chains (the conditional interaction parameters are calibrated, not
  set to the marginal log odds ratios), sampling, and by default raking of
  the sample weights to match the marginals and pairwise tables exactly. It
  reports constraint residuals before and after raking, R-hat, effective
  sample sizes and Monte Carlo errors; a fit that does not converge is
  reported as unresolved, not as infeasible. The global method, interaction
  offsets, Shapley allocation, the snapshot hazard model and
  `attributable_risk()` use it unchanged. Pass `backend = "sampled"` through
  `deconflate(..., method = "global")`, or pass the fitted joint.
* The pairwise methods (`"simultaneous"`, `"published"`) never enumerate
  combinations and work for any number of diseases (documented in
  `?deconflate`).

## Corrections (second review)

* **Importance-sampling support:** distributions now record their support as
  a set of intervals; a mixture's support is the union of the supports of
  its components with positive weight. A proposal must cover the input's
  whole support, so a mixture that spans the range but leaves a gap is
  rejected (it previously converged to a biased estimate with a high
  effective sample size). The same check applies to `cm_scenario()`.
* **Undefined results:** one finiteness check is shared by Monte Carlo,
  `compare_methods()`, the sensitivity screens, `compare_scenarios()` and
  `cm_threshold()`. A published result with a non-finite value is reported
  as failed ("undefined"); `compare_methods()` keeps it in `undefined` for
  inspection but not among the estimates. Screens stop if their baseline is
  undefined; `compare_scenarios()` gains a `failed` column.
* **Precision status:** a new stability status, `"insufficient_info"`, is
  used when precision cannot be assessed (fewer than two Latin hypercube
  blocks with positive weight, or a mean of zero with a positive standard
  error); a missing standard error is never reported as `"ok"`. The Latin
  hypercube standard error is the ratio-estimator error over all replicate
  blocks, counting blocks whose draws were all rejected or have zero weight
  (`n_blocks` is kept in the run).
* **Three-way terms:** the documentation now says that they leave additive
  results unchanged only when the pairs are constrained; with unknown pairs
  they can change them.
* `simulate_raw_impacts()` stops when its joint distribution (supplied or
  fitted) has not converged.
* `shapley_by_cell()` caches the loss of each combination by code, which
  makes `attributable_risk()` much faster for 10 or more diseases and works
  with sampled joints.

## Validation

* `inst/validation/reference_v030.py`: thresholds, the support counterexample
  and the reviewer's three-way counterexample.
* New tests: `test-threshold.R`, `test-joint-sampled.R` (against exact
  enumeration for 3 and 12 diseases, and a 24-disease population) and
  `test-review-v03.R` (each item of the second review).

# deconflate 0.2.0

Work towards 0.2.0. This version implements an external review of 0.1.0 and
the planned hazard-ratio, import and Monte Carlo work. It changes the
interface: models written for 0.1.0 need updating (see "Breaking changes").

## Breaking changes

* **One impact vector per analysis, in any units.** `cm_impacts()` no longer
  has `outcome`, `scale` or `direction`. Each analysis adjusts one set of
  additive impacts in the units supplied (kg, percent, days, euros, welfare
  scores); results come back in the same units, with optional `label` and
  `units` metadata. The engine does no unit conversion.
* **Population and analyses.** New `cm_population()` holds the diseases,
  associations and three-way terms. `cm_model(population, impacts,
  interactions)` adds one impact vector; `cm_analyses(population, yield = ...,
  fertility = ...)` bundles several analyses on one population, and
  `deconflate()` adjusts each of them (one joint fit is shared by the global
  method). `cm_model()` still accepts a `cm_diseases` object with
  `associations =` given by name.
* **Estimands.** Each impact declares its estimand: `"crude"` (default) or
  `"adjusted_linear"` (the coefficient of an additive regression adjusted
  for the diseases in `adjusted_for`, or `"all"`). Adjusted coefficients are
  mapped exactly by the population projection of the omitted diseases, with
  an identifiability check. `adjusted_for` without `"adjusted_linear"` is an
  error: the estimand is never inferred from it. The published method is
  limited to crude estimates.
* **Hazard ratios have their own adapter.** `cm_hazard_ratios()`,
  `cm_hr_model()` and `deconflate_hr()` (methods `"snapshot"`, a snapshot
  hazard-multiplier model solved exactly over the joint distribution,
  `"first_order"` and `"published"`), with stratified `adjusted_for` and
  `"all"`. `attributable_risk()` takes their results.
* **Removed:** `combine_impacts()`, `hr_to_risk()`, `excess_to_hr()`,
  `hr_conversion()`, `as_impacts()` and `adjusted_hr()`. The conversions used
  in the papers are kept only to reproduce them: hazard ratios as odds ratios
  and eq. 23 (internal, used by `example_uk_dairy_2022()` and
  `reproduce_rasmussen_2022()`), and HR - 1 adjusted additively (the
  `culling_hr_minus_1` analysis of `example_global_dairy(culling = TRUE)`,
  used by `reproduce_rasmussen_2024()`).
* **Gaps and values are helpers outside the engine.** `productivity_gap(result,
  observed, direction, effect)` and `value_losses(gap, unit_value,
  additional)`; the `valuation` argument of `summary()`,
  `contribution_table()` and `compare_methods()` replaces `economics` and
  must state `direction` and `effect`. Lump-sum costs (e.g. veterinary
  expenditure) are added with `additional`, never inside the adjustment.
* **Contributions.** `deconflate()` returns `contributions`: the closed-form
  Shapley value of each disease (its own term plus half of each interaction
  term it is involved in). They add up to the aggregate; shares are `NA` when
  the aggregate is zero. Gaps are attributed without dividing by the
  aggregate.
* **Covariate-adjusted associations** are rejected unless
  `adjusted_associations = "use_as_marginal"`.
* **Contingency tables:** empty tables are rejected; the zero-cell correction
  is explicit (`zero_cell = c("haldane", "error")`) and recorded in
  `corrected`.
* **Disease ids** must not contain `|`, `;` or `:`; `"all"` is reserved.
* **Monte Carlo keys** in `cm_sampler()` are `impact:<disease>` and
  `inter:<d1>:<d2>` (input files use `impact:<analysis>:<disease>` and
  `inter:<analysis>:<d1>:<d2>`); the sampler no longer correlates outcomes
  (`outcome_correlation` removed).
* `cm_read_inputs()` reads one impact table per analysis (see below).

## New features

* **Reading your own data:** `impacts_<analysis>.csv` files (or a named list
  of data frames), `interactions_<analysis>.csv`, `three_way.csv`,
  `hazard_ratios.csv`, and `uncertainty.csv` / `uncertainty_<analysis>.csv`.
  The result has `population`, `analyses`, `model` (one analysis),
  `hr_model` and a sampler (a batch sampler for several analyses). Every
  problem is reported with its table, row and column; `cm_read_inputs()`
  stops with class `deconflate_input_problems` when there are errors. The
  example folders in `inst/extdata` (`global_dairy_2024`,
  `example_with_errors`) use the new layout.
* `compare_methods()` works on models, analyses, hazard-ratio models and
  Monte Carlo runs. Methods that cannot be run, and valuations that cannot
  be evaluated for a method, are listed in `failed` with the reason.
* **Three-way scenarios:** `cm_three_way()` sets the ratio of conditional
  odds ratios for a triple (a log-linear three-way term in the global fit,
  with all pairwise tables still matched); `set_three_way()` and
  `screen_three_way()`.
* **Batch Monte Carlo:** `cm_batch_sampler()` runs several analyses on shared
  draws of the disease probabilities and associations.
* **Reproduction:** `reproduce_rasmussen_2022()` (Tables 8-10) and
  `reproduce_rasmussen_2024()` (Table 5), with the published values beside
  the package's.
* `cm_mc_gap()`: productivity gaps and values per Monte Carlo draw.
* `example_global_dairy_hr()`: the 2024 culling hazard ratios as a
  hazard-ratio model.

## Corrections and checks

* **Feasibility:** the pairwise methods screen every triple of diseases by
  default (`feasibility = c("screen", "lp", "none")`); an invertible
  conflation matrix does not mean the pairs are jointly feasible.
* **Joint distribution:** the IPF rescales every marginal in each sweep and
  recomputes the final residual; a supplied `joint` is validated against the
  model (diseases, probabilities, pairs and three-way terms).
* **Importance sampling:** proposals are checked against the support of each
  input (bounds and point masses), and the Monte Carlo standard error uses
  the self-normalised estimate `sqrt(sum(w^2 (x - mean)^2))`.
* **Latin hypercube sampling** runs in independent replicate blocks
  (`lhs_replicates`), and the standard error comes from the spread of the
  block means.
* **Stability statuses:** `no_mean` is replaced by `possible_pole`, flagged
  only when the published denominator `m + c` changes sign within the
  sampled inputs and differs in sign from `m`; the removable case (`c = 0`,
  where the formula reduces to `m`) is not flagged.
* A draw with a non-finite result is rejected for every method.
* **Conditions have classes:** `deconflate_infeasible`, `deconflate_singular`,
  `deconflate_nonconvergence`, `deconflate_unsupported` and
  `deconflate_nonfinite`.
* `shapley_by_cell()` reports the probability and loss of any skipped
  combinations, and `attributable_risk()` the unallocated part; the
  allocated and unallocated parts add up to the total.
* Sensitivity screens report failed scenarios in a `failed` column instead
  of dropping them, and reject unknown pairs or triples.

## Validation

* New Python references: `inst/validation/reference_v020.py` (estimands,
  interactions, feasibility, three-way terms) and `reference_v020_tests.py`.
* Tests cover the reviewer's list: known-truth recovery for crude, adjusted
  and fully adjusted estimands, with interactions, in arbitrary units and
  after rescaling; the projection counterexample; zero, cancelling, singular
  and non-finite cases; infeasible pairs with a well-conditioned matrix; a
  stale joint; removable poles; importance-sampling support and standard
  errors; repeated-run calibration of the standard errors; separate impact
  files on one population; allocation completeness; and the three-way term
  of the fitted joint.

## Known limitations

* The adjustment is exact for additive impacts. Adjusted estimands other than
  coefficients of additive regressions (e.g. matched or propensity-score
  estimates, or coefficients of non-linear models) are not supported.
* The probabilities and associations must describe the population the impact
  estimates come from.
* Pairwise associations do not identify three-way structure; the global
  method assumes maximum entropy unless three-way scenarios are given.
* The global method, the snapshot hazard model and `attributable_risk()`
  enumerate all 2^n disease combinations (about 20 diseases at most).
* The triple screen is a necessary condition for joint feasibility only; the
  exact check needs `lpSolve`.

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
