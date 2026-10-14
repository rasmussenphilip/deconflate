# deconflate (development version)

Event impacts now have the same results and tools as additive impacts
(ROADMAP v0.6, item 1; design note section 15). Reference values come from
an independent Python implementation
(`inst/validation/reference_v06_event_parity.py`). Everything that is not
the package's one model has been removed (item 2), and the sensitivity tools
are replaced by one function, `cm_sensitivity()` (item 3).

## New features: sensitivity

* **`cm_sensitivity()`** takes a `deconflate()` result (or a model) and
  varies one input at a time to a low and a high value, reporting the
  aggregate (for event impacts, the attributable risk), each disease's
  adjusted impact and how far the ranking moves. Entered inputs and inputs
  you did not enter are in one table, so missing estimates can be compared
  with the uncertainty of the ones you have:
  * entered inputs (`prob:<disease>`, `associations:<row>`,
    `impacts:<row>`, `risk`): one SD down and up when they have a
    distribution, otherwise a relative `change` (20% by default);
  * unknown inputs: pairs that no row mentions (`pair:<d1>:<d2>`, odds ratio
    0.5 to 3), interactions of pairs without a combination row
    (`interaction:<d1>;<d2>`, plus or minus half the smaller adjusted
    impact; hazard ratio 0.67 to 1.5 for event impacts) and, on request,
    three-way associations (`or_ratio:<d1>:<d2>|<d3>`, 0.5 to 2);
  * any range can be set with `ranges`; `keys` limits the inputs;
  * with draws, each uncertain input's share of the variance of the total
    (squared standardised regression coefficients on the normal scores,
    with R^2).

  `print()` shows the 10 largest swings (`print(x, n = Inf)` all), and
  `plot()` a tornado plot. Reference values come from an independent Python
  implementation (`inst/validation/reference_v06_sensitivity.py`).
* **Rank probabilities:** `deconflate()` results with draws have `$ranks`:
  how often each disease ranks 1st, 2nd, ... by contribution, and its mean
  rank; `summary()` shows them.
* **`cm_threshold()`** also takes a `deconflate()` result (its settings are
  used), and uses the same input keys as `cm_sensitivity()`.
* The diagnosis (`diagnose = TRUE`) runs on the same engine as
  `cm_sensitivity()`.

## Breaking changes: sensitivity tools

* `sensitivity_oat()`, `screen_associations()`, `screen_interactions()` and
  `screen_three_way()` are removed: use `cm_sensitivity()`.
* `set_association()`, `set_interaction()` and `set_three_way()` are no
  longer exported: add rows to the input tables instead (an interaction is
  a combination row adjusted for `all`; a three-way association an
  `OR_ratio` row), or vary them with `cm_sensitivity()`.
* `cm_threshold()` input keys: `"assoc:d1:d2"` is now `"pair:d1:d2"`,
  `"inter:d1:d2"` is `"interaction:d1;d2"`, `"three:d1:d2:d3"` is
  `"or_ratio:d1:d2|d3"` and `"impact:<disease>"` is `"impacts:<row>"`; the
  old keys give an error that says what to write. (Distribution keys in
  `cm_model()` are unchanged.)

## Breaking changes: one model only

* **Removed:** `compare_methods()`, the proportional approximation of
  Rasmussen et al. (2022, eq. 16) and the HR - 1 conversion of 2024, the
  first-order event method, `reproduce_rasmussen_2022()`,
  `reproduce_rasmussen_2024()`, `example_uk_dairy_2022()` and
  `vignette("reproducing-published")`, with the code that existed only for
  them (the Monte Carlo and batch sampler of the 2024 reproduction, and the
  historical hazard-ratio, gap and valuation conversions). The v0.5.0
  release is the record of the published calculations.
* **Removed helpers:** `attribute_burden()` (use `$contributions`, or
  `$attributable$by_disease` for event impacts, or `contribution_table()`),
  `plot_burden()`, `loss_additive()`, `loss_multiplicative()`,
  `combination_probs()` and `excess_matrix()` (both can be computed from
  the `$cells` and `$prob` of `fit_joint()`), and `cm_dist_table()`
  (distributions go in the input tables' `dist` columns).
* **No longer exported** (still used inside the model): `or_to_joint()`,
  `joint_to_or()`, `association_to_joint()`, `pair_tables()` and
  `shapley_by_cell()`.
* The error class `deconflate_unknown_pairs` is gone (no calculation needs
  a crude pair row for every pair any more), and the draw summaries no
  longer have the `possible_pole` stability label (it belonged to the
  published approximation).

## Breaking changes: event impacts

* **Raw totals have one definition for both models:** the result if the
  diseases were unrelated, that is, the same model with the same
  prevalences and impact rows but no association rows. `$totals$raw_sum`
  (the prevalence-weighted sum of the raw impacts) is renamed
  `$totals$raw_total` and follows this definition, so it is now also
  defined when a disease has several rows, subgroup rows or combination
  rows (above 20 diseases, only for single-disease rows; otherwise `NA`).
  For single-disease rows it is the same number as before. The same name is
  used in `summary()` and the draws.
* **Event results:** `change` and the sign check are on the measure of the
  term's raw row (HR, rate ratio, RR, OR or RD), so they are now defined for
  every row measure, not only hazard and rate ratios. The sign-change
  warning reads "on the other side of 1 (of 0 for risk differences)". The
  diagnosis (`diagnose = TRUE`) compares on the same measure, so risk
  difference rows are now also flagged when shrunk by more than 90%.
* `$attributable$by_disease` has the columns `disease`, `main`,
  `interaction`, `attributable` and `share`, and `$attributable$summary`
  gains `raw_attributable`.

## New features

* **Raw attributable risk:** event results report the attributable risk if
  the diseases were unrelated (`raw_attributable`) beside the adjusted one,
  in the print and `summary()` (with the reduction).
* **Own part and interaction share** of each disease's attributable risk
  (`main`, `interaction`), as for additive contributions.
* **Results on the risk scale:** `$adjusted` gains `risk_alone`, `rr_alone`,
  `rd_alone` and `or_alone`: the risk of an animal with only that disease
  (or combination), and its risk ratio, risk difference and odds ratio
  against a disease-free animal. `contribution_table()` shows them (except
  the odds ratio).
* **Local row equations of event impacts** in `$conflation$A`: how each
  row's value (on its measure's scale) changes with each log-hazard term, at
  the adjusted values.
* **`baseline_risk`:** the risk of disease-free animals, as an alternative
  to `overall_risk` (exactly one is given; a number or a distribution).
  The overall risk is then computed by the model. Works in `deconflate()`,
  `cm_threshold()`, `sensitivity_oat()`, the screens and
  `simulate_raw_impacts()`.
* **`simulate_raw_impacts(event_model = TRUE)`:** what each event row would
  report (HR, rate ratio, RR, OR or RD; crude, stratified or for a
  combination) given true hazard ratios and interactions, exactly or from
  `n` sampled animals (exponential event times; HR and rate-ratio rows are
  then estimated as rate ratios, an approximation).
* **Interactions of event impacts in the sensitivity tools:**
  `set_interaction()` adds an HR row adjusted for all;
  `screen_interactions()` and `cm_threshold()` (input `inter:<d1>:<d2>`)
  vary it.

# deconflate 0.5.0

The work planned as version 0.4 is released as 0.5.0 (there is no 0.4.0).
The interface is slimmed to de-conflation itself; everything runs through
one function, `deconflate()`, with one impact table per call; there is one
model and one input layout, in which a row can be about a combination of
diseases or a subgroup, crude or adjusted, and several rows can describe
the same thing; and `deconflate(diagnose = TRUE)` shows which inputs drive
problematic adjusted impacts. See `ROADMAP.md` for what comes next.

## New features: diagnosing the inputs

* `deconflate(diagnose = TRUE)` (off by default) adds `$problems`,
  printed with the result. It flags every disease whose adjusted impact is
  on the other side of zero from its raw impact (for event impacts, of 1),
  or has shrunk by more than 90%, and shows which inputs drive this:
  * for each input with a distribution (prevalences, associations, impacts
    and the overall risk), the change in each flagged impact when the input
    moves up by one SD, and the shift that alone would remove the flag
    (shown within 3 SD);
  * for each input without a distribution, the value at which the flag
    disappears (searched from a quarter to three times its value);
  * flagged impacts that no single input fixes, which points to several
    inputs together or a structural problem (such as units).

  A separate table lists rows that disagree with other rows for the same
  quantity (standardised residual of 2 or more), with the adjusted impact
  that changes most without each. Works for additive and event impacts
  (see `?deconflate`, "Diagnosing the inputs"). Reference values come from an
  independent Python implementation
  (`inst/validation/reference_v05_input_check.py`).
* The sign-change warning points to `diagnose = TRUE`.
* Not added, after a prototype (ROADMAP, item 6): an `anchor` argument,
  constraints on adjusted impacts, MCMC and importance sampling. A negative
  adjusted impact means the inputs contradict each other; a lower bound of
  zero only parks it at zero, and narrowing the distributions makes it more
  certain. Such conflicts are fixed in the inputs, which the check helps to
  find.

## Errors and warnings

* Errors and warnings about particular inputs hold them as fields, so code
  can use them without reading the message (see `?deconflate`, "Errors and
  warnings"): `table` and `row` or `rows` (with `residuals` for
  contradicting impact rows), `terms`, `triples`, `pairs` and `residual`.
  The error of `cm_read_inputs()` (class `deconflate_input_problems`) holds
  every problem found in `$problems`, and the warning of an incomplete
  Shapley allocation holds `skipped_mass`, `skipped_loss` and `total`.
* Every error the package raises itself has a class beside
  `deconflate_error`. The new class `deconflate_invalid_input` covers
  arguments, tables, rows and input values that are not valid (an unknown
  disease, a negative standard deviation, a missing table). Limits of the
  package (too many diseases for the exact backend, the LP check or the
  cell-wise Shapley allocation; the LP check without lpSolve) are
  `deconflate_unsupported`. When every method of `compare_methods()` fails,
  or the baseline of a screen or a threshold search cannot be used, the
  error has the class of the failure that caused it. Rejected draws and
  failed scenarios caused by invalid inputs are labelled `"invalid_input"`
  instead of `"error"`.
* Failures inside base R's root finding and linear solves are classed
  errors, so draws, `compare_methods()` and the screens handle them like
  any other failure instead of stopping. An `OR_ratio` row with a value
  that the other inputs do not allow gave an unclassed error from
  `uniroot()`, which stopped `deconflate()` with draws altogether; it is now
  `deconflate_infeasible` (the draw is rejected). The baseline hazard of
  non-finite hazard multipliers is `deconflate_nonfinite`, and a singular
  system in the first-order event approximation of `compare_methods()` is
  `deconflate_singular`.
* Messages are classed conditions (`deconflate_message`), so that code can
  hide them by class. The notes of `cm_read_inputs()` are also
  `deconflate_input_notes` and hold the notes table in `$problems`. The
  printed text is unchanged.
* An invalid choice for an argument (e.g. `sampling = "x"` in
  `deconflate()`, `pool`, `backend`, `outcome`, `conclusion`) is an error of
  class `deconflate_invalid_input` that names the argument, instead of the
  plain R error of `match.arg()`. Partial matching works as before. The
  problems found by `cm_read_inputs()` (class `deconflate_input_problems`)
  are also `deconflate_invalid_input`.

## Breaking changes: one model, one input layout

* **One model.** The `method` argument is removed from `deconflate()`,
  `cm_threshold()` and the screens (an error says so). `deconflate()`
  always gives the answer of the maximum-entropy model that matches every
  association row and reproduces every impact row. When the pairwise 2x2
  tables give the same answer (every association row a crude pair row,
  every pair covered, one-disease impact rows, additive impacts), they are
  used directly; up to 20 diseases the joint distribution is still fitted
  once, to check that the rows can hold together. `$diagnostics$method`
  says `"pairwise shortcut"` or `"full model"`. The `feasibility` argument
  is removed (use `check_feasibility()` for the separate checks).
  `compare_methods()` sets the model (`"deconflate"`) beside the published
  approximation (`"published"`, and for event impacts `"first_order"`).
* **One input layout.** Three tables, with the same core columns:
  - `diseases`: `disease`, `measure` (`prevalence`, `probability`,
    `incidence_rate`), `value`, and `label`, `time_horizon`,
    `reference_population`;
  - `associations`: `disease`, `associated_with`, `among`, `adjusted_for`,
    `measure`, `value`, `study`;
  - impacts: `disease`, `among`, `adjusted_for`, `measure` (event impacts
    only), `value`, `label`, `units`, `study`;
  plus `dist`, `p1`-`p4`, `source` and `note` everywhere. `disease` and
  `associated_with` take a disease or a combination (`"MAS;SCK"`); `among`
  takes a subgroup (`"SCK"`, `"no SCK"`, `"no others"`). `cm_diseases()`,
  `cm_associations()` and `cm_impacts()` take the same arguments. Old
  columns and arguments (`id`, `type`, `disease1`, `disease2`, `estimand`,
  `adjusted`, `adjusted_associations`) give a message saying what to write
  instead.
* **Association measures:** `OR`, `RR`, `RD`, `phi`, `prevalence` (replaces
  `cond_prob`) and `OR_ratio` (the odds ratio among a subgroup divided by
  that among its complement; replaces three-way terms). A row about a
  combination compares animals with all of its diseases against animals
  with none of them. The three-way and interactions tables, `cm_three_way()`
  and `cm_interactions()` are removed (`set_three_way()` and
  `set_interaction()` now add rows); old files are recognised and reported.
* **Estimands are decided by `adjusted_for`.** An impact or association row
  is crude when `adjusted_for` is empty, and adjusted for the diseases it
  names otherwise (`"all"`: every other disease and term). `estimand` is
  removed. Adjusted associations are odds ratios between two diseases,
  matched as logistic-regression coefficients; ids that are not diseases of
  the model (parity, herd) are noted, and a row adjusted only for them is
  used as crude.
* **Distribution keys** are `prob:<disease>`, `associations:<row>` and
  `impacts:<row>`; `assoc:<d1>:<d2>`, `impact:<disease>`, `inter:<d1>:<d2>`
  and `three:<d1>:<d2>:<d3>` still work when they name one row.
  `sensitivity_oat()` reports inputs by these keys.
* **Results:** `$adjusted` has one row per term (diseases, then
  interactions) with `kind`, `raw` (when the term has one row) and `rows`;
  `$conflation$A` has one row per impact row; `raw_sum` is `NA` when a
  disease has several rows. New: `$estimates` (every row with its fitted
  value and standardised residual), `$heterogeneity` and `$pooling`.

## New features: combinations, subgroups and several rows per item

* **Rows about combinations and subgroups.** Association rows can relate
  combinations, hold within subgroups and be adjusted for other diseases;
  the joint distribution matches them all (raking on events). Impact rows
  about a combination add an interaction (additive, or on the log-hazard
  scale for event impacts); rows within subgroups and adjusted rows also
  identify interactions. `deconflate()` checks that the rows determine
  every term and names those they do not. Contributions split each
  interaction equally among its diseases.
* **Several rows per item.** Association rows about the same quantity are
  pooled on the log scale, by random effects (DerSimonian-Laird; default)
  or fixed effect (`pool = "fixed"`), with Q, I^2 and tau^2 in
  `$heterogeneity`; draws keep the central weights and add the
  between-study spread. Impact rows that over-determine the impacts are
  fitted by weighted least squares (weights 1/SD^2; rows without a
  distribution are exact).
* `simulate_raw_impacts()` simulates rows of any kind (`rows`), with
  interactions named by combination (`c("LAM;MAS" = 2)`).
* The five-disease example has two studies of one odds ratio, an
  `OR_ratio` row, several rows for one disease, a subgroup row and
  interaction rows (additive and event). Reference values:
  `inst/validation/reference_v041.py` (independent Python implementation).

## Breaking changes: one function, one impact table

* **One function, one impact table.** `deconflate()` adjusts one table of
  impact estimates; for several outcomes, call it once per table.
  `cm_analyses()` and the multi-analysis results are removed.
* **Named input tables with any file name.** `cm_read_inputs()` and
  `cm_check_inputs()` take each table as its own argument (`diseases`,
  `associations`, `impacts`), as a path to a CSV file with any name or a data frame. The
  `dir` argument and the fixed file names are removed. `deconflate()` needs
  at least one association; the sensitivity tools also run without any.
  `cm_template()` writes one example of each table (`type` is removed).
* **Pairs without an association are unknown by default.** The model fills
  them in from the other associations, and the result lists
  them with the odds ratios the fit gave them (`$unknown_pairs`). An odds
  ratio of 1 states that two diseases are unrelated. The argument
  `missing_associations` and the association measures `"independent"`,
  `"unknown"` and `"table"` (with `n11`-`n00`) are removed; old tables get
  a message saying what to write instead. The published analyses treated
  unlisted pairs as independent; the `reproduce_*()` functions still do.
* With more than 20 diseases the sampled backend of `fit_joint()` is used
  automatically (for crude pair rows).
* **Event impacts run through `deconflate()`.** With `event_model = TRUE`
  and `overall_risk` (a proportion or a distribution), the impact table has
  a `measure` column: hazard ratios, rate ratios, risk ratios, odds ratios
  and risk differences, which can be mixed, all adjusted with the snapshot
  hazard model. The result has the adjusted hazard ratios and the risk
  attributable to disease, with its Shapley allocation. `deconflate_hr()`,
  `cm_hr_model()`, `cm_hazard_ratios()`, `attributable_risk()` and
  `example_global_dairy_hr()` are removed; `example_global_dairy("culling")`
  and `example_uk_dairy_2022("culling")` give the event impacts. In the
  example folders, `hazard_ratios.csv` is replaced by `culling.csv`.
* **Uncertainty is part of `deconflate()`.** When inputs have distributions
  (in the tables, in `cm_model(distributions = )`, or `overall_risk`),
  `deconflate()` runs `n_draws = 1000` draws by default and adds intervals
  to its tables (`$draws` holds the summary, rejections and the seed). With
  no distributions it gives point estimates, with a note. Rejected draws are
  counted and not replaced. `compare_methods()` takes `n_draws` too.
  Removed: `cm_monte_carlo()`, `cm_sampler()`, `cm_batch_sampler()`,
  `sampler_global_dairy()`, `cm_reweight()`, `cm_scenario()`,
  `compare_scenarios()`, `cm_diagnose()`, `cm_suggest_proposal()`,
  importance sampling and `summary.cm_mc()`.
* **The published approximation is no longer a method of `deconflate()`.**
  The approximation of Rasmussen et al. (2022, eq. 16), and the HR - 1
  approach of Rasmussen et al. (2024), are kept for comparison and
  reproduction only: in `compare_methods()` (whose default methods still
  include `"published"`) and `reproduce_rasmussen_2022()` /
  `reproduce_rasmussen_2024()`. `method = "published"` elsewhere gives an
  error of class `deconflate_unsupported`.
* **Valuation is removed from the package.** The package reports adjusted
  impacts and contributions in the units of the impacts. Removed:
  `productivity_gap()`, `value_losses()`, `cm_mc_gap()`,
  `uk_dairy_2022_economics()`, every `valuation` argument and the
  `unit_value` of the attributable risk. The gap and value calculation of
  the 2022 paper moved into `reproduce_rasmussen_2022()`, whose output is
  unchanged; `vignette("deconflate")` shows the calculation in base R.
* **Examples return one model.** `example_global_dairy(outcome, inputs)`
  and `example_uk_dairy_2022(outcome)` return the model of one outcome
  (`"yield"`, `"fertility"` or `"culling"`).
* **Sensitivity tools** (`screen_associations()`, `screen_interactions()`,
  `screen_three_way()`, `sensitivity_oat()`, `cm_threshold()`) adjust each
  scenario as `deconflate()` does, support event impacts (`event_model = TRUE`,
  `overall_risk`), and run on point estimates. `cm_threshold()` reports the
  methods it used and accepts the input `"risk"`.

## Documentation

* All vignettes, the README and the five-disease example (its README and
  `run_all_features.R`) use the new interface. The vignette on hazard ratios
  is now `vignette("event-impacts")`.

# deconflate 0.3.0

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

# deconflate 0.0.0.9001

Builds the remaining v0.1 features.

## Reproduction of the published analyses
* Added the UK dairy model of Rasmussen et al. (2022): `example_uk_dairy_2022()` and `uk_dairy_2022_economics()`. It reproduces the fertility and culling results of Tables 8-10 and documents why the yield results differ (the printed Table 4 differs from the inputs implied by Tables 8-10).
* Added `sampler_global_dairy()` with the input distributions of Rasmussen et al. (2024). Monte Carlo means of the adjusted yield impacts reproduce Table 5. The central value of the PERT distributions in Tables 2-4 is the mode.

## Culling and valuation
* Added `hr_conversion()`, `as_impacts()` and `adjusted_hr()`, which take culling hazard ratios through to adjusted hazard ratios. They support proportional-hazards and published (odds-ratio) conversions.
* `productivity_gap()` now supports impacts on the absolute scale.
* Added `value_losses()` to value gaps in money, and `combine_impacts()` to stack impact tables.

## Uncertainty
* Added input distributions: `dist_fixed()`, `dist_normal()` (optionally truncated), `dist_lognormal()`, `dist_lognormal_ci()`, `dist_beta()`, `dist_pert()`, `dist_pert_mean()`, `dist_uniform()` and `dist_mixture()`.
* Added `cm_sampler()`, which draws inputs on their entry scale. It supports optional correlation of impacts across outcomes (Gaussian copula).
* `cm_monte_carlo()` now:
  * accepts `economics` (observed means and unit values, fixed or uncertain) and records gaps and losses for each draw;
  * counts sampler failures as rejected draws;
  * takes `seed` and `progress` arguments.
* Added `cm_scenario()` for scenario reweighting from distributions. `summary.cm_mc()` gains `what = "loss"`, `"total"` and `"rejections"`, and reports Monte Carlo standard errors.
* Parameter keys are now `prob:<id>`, `assoc:<d1>:<d2>`, `impact:<outcome>:<disease>` and `inter:<outcome>:<d1>:<d2>`.

## Feasibility, attribution and sensitivity
* Added `check_feasibility()`. It runs a triple screen without extra packages and an exact linear-programming check when `lpSolve` is installed.
* Added `shapley_by_cell()` with `loss_additive()` and `loss_multiplicative()`, for Shapley attribution of general loss functions.
* Added sensitivity tools: `screen_associations()`, `screen_interactions()`, `sensitivity_oat()` and `compare_scenarios()`. Also added `set_association()` and `set_interaction()` for editing models.

## Reporting
* Added `contribution_table()`, `summary()` for adjustment results, and base-graphics plots: `plot()` for results, Monte Carlo runs, screens and one-at-a-time sensitivity, plus `plot_burden()`.
* Added four vignettes and a pkgdown configuration.

## Other changes
* `pair_tables()` is vectorised, which makes Monte Carlo runs faster.
* Hazard ratios of exactly 1 now give exactly zero excess risk.

# deconflate 0.0.0.9000

* Initial skeleton: input constructors, the three adjustment methods, the joint distribution, productivity gaps, Shapley attribution for pairwise interactions, hazard-ratio conversion and Monte Carlo with reweighting.
