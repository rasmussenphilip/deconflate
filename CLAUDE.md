# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`deconflate` is an R package (base R only) that adjusts ("de-conflates") disease impact estimates for associations between diseases, so impacts can be summed without double counting (Rasmussen et al. 2022, 2024). Version 0.5.0 is released (the work planned as 0.4); the next version is built item by item, following `ROADMAP.md` (v0.6); `NEWS.md` records each item's breaking changes. Read the current ROADMAP item before changing behaviour.

## Commands

```r
devtools::document()                 # regenerate NAMESPACE and man/ from roxygen (markdown on)
devtools::test()                     # full suite
devtools::test(filter = "item4")     # one file: tests/testthat/test-item4.R
devtools::check()                    # CI runs R CMD check with --no-manual --compact-vignettes=gs+qpdf
pkgdown::check_pkgdown()             # every exported topic must appear in _pkgdown.yml `reference:`
```

Tests call internal functions (`disease_cells()`, `all_present()`, `split_ids()`, ...), so run them through `devtools::test()` or after `devtools::load_all()`. Running `testthat::test_file()` after `library(deconflate)` fails, because the internals are not visible.

The test reference values come from independent Python implementations (numpy, scipy). Run them from their own directory, for example `cd inst/validation && python3 reference_v041.py`.

## Constraints that are easy to miss

- **Generated files:** do not edit `NAMESPACE` or `man/` by hand. Change the roxygen in `R/` and run `document()`. A new export also needs an entry in `_pkgdown.yml`.
- **No new dependencies:** `Imports` is base R only (`graphics`, `grDevices`, `stats`, `utils`). A new dependency needs a reason, and the README's Dependencies table must list it.
- **Reference values are independent:** expected values in tests come from the `inst/validation/reference_*.py` scripts, never from the package's own output. When model semantics change, update the Python reference and hard-code what it prints. The `prototype_v04_*.py` scripts are design prototypes; the `reference_*.py` scripts are the independent checks.
- **One model only:** the package has no other methods. The proportional approximation of 2022 (eq. 16), the HR - 1 conversion of 2024, the first-order event method, `compare_methods()` and the `reproduce_*()` functions were removed in v0.6 (the v0.5.0 tag is the record). Do not add alternative methods back; the pairwise shortcut and the sampled backend are ways of computing the same model.
- **Valuation stays out:** the package reports results in the impacts' own units. Productivity gaps and money values are left to the user.
- **Classed errors:** raise errors with `cm_abort(msg, class = ...)` and warnings with `cm_warn()`, using a class from the list above `cm_abort()` in `R/utils.R` (documented for users in `?deconflate-package`, Conventions; keep the two in step). Draws, the diagnosis and the sensitivity tools catch them by class, for example to count rejected draws. Retired arguments, columns and measures get an error or message that says what to write instead (see `retired_args()` and `old_args()`).
- **Finishing a roadmap item:** it updates the tests, vignettes, README, NEWS, the ROADMAP checkboxes, the five-disease example (`inst/extdata/five_diseases/`) and its `run_all_features.R`.
- **Style:** user-facing text and comments use British spelling, short plain sentences and the package's own vocabulary: row, combination, `among`, `adjusted for`, impact, interaction. The model assumes no species, outcome or unit. Printed tables and plots show only the user's disease codes, outcome label and units, although the examples and documentation use dairy cattle.
- **Output tables:** one set of column names for both models (`raw`, `adjusted`, `lower`, `upper`, `change`, `contribution`, `main`, `interaction`, `share`); both result classes have `$contributions`. `print()` is a short overview, `summary()` the full report in sections, `contribution_table()` the export table (fixed columns, full precision). Printed numbers go through `fmt_num()` (3 significant figures, each value on its own; fractions stay fractions) and `print_table()`; every printout fits 80 columns.
- **Input columns:** only `disease` and `value` are required; `note` is the only free-text column and is never used. `label` (diseases), `time_horizon`, `reference_population` and `source` are no longer read (a note, or a message for the arguments). `cm_model()` and `cm_population()` number rows joined with `rbind()` (`number_rows()`).

## Architecture

The package has one model. Every input row becomes a constraint on a maximum-entropy joint distribution of disease combinations. The impacts are the terms that reproduce every impact row. The pipeline:

1. **Inputs** (`inputs.R`, `import.R`): `cm_diseases()`, `cm_associations()`, `cm_impacts()`, `cm_population()` and `cm_model()`; `cm_read_inputs()` takes a CSV path or data frame per table and reports every problem at once (`cm_check_inputs()`). Association and impact rows share the core columns `disease`, `among`, `adjusted_for`, `measure` and `value`; associations also have `associated_with`. A `disease` can be a combination (`"MAS;SCK"`), and `among` a subgroup (`"SCK"`, `"no SCK"`, `"no others"`). Uncertainty goes in the same row (`dist`, `p1`-`p4`). Distribution keys are `prob:<disease>`, `associations:<row>` and `impacts:<row>`; the old `assoc:`, `impact:`, `inter:` and `three:` keys still work when they name one row.
2. **Row specs** (`rows.R`): each row is parsed into a spec: G = `disease`, H = `associated_with`, `among` as present/absent sets (`apos`/`aneg`), and the adjustment set. This file also pools several rows about one quantity (DerSimonian-Laird by default, or fixed effect) and turns association rows into event-probability targets for raking (`fit_events()`).
3. **Joint distribution** (`joint.R`): `fit_joint()` fits the maximum-entropy distribution by IPF/raking over all 2^n cells. Above 20 diseases, a sampled backend replaces it (`fit_joint_sampled()`). `event_joint()` in `hazard.R` is the wrapper the adjustment uses.
4. **Adjustment:**
   - Additive impacts (`adjust.R`): every impact row is linear in the disease impacts `b` and the interactions `delta`, giving `A theta = y`. One row per unknown is solved exactly; more rows go through weighted least squares (weights 1/SD², rows without a distribution exact). `identify_terms()` names the terms the rows do not determine.
   - Event impacts (`hazard.R`, `event_model = TRUE`): the snapshot hazard model with log-hazard terms, `h0` set from `overall_risk` or `baseline_risk` (exactly one; `check_overall_risk()` returns its kind, carried as `plan$risk_kind`), plus the attributable risk (`$attributable$summary`) and its Shapley split (`$contributions`: own part and interaction share), each term alone vs disease-free on the risk scale, and the local row equations (`$conflation$A`). Event impacts should keep the same results and tools as additive impacts.
   - Raw results (`raw_total()` in `deconflate.R`): the same model with the diseases unrelated (no association rows), for both models: `$totals$raw_total` (additive) and `$attributable$summary$raw_attributable` (event). Additive models with single-disease rows use a closed form (a weighted mean per disease); otherwise the model is refitted on the independent joint. The draws include it (`raw_total`); the diagnosis and the sensitivity tools skip it (`raw = FALSE`).
5. **Orchestration** (`deconflate.R`): `plan_model()` chooses between the full model and the **pairwise shortcut** (`pairwise.R`: 2x2 tables, eq. 14 of 2022). `shortcut_ok()` says when the shortcut is exact: additive impacts, a crude pair row for every pair, and single-disease impact rows. Any change to row semantics must keep the two paths giving the same answer. Even with the shortcut, the central estimate fits the joint once (up to 20 diseases) as a feasibility check. `finish_result()` adds `$estimates`, `$heterogeneity` and `$unknown_pairs`.
6. **Uncertainty** (`draws.R`, `distributions.R`): `n_draws` input sets go through `run_point()` with the pooling frozen at the central estimate. Rejected draws are counted, not replaced.
7. **Diagnosis** (`check_inputs.R`, `deconflate(diagnose = TRUE)`): flags diseases whose adjusted impact is on the other side of zero (of 1 for event impacts) from the raw impact or has shrunk by more than 90%, and reruns the model (`rerun_result()`) with one input changed at a time (pooling frozen, as in the draws): effects per SD of each uncertain input, threshold values of exact inputs, and rows with a standardised residual of 2 or more. Result: `$problems` (class `cm_diagnosis`; not to be confused with `cm_check_inputs()`, which checks the input files).
8. **Sensitivity** (`sensitivity.R`, `threshold.R`): `cm_sensitivity()` and `cm_threshold()` take a `deconflate()` result (its `$settings` hold the plan, risk and fit arguments) or a model. `rerun_result()` reruns a central result with inputs changed (pooling frozen) or on a changed model (replanned); the diagnosis uses it too. Input keys: `prob:`, `associations:`, `impacts:`, `risk`, and for inputs that are not entered `pair:d1:d2`, `interaction:d1;d2`, `or_ratio:d1:d2|d3` (built with the internal `set_association()`, `set_interaction()`, `set_three_way()`). Variance shares come from `$draws$uniforms`; `$ranks` from `draw_ranks()`.

The result classes are `cm_result` (additive impacts) and `cm_event_result` (event impacts). `example_supplement()` (three diseases) is the quickest model to experiment with. `tests/testthat/helper-models.R` has the shared fixtures `supp_model()`, `supp_hr_model()` and `read_five()`.
