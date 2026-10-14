# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`deconflate` is an R package (base R only) that adjusts ("de-conflates") disease impact estimates for associations between diseases, so impacts can be summed without double counting (Rasmussen et al. 2022, 2024). Version 0.4 is being built item by item, following `ROADMAP.md`; `NEWS.md` records each item's breaking changes. Read the current ROADMAP item before changing behaviour.

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
- **The published method is for comparison only:** the proportional approximation of 2022 (eq. 16) and the HR − 1 conversion of 2024 belong in `compare_methods()` and `reproduce_rasmussen_2022()`/`_2024()` only. Calling it from `deconflate()`, `cm_threshold()` or the screens must raise `deconflate_unsupported`. Several files exist only for reproduction: `legacy.R`, `montecarlo.R`, `sampler.R` and `adjust_published()`. In `examples.R`, the reproduction-only helpers are `*_analyses()`, `uk_dairy_2022_economics()`, `uk_dairy_2022_culling_valuation()` and `global_dairy_sampler()`. Leave `global_dairy_inputs()`, `global_dairy_dists()` and `uk_dairy_2022_inputs()` alone, because the exported `example_*()` functions use them. `reproduce_*()` must keep matching the published tables.
- **Valuation stays out:** the package reports results in the impacts' own units. Productivity gaps and money values are left to the user, except inside `reproduce_*()`.
- **Classed errors:** raise errors with `cm_abort(msg, class = ...)` and warnings with `cm_warn()`, using a class from the list above `cm_abort()` in `R/utils.R` (documented for users in `?deconflate-package`, Conventions; keep the two in step). Draws, `compare_methods()` and the screens catch them by class, for example to count rejected draws. Retired arguments, columns and measures get an error or message that says what to write instead (see `retired_args()` and `old_args()`).
- **Finishing a roadmap item:** it updates the tests, vignettes, README, NEWS, the ROADMAP checkboxes, the five-disease example (`inst/extdata/five_diseases/`) and its `run_all_features.R`.
- **Style:** user-facing text and comments use British spelling, short plain sentences and the package's own vocabulary: row, combination, `among`, `adjusted for`, impact, interaction. The model assumes no species, outcome or unit. Printed tables and plots show only the user's own labels and units, although the examples and documentation use dairy cattle.

## Architecture

The package has one model. Every input row becomes a constraint on a maximum-entropy joint distribution of disease combinations. The impacts are the terms that reproduce every impact row. The pipeline:

1. **Inputs** (`inputs.R`, `import.R`): `cm_diseases()`, `cm_associations()`, `cm_impacts()`, `cm_population()` and `cm_model()`; `cm_read_inputs()` takes a CSV path or data frame per table and reports every problem at once (`cm_check_inputs()`). Association and impact rows share the core columns `disease`, `among`, `adjusted_for`, `measure` and `value`; associations also have `associated_with`. A `disease` can be a combination (`"MAS;SCK"`), and `among` a subgroup (`"SCK"`, `"no SCK"`, `"no others"`). Uncertainty goes in the same row (`dist`, `p1`-`p4`). Distribution keys are `prob:<disease>`, `associations:<row>` and `impacts:<row>`; the old `assoc:`, `impact:`, `inter:` and `three:` keys still work when they name one row.
2. **Row specs** (`rows.R`): each row is parsed into a spec: G = `disease`, H = `associated_with`, `among` as present/absent sets (`apos`/`aneg`), and the adjustment set. This file also pools several rows about one quantity (DerSimonian-Laird by default, or fixed effect) and turns association rows into event-probability targets for raking (`fit_events()`).
3. **Joint distribution** (`joint.R`): `fit_joint()` fits the maximum-entropy distribution by IPF/raking over all 2^n cells. Above 20 diseases, a sampled backend replaces it (`fit_joint_sampled()`). `event_joint()` in `hazard.R` is the wrapper the adjustment uses.
4. **Adjustment:**
   - Additive impacts (`adjust.R`): every impact row is linear in the disease impacts `b` and the interactions `delta`, giving `A theta = y`. One row per unknown is solved exactly; more rows go through weighted least squares (weights 1/SD², rows without a distribution exact). `identify_terms()` names the terms the rows do not determine.
   - Event impacts (`hazard.R`, `event_model = TRUE`): the snapshot hazard model with log-hazard terms, `h0` set from `overall_risk`, plus the attributable risk and its Shapley split.
5. **Orchestration** (`deconflate.R`): `plan_model()` chooses between the full model and the **pairwise shortcut** (`pairwise.R`: 2x2 tables, eq. 14 of 2022). `shortcut_ok()` says when the shortcut is exact: additive impacts, a crude pair row for every pair, and single-disease impact rows. Any change to row semantics must keep the two paths giving the same answer. Even with the shortcut, the central estimate fits the joint once (up to 20 diseases) as a feasibility check. `finish_result()` adds `$estimates`, `$heterogeneity` and `$unknown_pairs`.
6. **Uncertainty** (`draws.R`, `distributions.R`): `n_draws` input sets go through `run_point()` with the pooling frozen at the central estimate. Rejected draws are counted, not replaced.
7. **Downstream tools:** these build on `deconflate()`. The sensitivity tools (`cm_threshold()`, `sensitivity_oat()` and the screens) work on point estimates, and the ROADMAP plans a redesign of them.

The result classes are `cm_result` (additive impacts) and `cm_event_result` (event impacts). `example_supplement()` (three diseases) is the quickest model to experiment with. `tests/testthat/helper-models.R` has the shared fixtures `supp_model()`, `supp_hr_model()` and `read_five()`.
