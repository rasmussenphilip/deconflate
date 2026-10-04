# Roadmap

## v0.1 (released as 0.1.0)

- [x] Input constructors with metadata.
- [x] Association measures: OR, RR, RD, conditional probability, phi, contingency table, independent and unknown.
- [x] Three adjustment methods: published, simultaneous (default) and global (IPF/MaxEnt with pairwise interactions).
- [x] Diagnostics: reconstruction residuals, sign changes, condition numbers, IPF convergence.
- [x] Feasibility check: a triple screen, plus an exact LP with `lpSolve`.
- [x] Productivity gaps on the proportion and absolute scales, monetary valuation, and Shapley attribution (closed form, and cell by cell for general losses).
- [x] Culling end to end: hazard ratio → excess risk → adjustment → hazard ratio.
- [x] Reproduction of Rasmussen et al. (2022) Tables 8–10 and (2024) Table 5, with documented differences.
- [x] Input distributions, samplers with outcome correlation, Monte Carlo with economics and rejection reporting, and scenario reweighting.
- [x] Sensitivity: one-at-a-time, association screening, interaction screening and scenario comparison.
- [x] Reporting and base-graphics plots.
- [x] Vignettes and pkgdown configuration.
- [x] Run `devtools::document()`, `test()` and `check()` on 0.0.0.9001 and fix any issues.
- [x] GitHub Actions: `usethis::use_github_action("check-standard")` and `usethis::use_pkgdown_github_pages()`.
- [x] Review the reproduction vignette against the original analysis files (2024: 1st-revision code and inputs; 2022: original files unavailable, corrigendum checked).
- [x] Tag v0.1.0.

## v0.2 (released as 0.2.0)

- [x] Culling hazard ratios: a separate multiplicative adapter (`deconflate_hr()`: snapshot, first-order and published) with stratified adjustment sets.
- [x] `attributable_risk()`: culling attributable to disease without double counting, with Shapley allocation.
- [x] Reading inputs from CSV files or data frames, with all problems reported at once; template files; one impact file per analysis.
- [x] Method comparison tables for models, analyses, hazard-ratio models and Monte Carlo runs.
- [x] Monte Carlo stability diagnostics with suggestions; Latin hypercube (replicate blocks) and importance sampling (support checks).
- [x] External review (Bob): one impact vector per analysis in arbitrary units; explicit estimands with the exact projection mapping for additive-regression coefficients; feasibility screening by default; condition classes; three-way scenarios; complete Shapley allocation; reproduction functions; the reviewer's regression tests.
- [x] Run `devtools::document()`, `test()` and `check()`, and fix any issues.
- [x] Tag v0.2.0.

## v0.3 (released as 0.3.0)

- [x] Second external review: importance-sampling support as intervals, one finiteness check everywhere, `insufficient_info` precision status, three-way documentation, converged joints for simulation, single-analysis template, explicit snapshot hazard-ratio estimands, historical conversions confined to reproduction.
- [x] Threshold searches (`cm_threshold()`): rank, sign, target and relative-change conclusions for one input at a time.
- [x] Sampled backend for the joint distribution (Monte Carlo moment matching, sampling, raking) beyond about 20 diseases.
- [x] Uncertainty given in the input tables (`dist`, `p1`-`p4`), with no separate uncertainty files; three-way ratios can be uncertain.
- [x] Five-disease example that uses every feature (`inst/extdata/five_diseases`, with `run_all_features.R`).
- [x] Run `devtools::document()`, `test()` and `check()` (0 errors, 0 warnings, 0 notes).
- [x] Tag v0.3.0.

## v0.4 (planned)

In this order:

1. **Test 0.3.0 in use.**
   - [ ] Run `run_all_features.R` end to end and review every printed result as a new user would.
   - [ ] Have someone else work through the CSV workflow with their own data, without help.
   - [ ] Optional: a short third review of the 0.3 additions (in-table uncertainty, `cm_threshold()`, the sampled backend).
2. **Probabilistic thresholds** (the one new feature): combine `cm_threshold()` with Monte Carlo runs, e.g. the input value at which a disease ranks first (or the aggregate exceeds a target) in a given share of draws, with Monte Carlo error.
3. **Simpler, more intuitive input tables.**
   - [ ] Accept missing optional columns everywhere with sensible defaults, and fewer required columns.
   - [ ] Disease labels (full names) next to the short ids, used in printed tables and plots.
   - [ ] Clearer column names and messages; review the template and the example folders accordingly.
4. **Clearer output and summary tables.** Consistent column names and order across `deconflate()`, `compare_methods()`, `summary()`, `contribution_table()` and the Monte Carlo summaries; readable rounding and units; labels.
5. **Better plots.** A consistent look across plot methods, labels and units on axes, and plots for uncertainty and thresholds.
6. **Dependencies.** Keep the package on base R (see the Dependencies section of the README); list any new dependency there with its reason.

Afterwards: a Zenodo DOI for each release, and possibly a CRAN submission and a short software paper.

## To discuss

- Pairwise-only bounds as a complementary sensitivity output for large systems; pseudo-likelihood fitting when individual-level disease records are available.
- Higher-order and non-linear interaction models in the adjustment itself.

Dropped (October 2026): advanced component accounting, and lost productive life as a culling outcome.
