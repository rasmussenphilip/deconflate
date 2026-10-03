# Roadmap

## v0.1 (released as 0.1.0)

Input constructors with metadata.

Association measures: OR, RR, RD, conditional probability, phi,
contingency table, independent and unknown.

Three adjustment methods: published, simultaneous (default) and global
(IPF/MaxEnt with pairwise interactions).

Diagnostics: reconstruction residuals, sign changes, condition numbers,
IPF convergence.

Feasibility check: a triple screen, plus an exact LP with `lpSolve`.

Productivity gaps on the proportion and absolute scales, monetary
valuation, and Shapley attribution (closed form, and cell by cell for
general losses).

Culling end to end: hazard ratio → excess risk → adjustment → hazard
ratio.

Reproduction of Rasmussen et al. (2022) Tables 8–10 and (2024) Table 5,
with documented differences.

Input distributions, samplers with outcome correlation, Monte Carlo with
economics and rejection reporting, and scenario reweighting.

Sensitivity: one-at-a-time, association screening, interaction screening
and scenario comparison.

Reporting and base-graphics plots.

Vignettes and pkgdown configuration.

Run `devtools::document()`, `test()` and `check()` on 0.0.0.9001 and fix
any issues.

GitHub Actions: `usethis::use_github_action("check-standard")` and
`usethis::use_pkgdown_github_pages()`.

Review the reproduction vignette against the original analysis files
(2024: 1st-revision code and inputs; 2022: original files unavailable,
corrigendum checked).

Tag v0.1.0.

## v0.2 (in progress)

Culling hazard ratios: a separate multiplicative adapter
([`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md):
snapshot, first-order and published) with stratified adjustment sets.

[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md):
culling attributable to disease without double counting, with Shapley
allocation.

Reading inputs from CSV files or data frames, with all problems reported
at once; template files; one impact file per analysis.

Method comparison tables for models, analyses, hazard-ratio models and
Monte Carlo runs.

Monte Carlo stability diagnostics with suggestions; Latin hypercube
(replicate blocks) and importance sampling (support checks).

External review (Bob): one impact vector per analysis in arbitrary
units; explicit estimands with the exact projection mapping for
additive-regression coefficients; feasibility screening by default;
condition classes; three-way scenarios; complete Shapley allocation;
reproduction functions; the reviewer’s regression tests.

Run `devtools::document()`, `test()` and `check()`, and fix any issues.

Tag v0.2.0.

## v0.3 (released as 0.3.0)

Second external review: importance-sampling support as intervals, one
finiteness check everywhere, `insufficient_info` precision status,
three-way documentation, converged joints for simulation,
single-analysis template, explicit snapshot hazard-ratio estimands,
historical conversions confined to reproduction.

Threshold searches
([`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)):
rank, sign, target and relative-change conclusions for one input at a
time.

Sampled backend for the joint distribution (Monte Carlo moment matching,
sampling, raking) beyond about 20 diseases.

Uncertainty given in the input tables (`dist`, `p1`-`p4`), with no
separate uncertainty files; three-way ratios can be uncertain.

Five-disease example that uses every feature
(`inst/extdata/five_diseases`, with `run_all_features.R`).

Run `devtools::document()`, `test()` and `check()` (0 errors, 0
warnings, 0 notes).

Tag v0.3.0.

## v0.4 (planned)

In this order:

**Test 0.3.0 in use.**

Run `run_all_features.R` end to end and review every printed result as a
new user would.

Have someone else work through the CSV workflow with their own data,
without help.

Optional: a short third review of the 0.3 additions (in-table
uncertainty,
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md),
the sampled backend).

**Probabilistic thresholds** (the one new feature): combine
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
with Monte Carlo runs, e.g. the input value at which a disease ranks
first (or the aggregate exceeds a target) in a given share of draws,
with Monte Carlo error.

**Simpler, more intuitive input tables.**

Accept missing optional columns everywhere with sensible defaults, and
fewer required columns.

Disease labels (full names) next to the short ids, used in printed
tables and plots.

Clearer column names and messages; review the template and the example
folders accordingly.

**Clearer output and summary tables.** Consistent column names and order
across
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)
and the Monte Carlo summaries; readable rounding and units; labels.

**Better plots.** A consistent look across plot methods, labels and
units on axes, and plots for uncertainty and thresholds.

**Dependencies.** Keep the package on base R (see the Dependencies
section of the README); list any new dependency there with its reason.

Afterwards: a Zenodo DOI for each release, and possibly a CRAN
submission and a short software paper.

## To discuss

- Pairwise-only bounds as a complementary sensitivity output for large
  systems; pseudo-likelihood fitting when individual-level disease
  records are available.
- Higher-order and non-linear interaction models in the adjustment
  itself.

Dropped (October 2026): advanced component accounting, and lost
productive life as a culling outcome.
