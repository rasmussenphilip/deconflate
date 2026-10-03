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

## v0.3 (in progress)

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

Run `devtools::document()`, `test()` and `check()`, and fix any issues.

Tag v0.3.0.

## To discuss

- Probabilistic thresholds (e.g. the input value at which a disease
  ranks first with 95% probability).
- Advanced component accounting.
- Lost productive life as a culling outcome.
- Pairwise-only bounds as a complementary sensitivity output for large
  systems; pseudo-likelihood fitting when individual-level disease
  records are available.
- Higher-order and non-linear interaction models in the adjustment
  itself.
