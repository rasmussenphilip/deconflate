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

## v0.2 (in progress)

- [x] Culling hazard ratios: a separate multiplicative adapter (`deconflate_hr()`: snapshot, first-order and published) with stratified adjustment sets.
- [x] `attributable_risk()`: culling attributable to disease without double counting, with Shapley allocation.
- [x] Reading inputs from CSV files or data frames, with all problems reported at once; template files; one impact file per analysis.
- [x] Method comparison tables for models, analyses, hazard-ratio models and Monte Carlo runs.
- [x] Monte Carlo stability diagnostics with suggestions; Latin hypercube (replicate blocks) and importance sampling (support checks).
- [x] External review (Bob): one impact vector per analysis in arbitrary units; explicit estimands with the exact projection mapping for additive-regression coefficients; feasibility screening by default; condition classes; three-way scenarios; complete Shapley allocation; reproduction functions; the reviewer's regression tests.
- [x] Run `devtools::document()`, `test()` and `check()`, and fix any issues.
- [ ] Tag v0.2.0.

## To discuss

- Threshold searches: input values at which conclusions or rankings change.
- Advanced component accounting.
- Lost productive life as a culling outcome.
- Scaling beyond about 20 diseases (the global method and the exact culling model enumerate all 2^n combinations), e.g. a pseudo-likelihood Ising fit or sampling.
- Higher-order and non-linear interaction models in the adjustment itself.
