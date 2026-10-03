# Roadmap

## v0.1 (in progress: version 0.0.0.9001)

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
- [ ] Run `devtools::document()`, `test()` and `check()` on 0.0.0.9001 and fix any issues.
- [ ] GitHub Actions: `usethis::use_github_action("check-standard")` and `usethis::use_pkgdown_github_pages()`.
- [ ] Review the reproduction vignette against the original spreadsheets (2022 yield inputs; 2024 handling of negative draws).
- [ ] Tag v0.1.0.

## After v0.1

- Advanced component accounting.
- Threshold searches: input values at which conclusions or rankings change.
- Higher-order and non-linear interaction models in the adjustment itself (a per-cell loss function and a non-linear solve).
- Scaling beyond about 20 diseases, e.g. a pseudo-likelihood Ising fit or sampling instead of the full 2^n table.
