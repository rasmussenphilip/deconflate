# Roadmap to v0.1

This list follows the agreed v0.1 scope (October 2026).

## Done in the skeleton

- [x] Input constructors with metadata: probability type and incidence conversion, time horizon, reference population, effect scale, units, direction, and prior adjustment of both impacts and associations.
- [x] Association measures: OR, RR, RD, conditional probability, phi, contingency tables, independent and unknown. Unknown is kept distinct from independence.
- [x] Three adjustment methods: published, simultaneous (default) and global (IPF/MaxEnt with pairwise interactions).
- [x] Diagnostics: reconstruction residuals, sign changes, condition numbers and IPF convergence.
- [x] Productivity gaps (decrease and increase outcomes) and closed-form Shapley attribution for pairwise interactions.
- [x] Hazard-ratio conversion: proportional hazards, plus the published OR approximation.
- [x] Monte Carlo with shared draws across outcomes and counting of rejected draws, plus importance reweighting for scenarios with an ESS.
- [x] Simulated-population validation helper and tests against Python reference values.

## To do for v0.1

- [ ] Run `devtools::document()`, `devtools::test()` and `devtools::check()`, then fix any issues.
- [ ] Linear-programming feasibility check of the pairwise tables before IPF (e.g. `lpSolve` in Suggests), with a diagnosis of the most conflicting pairs.
- [ ] Shapley allocation for higher-order additive terms (equal split among the diseases involved) and, for non-additive effects, per-cell computation over the joint table.
- [ ] Culling outcomes end to end: HR → excess risk → adjustment → HR (reproducing 2022 Table 8 with `or_approx`, then the corrected conversion).
- [ ] Full reproduction of 2022 Tables 8–10 and 2024 Table 5 (Monte Carlo means), with a vignette documenting where and why results differ.
- [ ] Sampler helpers for common input distributions (beta, PERT, normal and lognormal on the OR scale), including defensive mixtures for scenario reweighting.
- [ ] Optional correlations between outcome impacts within a disease.
- [ ] Sensitivity tools: pairwise screening of uncertain or missing associations and interaction scenarios, comparing totals and rankings.
- [ ] Reporting: contribution tables, uncertainty intervals and optional plots (in Suggests only).
- [ ] Vignettes: getting started, method comparison, interactions and assumptions, uncertainty and scenarios.
- [ ] GitHub Actions R-CMD-check workflow and a pkgdown site.

## After v0.1

- Advanced component accounting, threshold searches, and higher-order or non-linear interaction models.
- Scaling beyond about 20 diseases, e.g. a pseudo-likelihood Ising fit or sampling instead of the full 2^n table.
