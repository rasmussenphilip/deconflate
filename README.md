# deconflate

Comorbidity adjustment ("de-conflation") of disease impact estimates, so that impacts of multiple associated diseases can be aggregated without double counting. Based on Rasmussen et al. (2022, *Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.* 107:6945–6970).

> **Status:** development skeleton (0.0.0.9000). The core methods have been implemented, and the test reference values have been cross-checked against an independent Python implementation (`inst/validation/python_reference.py`). The R code has not yet been run: see *First run* below.

## Installation

```r
# install.packages("remotes")
remotes::install_github("<user>/deconflate")
```

## First run (development)

```r
# install.packages(c("devtools", "roxygen2", "testthat"))
devtools::document()   # generates man/ and refreshes NAMESPACE
devtools::test()       # runs the test suite
devtools::check()      # full R CMD check
```

## Three adjustment methods

| Method | What it does |
|---|---|
| `"simultaneous"` (default) | Exact solution of the additive impact equations, `(I + Eᵀ) m = m_raw`, using pairwise 2×2 tables. Reports sign changes and reconstruction residuals. |
| `"published"` | Proportional approximation of Rasmussen et al. (2022), eqs. 15–16. Use it to reproduce the published results and to compare methods. |
| `"global"` | Maximum-entropy distribution of disease combinations fitted by iterative proportional fitting (IPF), then the full equations including specified pairwise interactions. |

"Exact" refers to the specified equations and assumptions, solved within numerical tolerance.

## Example

```r
library(deconflate)

m <- example_supplement()          # Supplementary File, Rasmussen et al. (2022)
compare_methods(m)

res <- deconflate(m)               # simultaneous (default)
productivity_gap(res, c(yield = 10000))

# Validation: simulate raw impacts from known truth and recover it
m$impacts <- simulate_raw_impacts(m, c(d1 = 0.02, d2 = 0.04, d3 = 0.06), outcome = "yield")
deconflate(m)$adjusted

# Global inputs from Rasmussen et al. (2024): the exact solution flags
# incompatible inputs that the published approximation hides
deconflate(example_global_dairy())$diagnostics
```

## Main functions

- **Inputs:** `cm_diseases()`, `cm_associations()` (OR, RR, RD, conditional probability, phi, contingency table, independent, unknown), `cm_impacts()`, `cm_interactions()` and `cm_model()`. Metadata covers the probability type, time horizon, effect scale, reference population and prior adjustment.
- **Pairwise algebra:** `or_to_joint()`, `association_to_joint()`, `pair_tables()` and `excess_matrix()`.
- **Joint distribution:** `fit_joint()` and `combination_probs()`.
- **Adjustment:** `deconflate()` and `compare_methods()`.
- **Aggregation:** `attribute_burden()` (closed-form Shapley) and `productivity_gap()`.
- **Conversions:** `hr_to_risk()` (proportional hazards, or the published OR approximation) and `excess_to_hr()`.
- **Uncertainty:** `cm_monte_carlo()` (shared draws across outcomes, with rejection reporting) and `cm_reweight()` (importance reweighting for scenarios).
- **Validation:** `simulate_raw_impacts()`.

See `ROADMAP.md` for the v0.1 features still to be built.
