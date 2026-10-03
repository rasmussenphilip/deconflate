# deconflate

Comorbidity adjustment ("de-conflation") of disease impact estimates, so that the impacts of multiple associated diseases can be aggregated without double counting. Based on Rasmussen et al. (2022, *Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.* 107:6945–6970).

> **Status:** development version (0.0.0.9001), working towards v0.1.0. The test reference values were computed independently in Python (`inst/validation/python_reference.py`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rasmussenphilip/deconflate", build_vignettes = TRUE)
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

# UK dairy example, Rasmussen et al. (2022)
m <- example_uk_dairy_2022()
eco <- uk_dairy_2022_economics()
res <- deconflate(m, method = "published")
summary(res, economics = eco)
plot_burden(res, eco)

# Exact solution, with diagnostics for inconsistent inputs
compare_methods(m)

# Uncertainty: Monte Carlo of the global inputs, Rasmussen et al. (2024).
# Means of the adjusted impacts reproduce Table 5 (culling: add 1 to get HRs).
mc <- cm_monte_carlo(sampler_global_dairy(), 1000, method = "published", seed = 1)
summary(mc)
```

## Features

- **Inputs** with metadata (probability type and time horizon, effect scale, reference population, prior adjustment). Association measures: OR, RR, RD, conditional probability, phi, contingency tables, independent and unknown.
- **Adjustment:** `deconflate()` and `compare_methods()`. Interactions between diseases are supported with the global method.
- **Joint distribution and feasibility:** `fit_joint()`, `combination_probs()` and `check_feasibility()`.
- **Gaps, losses and attribution:** `productivity_gap()`, `value_losses()`, `attribute_burden()`, `shapley_by_cell()` and `contribution_table()`.
- **Culling:** `hr_conversion()`, `as_impacts()` and `adjusted_hr()` (proportional hazards, or the published approaches of 2022 and 2024).
- **Uncertainty:** `dist_*()` distributions, `cm_sampler()` (with optional correlation across outcomes), `cm_monte_carlo()` (with rejection reporting) and `cm_scenario()` (importance reweighting).
- **Sensitivity:** `sensitivity_oat()`, `screen_associations()`, `screen_interactions()` and `compare_scenarios()`.
- **Plots** (base graphics) and four vignettes. Start with `vignette("deconflate")`.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
