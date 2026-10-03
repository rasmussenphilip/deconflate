# deconflate

Comorbidity adjustment (“de-conflation”) of disease impact estimates, so
that the impacts of multiple associated diseases can be aggregated
without double counting. Based on Rasmussen et al. (2022, *Prev. Vet.
Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.*
107:6945–6970).

> **Status:** development version (0.1.0.9000), working towards v0.2.0.
> The test reference values were computed independently in Python
> (`inst/validation/python_reference.py`).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("rasmussenphilip/deconflate", build_vignettes = TRUE)
```

## Three adjustment methods

| Method | What it does |
|----|----|
| `"simultaneous"` (default) | Exact solution of the additive impact equations, `(I + Eᵀ) m = m_raw`, using pairwise 2×2 tables. Reports sign changes and reconstruction residuals. |
| `"published"` | Proportional approximation of Rasmussen et al. (2022), eqs. 15–16. Use it to reproduce the published results and to compare methods. |
| `"global"` | Maximum-entropy distribution of disease combinations fitted by iterative proportional fitting (IPF), then the full equations including specified pairwise interactions. |

“Exact” refers to the specified equations and assumptions, solved within
numerical tolerance.

## Example

``` r

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

- **Inputs** with metadata (probability type and time horizon, effect
  scale, reference population, prior adjustment). Association measures:
  OR, RR, RD, conditional probability, phi, contingency tables,
  independent and unknown.
- **Your own data:**
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  reads CSV files or data frames and reports every problem at once;
  [`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
  writes example files.
- **Adjustment:**
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
  and
  [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
  for side-by-side tables of the methods (also for Monte Carlo runs).
  Interactions between diseases are supported with the global method.
- **Joint distribution and feasibility:**
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md),
  [`combination_probs()`](https://rasmussenphilip.github.io/deconflate/reference/combination_probs.md)
  and
  [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md).
- **Gaps, losses and attribution:**
  [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md),
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md),
  [`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md),
  [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
  and
  [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md).
- **Culling:** hazard ratios adjusted under a multiplicative model
  (exact or first-order),
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  for the culling attributable to disease, and the published approaches
  of 2022 and 2024
  ([`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md),
  [`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md),
  [`adjusted_hr()`](https://rasmussenphilip.github.io/deconflate/reference/adjusted_hr.md)).
- **Uncertainty:** `dist_*()` distributions,
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
  (with optional correlation across outcomes),
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  (several methods on the same draws, rejection reporting, Latin
  hypercube and importance sampling), stability checks with suggestions
  ([`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md))
  and
  [`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
  (importance reweighting).
- **Sensitivity:**
  [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md),
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
  [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
  and
  [`compare_scenarios()`](https://rasmussenphilip.github.io/deconflate/reference/compare_scenarios.md).
- **Plots** (base graphics) and six vignettes. Start with
  [`vignette("deconflate")`](https://rasmussenphilip.github.io/deconflate/articles/deconflate.md).

## Development

``` r

devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
