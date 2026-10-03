# deconflate

Comorbidity adjustment ("de-conflation") of disease impact estimates, so that the impacts of multiple associated diseases can be aggregated without double counting. Based on Rasmussen et al. (2022, *Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.* 107:6945–6970).

> **Status:** development version (0.1.0.9000), working towards v0.2.0, which changes the interface (see `NEWS.md`). The test reference values were computed independently in Python (`inst/validation/`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rasmussenphilip/deconflate", build_vignettes = TRUE)
```

## How it works

A **population** holds the disease probabilities and their associations. Each **analysis** adds one vector of additive impact estimates, in any units (kg of milk, percent of yield, days, euros, welfare scores); results come back in the same units. Each estimate declares its estimand: a crude comparison (`"crude"`) or the coefficient of an additive regression adjusted for named diseases (`"adjusted_linear"`).

| Method | What it does |
|---|---|
| `"simultaneous"` (default) | Exact solution of the additive impact equations `raw = A b`, using pairwise 2×2 tables. Reports sign changes, reconstruction residuals and the conditioning of `A`, and screens the pairs for joint feasibility. |
| `"published"` | Proportional approximation of Rasmussen et al. (2022), eq. 16, for crude estimates. Use it to reproduce the published results and to compare methods. |
| `"global"` | Maximum-entropy distribution of disease combinations fitted by iterative proportional fitting, then the full equations including specified pairwise interactions (and optional three-way scenarios). |

Hazard ratios (e.g. of culling) combine multiplicatively and have their own adapter, `deconflate_hr()`, with `attributable_risk()`.

## Example

```r
library(deconflate)

# Global dairy inputs, Rasmussen et al. (2024): yield and fertility analyses
# on one population
gd <- example_global_dairy()
res <- deconflate(gd)
res$yield$contributions
compare_methods(gd)

# Culling hazard ratios
hr <- deconflate_hr(example_global_dairy_hr())
attributable_risk(hr, overall_risk = 0.2366)

# Your own data: one impact file per analysis
dir <- file.path(tempdir(), "my-inputs")
cm_template(dir)
inp <- cm_read_inputs(dir = dir)
deconflate(inp$analyses)

# Uncertainty, with shared population draws across analyses
mc <- cm_monte_carlo(sampler_global_dairy(), 1000, method = "published", seed = 1)
summary(mc)

# The published tables
reproduce_rasmussen_2022()
```

## Features

- **Inputs:** `cm_diseases()`, `cm_associations()` (OR, RR, RD, conditional probability, phi, contingency tables, independent, unknown), `cm_three_way()`, `cm_population()`, `cm_impacts()`, `cm_interactions()`, `cm_model()` and `cm_analyses()`.
- **Your own data:** `cm_read_inputs()` reads CSV files or data frames and reports every problem at once; `cm_template()` writes example files.
- **Adjustment:** `deconflate()` and `compare_methods()`; contributions of each disease (closed-form Shapley values) in `attribute_burden()` and `contribution_table()`.
- **Joint distribution and feasibility:** `fit_joint()`, `combination_probs()` and `check_feasibility()`.
- **Optional helpers:** `productivity_gap()` and `value_losses()` turn an aggregate into a gap and a value; `shapley_by_cell()` attributes non-additive losses.
- **Hazard ratios:** `cm_hazard_ratios()`, `cm_hr_model()`, `deconflate_hr()` and `attributable_risk()`.
- **Uncertainty:** `dist_*()` distributions, `cm_sampler()` and `cm_batch_sampler()`, `cm_monte_carlo()` (several methods on the same draws, rejection reporting, Latin hypercube and importance sampling), stability checks with suggestions (`cm_diagnose()`, `cm_suggest_proposal()`), reweighting (`cm_scenario()`, `cm_reweight()`) and `cm_mc_gap()`.
- **Sensitivity:** `sensitivity_oat()`, `screen_associations()`, `screen_interactions()`, `screen_three_way()` and `compare_scenarios()`.
- **Reproduction:** `reproduce_rasmussen_2022()` and `reproduce_rasmussen_2024()`.
- **Plots** (base graphics) and six vignettes. Start with `vignette("deconflate")`.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
