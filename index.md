# deconflate

Comorbidity adjustment (“de-conflation”) of disease impact estimates, so
that the impacts of multiple associated diseases can be aggregated
without double counting. Based on Rasmussen et al. (2022, *Prev. Vet.
Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.*
107:6945–6970).

> **Status:** development version (0.1.0.9000), working towards v0.2.0,
> which changes the interface (see `NEWS.md`). The test reference values
> were computed independently in Python (`inst/validation/`).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("rasmussenphilip/deconflate", build_vignettes = TRUE)
```

## How it works

A **population** holds the disease probabilities and their associations.
Each **analysis** adds one vector of additive impact estimates, in any
units (kg of milk, percent of yield, days, euros, welfare scores);
results come back in the same units. Each estimate declares its
estimand: a crude comparison (`"crude"`) or the coefficient of an
additive regression adjusted for named diseases (`"adjusted_linear"`).

| Method | What it does |
|----|----|
| `"simultaneous"` (default) | Exact solution of the additive impact equations `raw = A b`, using pairwise 2×2 tables. Reports sign changes, reconstruction residuals and the conditioning of `A`, and screens the pairs for joint feasibility. |
| `"published"` | Proportional approximation of Rasmussen et al. (2022), eq. 16, for crude estimates. Use it to reproduce the published results and to compare methods. |
| `"global"` | Maximum-entropy distribution of disease combinations fitted by iterative proportional fitting, then the full equations including specified pairwise interactions (and optional three-way scenarios). |

Hazard ratios (e.g. of culling) combine multiplicatively and have their
own adapter,
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md),
with
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).

## Example

``` r

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

- **Inputs:**
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md),
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
  (OR, RR, RD, conditional probability, phi, contingency tables,
  independent, unknown),
  [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md),
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md),
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md),
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md),
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  and
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md).
- **Your own data:**
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  reads CSV files or data frames and reports every problem at once;
  [`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
  writes example files.
- **Adjustment:**
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  and
  [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md);
  contributions of each disease (closed-form Shapley values) in
  [`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md)
  and
  [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md).
- **Joint distribution and feasibility:**
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md),
  [`combination_probs()`](https://rasmussenphilip.github.io/deconflate/reference/combination_probs.md)
  and
  [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md).
- **Optional helpers:**
  [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
  and
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)
  turn an aggregate into a gap and a value;
  [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
  attributes non-additive losses.
- **Hazard ratios:**
  [`cm_hazard_ratios()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hazard_ratios.md),
  [`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md),
  [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
  and
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).
- **Uncertainty:** `dist_*()` distributions,
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
  and
  [`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md),
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  (several methods on the same draws, rejection reporting, Latin
  hypercube and importance sampling), stability checks with suggestions
  ([`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md),
  [`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md)),
  reweighting
  ([`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md),
  [`cm_reweight()`](https://rasmussenphilip.github.io/deconflate/reference/cm_reweight.md))
  and
  [`cm_mc_gap()`](https://rasmussenphilip.github.io/deconflate/reference/cm_mc_gap.md).
- **Sensitivity:**
  [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md),
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
  [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md),
  [`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md)
  and
  [`compare_scenarios()`](https://rasmussenphilip.github.io/deconflate/reference/compare_scenarios.md).
- **Reproduction:**
  [`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
  and
  [`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).
- **Plots** (base graphics) and six vignettes. Start with
  [`vignette("deconflate")`](https://rasmussenphilip.github.io/deconflate/articles/deconflate.md).

## Development

``` r

devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
