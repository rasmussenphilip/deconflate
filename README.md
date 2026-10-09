# deconflate

Comorbidity adjustment ("de-conflation") of disease impact estimates, so that the impacts of multiple associated diseases can be aggregated without double counting. Based on Rasmussen et al. (2022, *Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.* 107:6945–6970).

> **Status:** latest release [v0.3.0](https://github.com/rasmussenphilip/deconflate/releases); changes are listed in `NEWS.md` and plans in `ROADMAP.md`. The test reference values were computed independently in Python (`inst/validation/`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("rasmussenphilip/deconflate", build_vignettes = TRUE)
```

## Dependencies

The package itself needs only R (>= 4.1.0) and the base packages that come
with it: `graphics`, `grDevices`, `stats` and `utils`. Nothing else is
installed with it.

Optional packages (`Suggests`):

| Package | Used for |
|---|---|
| `lpSolve` | The exact feasibility check (`check_feasibility(method = "lp")`, `deconflate(feasibility = "lp")`). Without it, the triple screen is used. |
| `knitr`, `rmarkdown` | Building the vignettes. |
| `testthat` (>= 3.0.0) | Running the tests. |
| `remotes` | Installing from GitHub (see above). |

For development: `devtools`, `roxygen2` (documentation), `usethis`
(versions and releases) and `pkgdown` (the website). The reference values
in `inst/validation/` were computed independently in Python 3 with `numpy`
and `scipy`; they are not needed to use or test the package.

## How it works

A **population** holds the disease probabilities and their associations. Each **analysis** adds one vector of additive impact estimates, in any units (kg of milk, percent of yield, days, euros, welfare scores); results come back in the same units. Each estimate declares its estimand: a crude comparison (`"crude"`) or the coefficient of an additive regression adjusted for named diseases (`"adjusted_linear"`).

| Method | What it does |
|---|---|
| `"simultaneous"` (default) | Exact solution of the additive impact equations `raw = A b`, using pairwise 2×2 tables. Reports sign changes, reconstruction residuals and the conditioning of `A`, and screens the pairs for joint feasibility. |
| `"global"` | Maximum-entropy distribution of disease combinations fitted by iterative proportional fitting, then the full equations including specified pairwise interactions (and optional three-way scenarios). |

Use `"simultaneous"` when impacts are additive and every pair of diseases has an association estimate or a defensible independence assumption; use `"global"` otherwise (impact interactions, unknown pairs). The proportional approximation of Rasmussen et al. (2022, eq. 16) is not a method of `deconflate()`: it is kept as `"published"` in `compare_methods()`, `cm_monte_carlo()` and the `reproduce_*()` functions, for comparison and reproduction only.

The package reports adjusted impacts and contributions in the units of the impacts. Converting them into a productivity gap or a monetary value is a few lines of base R (see `vignette("deconflate")`).

Hazard ratios (e.g. of culling) combine multiplicatively and have their own model, the snapshot hazard-multiplier model (`deconflate_hr()`, with explicitly named estimands; the default and recommended method), with `attributable_risk()`.

The simultaneous method works for any number of diseases. The global method enumerates disease combinations up to about 20 diseases; beyond that, `fit_joint(backend = "sampled")` fits the same model by Monte Carlo.

## Example

```r
library(deconflate)

# Global dairy inputs, Rasmussen et al. (2024): yield and fertility analyses
# on one population
gd <- example_global_dairy()
res <- deconflate(gd)
res$yield$contributions
compare_methods(gd)   # includes the published approximation, for comparison

# At what odds ratio between d2 and d3 does the ranking of d1 and d2 change?
cm_threshold(example_supplement(), "assoc:d2:d3", c(0.1, 100), conclusion = "rank")

# Culling hazard ratios (snapshot model)
hr <- deconflate_hr(example_global_dairy_hr())
attributable_risk(hr, overall_risk = 0.2366)

# Your own data
dir <- file.path(tempdir(), "my-inputs")
cm_template(dir)
inp <- cm_read_inputs(dir = dir)
deconflate(inp$model)

# Uncertainty, with shared population draws across analyses
mc <- cm_monte_carlo(sampler_global_dairy(), 1000, seed = 1)
summary(mc)

# The published tables
reproduce_rasmussen_2022()
```

## Features

- **Inputs:** `cm_diseases()`, `cm_associations()` (OR, RR, RD, conditional probability, phi, contingency tables, independent, unknown), `cm_three_way()`, `cm_population()`, `cm_impacts()`, `cm_interactions()`, `cm_model()` and `cm_analyses()`.
- **Your own data:** `cm_read_inputs()` reads CSV files or data frames and reports every problem at once; uncertain values carry their distribution in their own row (`dist`, `p1`-`p4`), mixed freely with point values. `cm_template()` writes example files, and `system.file("extdata", "five_diseases", package = "deconflate")` holds a five-disease example that uses every feature, with a script (`run_all_features.R`) that runs them all.
- **Adjustment:** `deconflate()` and `compare_methods()` (which also shows the published approximation, for comparison); contributions of each disease (closed-form Shapley values) in `attribute_burden()` and `contribution_table()`.
- **Joint distribution and feasibility:** `fit_joint()` (exact, or sampled for many diseases), `combination_probs()` and `check_feasibility()`.
- **Attribution of non-additive losses:** `shapley_by_cell()`.
- **Hazard ratios:** `cm_hazard_ratios()`, `cm_hr_model()`, `deconflate_hr()` and `attributable_risk()`.
- **Uncertainty:** `dist_*()` distributions, `cm_sampler()` and `cm_batch_sampler()`, `cm_monte_carlo()` (several methods on the same draws, rejection reporting, Latin hypercube and importance sampling), stability checks with suggestions (`cm_diagnose()`, `cm_suggest_proposal()`), and reweighting (`cm_scenario()`, `cm_reweight()`).
- **Sensitivity and thresholds:** `cm_threshold()` (where a ranking, a sign or the aggregate changes as one input varies), `sensitivity_oat()`, `screen_associations()`, `screen_interactions()`, `screen_three_way()` and `compare_scenarios()`.
- **Reproduction:** `reproduce_rasmussen_2022()` and `reproduce_rasmussen_2024()`.
- **Plots** (base graphics) and seven vignettes. Start with `vignette("deconflate")`.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
