# deconflate

Comorbidity adjustment ("de-conflation") of disease impact estimates, so that the impacts of multiple associated diseases can be aggregated without double counting. Based on Rasmussen et al. (2022, *Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy Sci.* 107:6945–6970).

> **Status:** latest release [v0.5.0](https://github.com/rasmussenphilip/deconflate/releases); changes are listed in `NEWS.md` and plans in `ROADMAP.md`. The test reference values were computed independently in Python (`inst/validation/`).

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
| `lpSolve` | The exact feasibility check of the pairwise tables (`check_feasibility(method = "lp")`). Without it, the triple screen is used. |
| `knitr`, `rmarkdown` | Building the vignettes. |
| `testthat` (>= 3.0.0) | Running the tests. |
| `remotes` | Installing from GitHub (see above). |

For development: `devtools`, `roxygen2` (documentation), `usethis`
(versions and releases) and `pkgdown` (the website). The reference values
in `inst/validation/` were computed independently in Python 3 with `numpy`
and `scipy`; they are not needed to use or test the package.

## How it works

A model has three tables: the diseases, the association estimates and **one** table of impact estimates. You give each table as a CSV file (with any name) or a data frame, and adjust one impact table per call; for several outcomes, repeat the call. The package assumes no species, outcome or unit.

Every association and impact row has the same core columns: `disease`, `among` (a subgroup), `adjusted_for`, `measure` and `value`. A row can be about one disease or a combination (`"MAS;SCK"`), hold in the whole population or within a subgroup (`"SCK"`, `"no SCK"`, `"no others"`), and be crude or adjusted for other diseases; several rows can describe the same thing (several studies), and are pooled or fitted together.

* **Additive impacts** (the default) are in any units (kg of milk, percent of yield, days, euros, welfare scores), the same for every row; results come back in those units. A row about a combination adds an interaction.
* **Event impacts** (`event_model = TRUE`, with `overall_risk` or the risk of disease-free animals, `baseline_risk`) are hazard ratios, rate ratios, risk ratios, odds ratios or risk differences of an event such as culling or death, which can be mixed. They are adjusted with the snapshot hazard model, and the result includes the risk attributable to disease, split by disease, and each adjusted result on the risk scale. Event impacts have the same tools as additive impacts.

There is one model: the maximum-entropy distribution of disease combinations that matches every association row, and the impacts (and interactions) that reproduce every impact row. Pairs without an association are unknown and are filled in by the fit; an odds ratio of 1 states that two diseases are unrelated. When the pairwise 2×2 tables give the same answer (additive impacts, crude pair rows for every pair and nothing else, impact rows about single diseases in the whole population), `deconflate()` uses them directly.

Inputs can carry distributions; `deconflate()` then runs 1000 draws by default (`n_draws`) and reports intervals. The proportional approximation used in Rasmussen et al. (2022, eq. 16; 2024) is not part of the package.

The package reports adjusted impacts and contributions in the units of the impacts. Converting them into a productivity gap or a monetary value is a few lines of base R (see `vignette("deconflate")`).

## Example

```r
library(deconflate)

# Your own data: one call per impact table, file names of your choice
dir <- file.path(tempdir(), "my-inputs")
cm_template(dir)
m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                    associations = file.path(dir, "associations.csv"),
                    impacts = file.path(dir, "yield.csv"))
res <- deconflate(m)          # intervals from 1000 draws of the uncertain inputs
res$contributions

# Event impacts (culling): hazard ratios and a risk ratio
cull <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                       associations = file.path(dir, "associations.csv"),
                       impacts = file.path(dir, "culling.csv"))
deconflate(cull, event_model = TRUE, overall_risk = 0.25)

# Global dairy inputs, Rasmussen et al. (2024); unlisted pairs are unknown
gd <- example_global_dairy("yield")
deconflate(gd, n_draws = 0)

# At what odds ratio between d2 and d3 does the ranking of d1 and d2 change?
cm_threshold(example_supplement(), "assoc:d2:d3", c(0.1, 100), conclusion = "rank")
```

## Features

- **Inputs:** `cm_read_inputs()` reads CSV files or data frames (each table a named argument) and reports every problem at once (`cm_check_inputs()`); uncertain values carry their distribution in their own row (`dist`, `p1`-`p4`), mixed freely with point values. `cm_template()` writes example files, and `system.file("extdata", "five_diseases", package = "deconflate")` holds a five-disease example that uses every feature, with a script (`run_all_features.R`) that runs them all. In R: `cm_diseases()`, `cm_associations()` (OR, RR, RD, phi, prevalence, OR ratio; about pairs, combinations and subgroups, crude or adjusted), `cm_population()`, `cm_impacts()` and `cm_model()`.
- **Adjustment:** `deconflate()` for additive and event impacts, with interactions, several rows per item (pooling with random or fixed effects, weighted least squares; `$heterogeneity`), intervals (`n_draws`), and a diagnosis (`diagnose = TRUE`) that shows which inputs drive adjusted impacts that change sign or nearly vanish, and which rows disagree with others; contributions of each disease (closed-form Shapley values), split into the disease's own part and its share of interactions, in `$contributions`, `contribution_table()` and `summary()`; the attributable risk of event impacts; the raw total (the result if the diseases were unrelated) beside the adjusted one.
- **Joint distribution and feasibility:** `fit_joint()` (exact, or sampled for many diseases) and `check_feasibility()`.
- **Sensitivity and thresholds:** `cm_threshold()` (where a ranking, a sign or the aggregate changes as one input varies), `sensitivity_oat()`, `screen_associations()`, `screen_interactions()` and `screen_three_way()`; they also run without any association estimates, to show which associations would matter.
- **Validation:** `simulate_raw_impacts()` shows what studies would report given true impacts, so you can check that the model recovers them.
- **Plots** (base graphics) and six vignettes. Start with `vignette("deconflate")`.

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
