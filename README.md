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
                    impacts = file.path(dir, "yield.csv"),
                    outcome = "milk yield loss", units = "% of yield")   # both optional
res <- deconflate(m)          # intervals from 1000 draws of the uncertain inputs
res$contributions

# Event impacts (culling): hazard ratios and a risk ratio
cull <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                       associations = file.path(dir, "associations.csv"),
                       impacts = file.path(dir, "culling.csv"), outcome = "culling")
deconflate(cull, event_model = TRUE, overall_risk = 0.25)

# Global dairy inputs, Rasmussen et al. (2024); unlisted pairs are unknown
gd <- example_global_dairy("yield")
res_gd <- deconflate(gd, n_draws = 0)
# Which inputs matter, including the associations nobody estimated?
cm_sensitivity(res_gd)

# At what odds ratio between d2 and d3 does the ranking of d1 and d2 change?
cm_threshold(example_supplement(), "pair:d2:d3", c(0.1, 100), conclusion = "rank")
# Which inputs could push the yield loss of LAM and SCK 10% below its current value, and where?
cm_threshold(res, conclusion = "total", diseases = c("LAM", "SCK"), target = 0.9, relative_to = "adjusted")

# Plots, and the data of each plot to draw it with other tools
plot(res, type = "contributions")
d <- plot_data(res, "contributions")
```

The plots use base graphics. With ggplot2 (not needed by the package), the
same data give your own version:

```r
library(ggplot2)
ggplot(d, aes(contribution, reorder(disease, contribution))) +
  geom_col(fill = "#2a78d6") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
  geom_point(aes(raw_contribution), shape = 21, fill = "white") +
  labs(x = "Contribution (% of yield)", y = NULL)
```

## Features

- **Inputs:** `cm_read_inputs()` reads CSV files or data frames (each table a named argument) and reports every problem at once (`cm_check_inputs()`); uncertain values carry their distribution in their own row (`dist`, `p1`-`p4`), mixed freely with point values. `cm_template()` writes example files, and `system.file("extdata", "five_diseases", package = "deconflate")` holds a five-disease example that uses every feature, with a script (`run_all_features.R`) that runs them all; the folder `pigs` holds an example outside dairy cattle. In R: `cm_diseases()`, `cm_associations()` (OR, RR, RD, phi, prevalence, OR ratio; about pairs, combinations and subgroups, crude or adjusted), `cm_population()`, `cm_impacts()` and `cm_model()`.
- **Adjustment:** `deconflate()` for additive and event impacts, with interactions, several rows per item (pooling with random or fixed effects, weighted least squares; `$heterogeneity`), intervals (`n_draws`), and a diagnosis (`diagnose = TRUE`) that shows which inputs drive adjusted impacts that change sign or nearly vanish, and which rows disagree with others; contributions of each disease (closed-form Shapley values), split into the disease's own part and its share of interactions, in `$contributions`, `contribution_table()` and `summary()`; the attributable risk of event impacts; the raw total (the result if the diseases were unrelated) beside the adjusted one.
- **Groups of diseases:** `deconflate(groups = list(ketosis = c("SCK", "CK")))` reports the forms of one disease (or any set) together, with intervals from the draws, while the model keeps them separate.
- **Joint distribution and feasibility:** `fit_joint()` (exact, or sampled for many diseases) and `check_feasibility()`.
- **Sensitivity and thresholds:** `cm_sensitivity()` varies every input one at a time, entered inputs (by their own SD, or by a relative change) and inputs you did not enter (unknown associations, interactions, three-way associations) in one table, with each input's share of the variance over the draws; `$ranks` gives how often each disease ranks first, second, ... over the draws; `cm_threshold()` finds where a ranking, a sign or the loss of any set of diseases (against a value, or a fraction of its adjusted or raw value) changes as an input varies, for named inputs or the ones that matter most; `cm_threshold_prob()` gives the probability of the conclusion over the draws of the other uncertain inputs, and where it crosses chosen levels. They also run without any association estimates, to show which associations would matter.
- **Validation:** `simulate_raw_impacts()` shows what studies would report given true impacts, so you can check that the model recovers them.
- **Plots** (base graphics): `plot()` of a result draws the raw and adjusted impacts, the contributions (with the raw ones), the ranks over the draws or the impact rows and their fit (`type`); of `cm_sensitivity()` the swings; of the threshold functions the compared quantity or the probability over each input's range. Titles and axes use the optional `outcome` and `units` given when the impacts are read (`cm_read_inputs()`) or built (`cm_impacts()`). `plot_data()` returns the data frame each plot draws, to draw it with ggplot2 or any other tool.
- **Vignettes:** seven. Start with `vignette("deconflate")`; `vignette("other-settings")` shows the package outside dairy cattle (PRRS and the respiratory diseases of growing pigs, with average daily gain, emissions and welfare; and chronic conditions and hospitalisation in people).

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

See `ROADMAP.md` for planned work and `NEWS.md` for changes.
