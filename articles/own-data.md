# Using your own data

``` r

library(deconflate)
```

Inputs are tables: the diseases, their associations, one impact table
and, if needed, interactions and three-way terms. Each table can be a
CSV file (with any name) or a data frame typed in R.
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
checks every table, reports all problems at once, and returns a model
([`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md))
that
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
adjusts.

## The tables

Each table is a named argument of
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):

| Argument | Contents | Required columns |
|----|----|----|
| `diseases` | one row per disease | `id`, `value` |
| `associations` | one row per pair with an association estimate | `disease1`, `disease2`, `value` |
| `impacts` | the raw impact estimates of one outcome, one row per disease | `disease`, `value` |
| `interactions` | optional pairwise impact interactions (additive impacts only) | `disease1`, `disease2`, `value` |
| `three_way` | optional three-way association terms | `disease1`, `disease2`, `disease3`, `ratio` |

`diseases` and `impacts` are required.
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
also needs at least one association; without any, the sensitivity tools
([`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md))
show how much associations could matter. The optional columns are
described in
[`?cm_read_inputs`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):

- `diseases`: `type` (`prevalence` (default), `probability` or
  `incidence_rate`), `time_horizon`, `reference_population`, `source`;
- `associations`: `measure` (`OR` (default), `RR`, `RD`, `cond_prob` or
  `phi`), `adjusted`, `adjusted_for`, `source`;
- `impacts`: `estimand` (`crude` (default) or `adjusted_linear`),
  `adjusted_for`, `source`, and one `label` and one `units` per table; a
  `measure` column makes them event impacts (see below);
- `interactions` and `three_way`: `source`.

Every table can also have a free-text `note` column, and the columns
`dist` and `p1`-`p4` for uncertain values (see Uncertainty below).
Column names are not case-sensitive.

Every disease needs a row in the impact table (use 0 for no impact).
Additive impacts can be in any units, the same for every row; results
come back in those units.

## Start from the template

[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set of files to a folder:

``` r

dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 6 files to /tmp/RtmpATEobY/my-inputs
list.files(dir)
#> [1] "associations.csv"       "culling.csv"            "diseases.csv"          
#> [4] "three_way.csv"          "yield_interactions.csv" "yield.csv"
read.csv(file.path(dir, "yield.csv"))
#>   disease value estimand adjusted_for           label      units
#> 1     LAM  4.81    crude           NA milk yield loss % of yield
#> 2     SCK  8.40    crude           NA milk yield loss % of yield
#> 3     MET  5.61    crude           NA milk yield loss % of yield
#>                source   dist   p1   p2 p3 p4             note
#> 1 Illustrative values normal 4.81 0.87 NA NA normal(mean, sd)
#> 2 Illustrative values normal 8.40 1.19 NA NA normal(mean, sd)
#> 3 Illustrative values          NA   NA NA NA
```

The files are three diseases (`diseases.csv`), their odds ratios
(`associations.csv`), additive impacts on milk yield in percent
(`yield.csv`), event impacts on culling (`culling.csv`, see below), and
an empty interactions table (`yield_interactions.csv`) and three-way
table (`three_way.csv`), with a header only. Rename them, or add files
for other outcomes, as you like: the file names carry no meaning.

## Read the files

Pass each file to its argument:

``` r

yield <- cm_read_inputs(
  diseases = file.path(dir, "diseases.csv"),
  associations = file.path(dir, "associations.csv"),
  impacts = file.path(dir, "yield.csv"),
  interactions = file.path(dir, "yield_interactions.csv"),
  three_way = file.path(dir, "three_way.csv")
)
yield
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 (3 with an association, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Estimands: crude: 3
#>   Uncertain inputs (with a distribution): 4
deconflate(yield, n_draws = 0)
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0
```

Empty tables (here the interactions and three-way terms) are the same as
leaving the argument out. Some values in the template have a
distribution, so
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
would by default also run 1000 draws for intervals; `n_draws = 0` gives
point estimates only (see Uncertainty below).

## Associations

### Unknown pairs and odds ratios of 1

A pair of diseases with no row in the association table is **unknown**.
Then
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
uses the global method, which fills in the association of unknown pairs
from the other associations (the maximum-entropy distribution of disease
combinations), and says so in a note. To state that two diseases are
unrelated, give the pair an odds ratio of 1 (with no distribution).

Here the MET:SCK row is left out:

``` r

a <- read.csv(file.path(dir, "associations.csv"))
a[, c("disease1", "disease2", "value", "measure")]
#>   disease1 disease2 value measure
#> 1      LAM      SCK  2.01      OR
#> 2      MET      SCK  1.94      OR
#> 3      MET      LAM  6.10      OR
y_unknown <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                            associations = a[-2, ],
                            impacts = file.path(dir, "yield.csv"))
y_unknown
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 (2 with an association, 1 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Estimands: crude: 3
#>   Uncertain inputs (with a distribution): 4
r_unknown <- deconflate(y_unknown, n_draws = 0)
r_unknown
#> <cm_result> method: global; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted  change
#>      LAM 4.81    2.696 -0.4396
#>      SCK 8.40    7.935 -0.0554
#>      MET 5.61    3.964 -0.2935
#> 
#> Raw sum: 4.966; adjusted total: 4.095
#> Diagnostics: residual 8.88e-16, condition number 1.99, sign changes 0
#> Unknown pairs: 1 without an association; the global fit gave it an odds ratio of 1.33 (see $unknown_pairs).
#> 
#> Notes:
#> * The global method was used because of 1 pair without an association
#>   (unknown).
r_unknown$unknown_pairs
#>   disease1 disease2 fitted_or
#> 1      SCK      MET  1.330169
```

`$unknown_pairs` lists the pairs without an association and the odds
ratio the global fit gave them. With an odds ratio of 1 for MET:SCK
instead, every pair has an association and the simultaneous method is
used:

``` r

a_indep <- a
a_indep$value[2] <- 1
y_indep <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = a_indep,
                          impacts = file.path(dir, "yield.csv"))
deconflate(y_indep, n_draws = 0)
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.557 -0.46844
#>      SCK 8.40    8.057 -0.04079
#>      MET 5.61    4.566 -0.18604
#> 
#> Raw sum: 4.966; adjusted total: 4.167
#> Diagnostics: residual 8.88e-16, condition number 2.01, sign changes 0
```

The two are different statements: an unknown pair is filled in from the
other associations, while an odds ratio of 1 is an estimate of no
association. Some published analyses treated unlisted pairs as
independent; to do the same, add those pairs with an odds ratio of 1.

### Measures

The measures are `OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 \|
disease2)) and `phi`. Files written for earlier versions may use the
measures `unknown`, `independent` or `table`, or the count columns
`n11`, `n10`, `n01` and `n00`; these are reported with what to do
instead:

``` r

cm_check_inputs(
  diseases = data.frame(id = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease1 = c("d1", "d1", "d2"), disease2 = c("d2", "d3", "d3"),
                            value = c(NA, 1, NA), measure = c("unknown", "independent", "table")),
  impacts = data.frame(disease = c("d1", "d2", "d3"), value = c(2.5, 5, 7.5))
)
#> Found 3 problem(s) in the inputs:
#>   associations, row 1, column 'measure': measure = 'unknown' is no longer used: leave the pair out. Pairs without a row are unknown, and the global model fills in their association.
#>   associations, row 2, column 'measure': measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).
#>   associations, row 3, column 'measure': measure = 'table' is no longer used: compute the odds ratio from the counts, n11 * n00 / (n10 * n01), and enter it with measure OR.
```

For a study’s 2x2 table, compute the odds ratio yourself,
`n11 * n00 / (n10 * n01)` (both diseases, only the first, only the
second, neither), and enter it with measure `OR`. If a cell is zero,
decide how to handle it (e.g. add 0.5 to every cell) before computing
the odds ratio.

### Adjusted association measures

The 2x2 algebra needs marginal (crude) association measures. An odds
ratio from a multivariable model is conditional on its covariates, so it
is rejected unless you choose to use it as if it were marginal, which is
an approximation:

``` r

d12 <- cm_diseases(c("d1", "d2"), c(0.10, 0.15))
a_adj <- cm_associations("d1", "d2", 1.8, adjusted = TRUE, adjusted_for = "parity; herd")
tryCatch(cm_population(d12, a_adj), deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "Covariate-adjusted associations (d1:d2) are not marginal 2x2 associations. Use crude measures, or set adjusted_associations = 'use_as_marginal' to use them as an approximation."
cm_population(d12, a_adj, adjusted_associations = "use_as_marginal")
#> <cm_population>
#>   Diseases: 2 (d1, d2)
#>   Disease pairs: 1 (1 with an association, 0 unknown)
#>   Covariate-adjusted association measures used as marginal: 1
```

In files, set the column `adjusted` to `TRUE` and call
`cm_read_inputs(..., adjusted_associations = "use_as_marginal")`; each
such row is then reported as a note.

## Several outcomes

An impact table holds one outcome, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
adjusts one impact table per call. For several outcomes, keep one impact
table per outcome and repeat the call; the disease and association
tables are shared. The installed folder `five_diseases` has impact
tables for milk yield (in percent) and calving interval (in days), among
others:

``` r

five <- system.file("extdata", "five_diseases", package = "deconflate")
list.files(five)
#>  [1] "associations.csv"             "culling.csv"                 
#>  [3] "diseases.csv"                 "impacts_calving_interval.csv"
#>  [5] "impacts_welfare.csv"          "impacts_yield.csv"           
#>  [7] "interactions_welfare.csv"     "README.md"                   
#>  [9] "run_all_features.R"           "three_way.csv"
outcomes <- c(yield = "impacts_yield.csv", calving_interval = "impacts_calving_interval.csv")
results <- lapply(outcomes, function(f) {
  m <- cm_read_inputs(diseases = file.path(five, "diseases.csv"),
                      associations = file.path(five, "associations.csv"),
                      three_way = file.path(five, "three_way.csv"),
                      impacts = file.path(five, f))
  deconflate(m, n_draws = 0)
})
results$calving_interval
#> <cm_result> method: global; calving interval increase [days]
#> 
#>  disease raw adjusted  change        estimand
#>      LAM  12   10.280 -0.1434           crude
#>      MAS   6    4.410 -0.2651 adjusted_linear
#>      MET  18   15.616 -0.1324           crude
#>      SCK   4    4.000  0.0000 adjusted_linear
#>       RP  10    5.067 -0.4933           crude
#> 
#> Raw sum: 8.355; adjusted total: 6.978
#> Diagnostics: residual 3.55e-15, condition number 1.67, sign changes 0
#> 
#> Notes:
#> * The global method was used because of three-way terms.
```

The results can then be combined in base R, e.g. into one long table:

``` r

do.call(rbind, lapply(names(results), function(o) {
  data.frame(outcome = o, units = results[[o]]$units,
             results[[o]]$adjusted[, c("disease", "raw", "adjusted")])
}))
#>             outcome      units disease  raw  adjusted
#> 1             yield % of yield     LAM  4.8  4.139274
#> 2             yield % of yield     MAS  3.3  2.355802
#> 3             yield % of yield     MET  5.6  4.514058
#> 4             yield % of yield     SCK  2.5  1.366295
#> 5             yield % of yield      RP  4.2  2.665420
#> 6  calving_interval       days     LAM 12.0 10.279625
#> 7  calving_interval       days     MAS  6.0  4.409635
#> 8  calving_interval       days     MET 18.0 15.616227
#> 9  calving_interval       days     SCK  4.0  4.000000
#> 10 calving_interval       days      RP 10.0  5.067130
```

## Event impacts

Impacts on an event such as culling or death are event impacts: hazard
ratios, rate ratios, risk ratios, odds ratios or risk differences. Their
table has a `measure` column (`HR`, `rate_ratio`, `RR`, `OR` or `RD`,
which can be mixed, one row per disease) and an `estimand`
(`snapshot_crude`, or `snapshot_stratified` with `adjusted_for`). They
are adjusted with `event_model = TRUE`, which needs the overall risk of
the event in the population over the period; risk ratios, odds ratios
and risk differences must refer to the same period. The template’s
`culling.csv`:

``` r

read.csv(file.path(dir, "culling.csv"))[, c("disease", "value", "measure", "estimand")]
#>   disease value measure       estimand
#> 1     LAM  1.74      HR snapshot_crude
#> 2     SCK  1.92      HR snapshot_crude
#> 3     MET  1.45      RR snapshot_crude
culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = file.path(dir, "associations.csv"),
                          impacts = file.path(dir, "culling.csv"))
deconflate(culling, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure  raw adjusted_hr
#>      LAM      HR 1.74       1.514
#>      SCK      HR 1.92       1.782
#>      MET      RR 1.45       1.201
#> 
#> Overall risk 0.25; disease-free risk 0.1767; attributable to disease 0.07328 (29.3% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable   share
#>      LAM     0.022295 0.30424
#>      SCK     0.047051 0.64205
#>      MET     0.003936 0.05371
#> 
#> Diagnostics: residual 2.22e-16, condition number 2.2, sign changes 0
```

The adjusted hazard ratios are each disease’s own effect; the result
also gives the risk attributable to disease and its split by disease. If
an adjusted hazard ratio falls on the other side of 1 from its raw
estimate,
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
warns: the raw estimate is weaker than the disease’s associations alone
would produce, so check the estimate and its source population.

A table with a `measure` column needs `event_model = TRUE`, and
`event_model = TRUE` needs a `measure` column:

``` r

tryCatch(deconflate(culling, n_draws = 0), deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "These are event impacts (the impact table has a measure column): use event_model = TRUE, with overall_risk."
```

See
[`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md)
for the model and its results.

## Uncertainty

There is no separate uncertainty file: an uncertain value gets its
distribution in its own row, in the columns `dist` and `p1`-`p4`, so
point values and uncertain values sit side by side. A row with an empty
`dist` is a point value. In the template, the LAM prevalence, the
LAM:SCK odds ratio and two yield impacts have a distribution:

``` r

read.csv(file.path(dir, "diseases.csv"))[, c("id", "value", "type", "dist", "p1", "p2")]
#>    id value           type dist    p1     p2
#> 1 LAM  0.25     prevalence beta 78.29 227.42
#> 2 SCK  0.48 incidence_rate         NA     NA
#> 3 MET  0.10     prevalence         NA     NA
read.csv(file.path(dir, "associations.csv"))[, c("disease1", "disease2", "value", "dist", "p1", "p2", "p3")]
#>   disease1 disease2 value   dist   p1  p2 p3
#> 1      LAM      SCK  2.01 normal 2.01 0.2  0
#> 2      MET      SCK  1.94          NA  NA NA
#> 3      MET      LAM  6.10          NA  NA NA
```

[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
keeps them in the model, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
draws them (1000 draws by default) and adds 95% intervals to the
results:

``` r

deconflate(yield, n_draws = 200, seed = 1)
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted lower  upper   change
#>      LAM 4.81    2.871 0.600  4.697 -0.40320
#>      SCK 8.40    7.819 5.475 10.065 -0.06915
#>      MET 5.61    3.170 2.302  4.010 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015 (95% interval 3.111 to 4.932)
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
```

The central estimate always uses the point values (`value`); the draws
give the intervals. Values are drawn on the scale they were entered on:

- a disease’s `value` (e.g. a prevalence or an incidence rate);
- an association’s `value`, on its own measure (e.g. an odds ratio);
- a three-way `ratio`;
- an impact or an interaction, in the units of its table (for event
  impacts, on the scale of the measure, e.g. a hazard ratio).

`dist` and the parameters `p1`-`p4`:

| `dist` | `p1` | `p2` | `p3` | `p4` |
|----|----|----|----|----|
| `fixed` | value |  |  |  |
| `normal` | mean | sd | lower bound (optional) | upper bound (optional) |
| `lognormal` | meanlog | sdlog |  |  |
| `lognormal_ci` | estimate | lower CI | upper CI | level (default 0.95) |
| `beta` | shape1 | shape2 | min (default 0) | max (default 1) |
| `pert` | min | mode | max | lambda (default 4) |
| `pert_mean` | min | mean | max | lambda (default 4) |
| `uniform` | min | max |  |  |

The checks report a distribution with missing or invalid parameters and
parameters without a `dist`; a point value outside its distribution’s
range is noted. In the model, the distributions are kept in
`$distributions`, keyed `prob:<disease>`, `assoc:<d1>:<d2>`,
`three:<d1>:<d2>:<d3>`, `impact:<disease>` and `inter:<d1>:<d2>`, the
same keys as `cm_model(distributions = )`:

``` r

names(yield$distributions)
#> [1] "prob:LAM"      "assoc:LAM:SCK" "impact:LAM"    "impact:SCK"
```

[`cm_dist_table()`](https://rasmussenphilip.github.io/deconflate/reference/cm_dist_table.md)
builds distributions from a table with columns `key`, `dist` and
`p1`-`p4`. See
[`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md)
for the draws, their intervals and rejected draws.

## When something is wrong

All problems are listed with their table, row (counting from the first
row below the header) and column, so that they can be fixed in one pass.
[`cm_check_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_check_inputs.md)
returns them without stopping:

``` r

cm_check_inputs(
  diseases = data.frame(id = c("LAM", "SCK", "SCK"), value = c(0.25, 1.4, 0.3)),
  associations = data.frame(disease1 = "LAM", disease2 = "MET", value = 2.0),
  impacts = data.frame(disease = c("LAM", "SCK"), value = c("4.81", "8,40"))
)
#> Found 4 problem(s) in the inputs:
#>   diseases, row 3, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 2, column 'value': Gives a probability of 1.4; it must be strictly between 0 and 1.
#>   associations, row 1, column 'disease2': Unknown disease 'MET' (not in the diseases table).
#>   impacts, row 2, column 'value': '8,40' is not a number.
```

[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
stops with the same list (an error of class `deconflate_input_problems`)
when there are errors, and prints notes otherwise.

The installed folder `example_with_errors` has a disease table, an
association table, an additive impact table (`impacts_yield.csv`) and an
event impact table (`culling.csv`), each with mistakes explained in its
`note` column: a missing distribution parameter, a prevalence above 1, a
duplicate disease or pair, a missing association value, an unknown
disease, a covariate-adjusted odds ratio, measures and count columns
from earlier versions, a decimal comma, two units in one table, an
unsupported distribution, and a point value outside its distribution:

``` r

err <- system.file("extdata", "example_with_errors", package = "deconflate")
list.files(err)
#> [1] "associations.csv"  "culling.csv"       "diseases.csv"     
#> [4] "impacts_yield.csv"
cm_check_inputs(diseases = file.path(err, "diseases.csv"),
                associations = file.path(err, "associations.csv"),
                impacts = file.path(err, "impacts_yield.csv"))
#> Found 18 problem(s) in the inputs:
#>   diseases, row 5, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   diseases, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   diseases, row 4, column 'dist': Parameters p1-p4 are given without a distribution (dist).
#>   associations: Contingency-table counts (n11, n10, n01, n00) are no longer read: compute each table's odds ratio, n11 * n00 / (n10 * n01), and enter it in `value` with measure OR.
#>   associations, row 3, column 'disease2': Unknown disease 'CM' (not in the diseases table).
#>   associations, row 9, column 'disease2': Unknown disease 'MAS' (not in the diseases table).
#>   associations, row 4: Duplicate pair SCK:LAM.
#>   associations, row 7, column 'measure': measure = 'table' is no longer used: compute the odds ratio from the counts, n11 * n00 / (n10 * n01), and enter it with measure OR.
#>   associations, row 8, column 'measure': measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).
#>   associations, row 9, column 'measure': measure = 'unknown' is no longer used: leave the pair out. Pairs without a row are unknown, and the global model fills in their association.
#>   associations, row 5, column 'value': Missing value. Leave the pair out (no row) if its association is unknown; an odds ratio of 1 states that the diseases are unrelated.
#>   associations, row 6, column 'adjusted': A covariate-adjusted measure is not a marginal 2x2 association. Use the crude measure, or set adjusted_associations = 'use_as_marginal' to use it as an approximation.
#>   associations, row 2, column 'dist': PERT needs min <= mode <= max and min < max.
#>   impacts, row 2, column 'value': '8,40' is not a number.
#>   impacts, row 3, column 'adjusted_for': adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.
#>   impacts, column 'units': One units per table is allowed (found: % decrease, days).
#>   impacts, row 3, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
#>   impacts, row 4, column 'dist' (note): The point value 5.61 lies outside the distribution's support [6, 8].
```

The event impact table has a missing estimand, a hazard ratio of 0, and
a stratified estimand without `adjusted_for` (the disease table’s
problems, listed above, are left out here):

``` r

p <- cm_check_inputs(diseases = file.path(err, "diseases.csv"),
                     impacts = file.path(err, "culling.csv"))
p[p$table == "impacts", ]
#> Found 3 problem(s) in the inputs:
#>   impacts, row 2, column 'estimand': Missing estimand: event impacts need snapshot_crude or snapshot_stratified (see ?cm_impacts).
#>   impacts, row 4, column 'adjusted_for': estimand = snapshot_stratified needs adjusted_for (disease ids or all).
#>   impacts, row 3, column 'value': Ratios must be positive.
```

## Tables typed in R

Any table can be a data frame instead of a file, with the same columns,
including the distribution columns:

``` r

typed <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease1 = c("d1", "d2"), disease2 = c("d2", "d3"),
                            value = c(2, 3), measure = "OR",
                            dist = c("lognormal_ci", NA), p1 = c(2, NA), p2 = c(1.4, NA),
                            p3 = c(2.9, NA)),
  impacts = data.frame(disease = c("d1", "d2", "d3"), value = c(2.5, 5, 7.5),
                       units = "% of yield", dist = c("normal", NA, NA),
                       p1 = c(2.5, NA, NA), p2 = c(0.5, NA, NA))
)
typed
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (2 with an association, 1 unknown)
#>   Impacts: (unlabelled) [% of yield] (additive)
#>   Estimands: crude: 3
#>   Uncertain inputs (with a distribution): 2
deconflate(typed, n_draws = 0)$adjusted
#>   disease raw adjusted      change estimand adjusted_for
#> 1      d1 2.5 1.988942 -0.20442313    crude         <NA>
#> 2      d2 5.0 3.404334 -0.31913313    crude         <NA>
#> 3      d3 7.5 6.906626 -0.07911658    crude         <NA>
```

Here the pair d1:d3 has no row, so it is unknown and the global method
is used. The R constructors
([`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md),
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md),
[`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md),
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md),
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md))
build the same model directly.

## Disease ids

Ids are used in keys such as `assoc:d1:d2` and in lists such as
`"d1; d2"`, so they must not contain `|`, `;` or `:`, nor leading or
trailing spaces. `"all"` is reserved (it means every other disease in
`adjusted_for`):

``` r

tryCatch(cm_diseases(c("LAM", "SCK:early"), c(0.25, 0.40)),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "Disease ids must not contain '|', ';' or ':' or surrounding spaces (check: SCK:early)."
tryCatch(cm_diseases(c("LAM", "all"), c(0.25, 0.40)),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "'all' is reserved and cannot be a disease id."
```

## Fuller examples

The folder `five_diseases` uses every feature on five diseases: all
probability types and association measures, a three-way term, crude and
adjusted impacts, interactions (`interactions_welfare.csv`, for the
welfare impacts), event impacts with mixed measures (`culling.csv`), and
every distribution type mixed with point values. Its
`run_all_features.R` script runs every part of the package on it.

The folder `global_dairy_2024` holds the inputs of Rasmussen et
al. (2024) for 12 diseases: `diseases.csv`, `associations.csv` (the 38
odds ratios of Table 3), `impacts_yield.csv`, `impacts_fertility.csv`
and `culling.csv` (hazard ratios entered as `snapshot_crude`), with the
input distributions of the paper’s analysis in the `dist` columns. Pairs
not in Table 3 have no row, so they are unknown:

``` r

gd <- system.file("extdata", "global_dairy_2024", package = "deconflate")
cm_read_inputs(diseases = file.path(gd, "diseases.csv"),
               associations = file.path(gd, "associations.csv"),
               impacts = file.path(gd, "impacts_yield.csv"))
#> <cm_model>
#>   Diseases: 12 (CK, CM, DA, DYS, LAM, MET, MF, OC, PTB, RP, SCK, SCM)
#>   Disease pairs: 66 (38 with an association, 28 unknown)
#>   Impacts: milk yield loss [% decrease] (additive)
#>   Estimands: crude: 12
#>   Uncertain inputs (with a distribution): 46
```

[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)
returns the same inputs as a model, and
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)
compares it with the published analysis, which treated these pairs as
independent.
