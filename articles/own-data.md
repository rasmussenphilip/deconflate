# Using your own data

``` r

library(deconflate)
```

Inputs are three tables: the diseases, their associations and one impact
table. Each table can be a CSV file (with any name) or a data frame
typed in R.
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
| `diseases` | one row per disease | `disease`, `value` |
| `associations` | association estimates, one row per estimate | `disease`, `associated_with`, `value` |
| `impacts` | the raw impact estimates of one outcome, one row per estimate | `disease`, `value` |

`diseases` and `impacts` are required.
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
also needs at least one association; without any, the sensitivity tools
([`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md))
show how much associations could matter. The optional columns are
described in
[`?cm_read_inputs`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):

- `diseases`: `measure` (`prevalence` (default), `probability` or
  `incidence_rate`), `label`, `time_horizon`, `reference_population`;
- `associations`: `among`, `adjusted_for`, `measure` (`OR` (default),
  `RR`, `RD`, `phi`, `prevalence` or `OR_ratio`), `study`;
- `impacts`: `among`, `adjusted_for`, `study`, and one `label` and one
  `units` per table; a `measure` column makes them event impacts (see
  below).

Every table can also have the free-text columns `source` and `note`, and
the columns `dist` and `p1`-`p4` for uncertain values (see Uncertainty
below). Column names are not case-sensitive.

A row is about a disease or a combination of diseases (ids separated by
`";"`), holds in the whole population or within a subgroup (`among`: a
disease, `"no "` and a disease, or `"no others"`), and is crude unless
`adjusted_for` names other diseases (ids separated by `";"`, or
`"all"`). Several rows can describe the same thing, e.g. estimates from
several studies.
[`?cm_associations`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
and
[`?cm_impacts`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
define every kind of row;
[`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)
shows them at work.

Every disease needs a row of its own in the impact table (use 0 for no
impact). Additive impacts can be in any units, the same for every row;
results come back in those units.

## Start from the template

[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set of files to a folder:

``` r

dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 4 files to /tmp/RtmpSIAJcF/my-inputs
list.files(dir)
#> [1] "associations.csv" "culling.csv"      "diseases.csv"     "yield.csv"
read.csv(file.path(dir, "yield.csv"))
#>   disease among adjusted_for value           label      units study   dist   p1
#> 1     LAM    NA           NA  4.81 milk yield loss % of yield    NA normal 4.81
#> 2     SCK    NA           NA  8.40 milk yield loss % of yield    NA normal 8.40
#> 3     MET    NA           NA  5.61 milk yield loss % of yield    NA          NA
#>     p2 p3 p4              source             note
#> 1 0.87 NA NA Illustrative values normal(mean, sd)
#> 2 1.19 NA NA Illustrative values normal(mean, sd)
#> 3   NA NA NA Illustrative values
```

The files are three diseases (`diseases.csv`), their odds ratios
(`associations.csv`), additive impacts on milk yield in percent
(`yield.csv`) and event impacts on culling (`culling.csv`, see below).
Rename them, or add files for other outcomes, as you like: the file
names carry no meaning.

## Read the files

Pass each file to its argument:

``` r

yield <- cm_read_inputs(
  diseases = file.path(dir, "diseases.csv"),
  associations = file.path(dir, "associations.csv"),
  impacts = file.path(dir, "yield.csv")
)
yield
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Association rows: 3 (OR: 3)
#>   Disease pairs: 3 (3 with a crude pair row, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Impact rows: 3 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 4
deconflate(yield, n_draws = 0)
#> <cm_result> pairwise shortcut; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 0.00e+00, condition number 2, sign changes 0
```

Some values in the template have a distribution, so
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
would by default also run 1000 draws for intervals; `n_draws = 0` gives
point estimates only (see Uncertainty below).

## Associations

### Unknown pairs and odds ratios of 1

A pair of diseases that no association row mentions is **unknown**.
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
fills in its association from the other associations (the
maximum-entropy distribution of disease combinations), and lists it in
`$unknown_pairs`. To state that two diseases are unrelated, give the
pair an odds ratio of 1 (with no distribution).

Here the MET:SCK row is left out:

``` r

a <- read.csv(file.path(dir, "associations.csv"))
a[, c("disease", "associated_with", "value", "measure")]
#>   disease associated_with value measure
#> 1     LAM             SCK  2.01      OR
#> 2     MET             SCK  1.94      OR
#> 3     MET             LAM  6.10      OR
y_unknown <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                            associations = a[-2, ],
                            impacts = file.path(dir, "yield.csv"))
y_unknown
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Association rows: 2 (OR: 2)
#>   Disease pairs: 3 (2 with a crude pair row, 1 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Impact rows: 3 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 4
r_unknown <- deconflate(y_unknown, n_draws = 0)
r_unknown
#> <cm_result> full model; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted  change
#>      LAM 4.81    2.696 -0.4396
#>      SCK 8.40    7.935 -0.0554
#>      MET 5.61    3.964 -0.2935
#> 
#> Raw sum: 4.966; adjusted total: 4.095
#> Diagnostics: residual 0.00e+00, condition number 1.99, sign changes 0
#> Unknown pairs: 1 without an association; the fit gave it an odds ratio of 1.33 (see $unknown_pairs).
r_unknown$unknown_pairs
#>   disease1 disease2 fitted_or
#> 1      SCK      MET  1.330169
```

`$unknown_pairs` lists the pairs without an association and the odds
ratio the fit gave them. With an odds ratio of 1 for MET:SCK instead,
every pair has an association, and the pairwise tables are used
directly:

``` r

a_indep <- a
a_indep$value[2] <- 1
y_indep <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = a_indep,
                          impacts = file.path(dir, "yield.csv"))
deconflate(y_indep, n_draws = 0)
#> <cm_result> pairwise shortcut; milk yield loss [% of yield]
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

The measures are `OR` (default), `RR`, `RD`, `phi`, `prevalence` (the
prevalence of `disease` among animals with `associated_with`, or within
`among`) and `OR_ratio` (how much a third disease, in `among`, changes
an odds ratio). Files written for earlier versions may use the columns
`disease1` and `disease2`, the measures `cond_prob`, `unknown`,
`independent` or `table`, the count columns `n11`, `n10`, `n01` and
`n00`, or a separate three-way or interactions table; these are reported
with what to do instead:

``` r

cm_check_inputs(
  diseases = data.frame(disease = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease = c("d1", "d1", "d2"), associated_with = c("d2", "d3", "d3"),
                            value = c(0.3, 1, NA), measure = c("cond_prob", "independent", "table")),
  impacts = data.frame(disease = c("d1", "d2", "d3"), value = c(2.5, 5, 7.5))
)
#> Found 3 problem(s) in the inputs:
#>   associations, row 1, column 'measure': measure = 'cond_prob' is now 'prevalence': P(disease | associated_with) is the prevalence of `disease` among animals with `associated_with`.
#>   associations, row 2, column 'measure': measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).
#>   associations, row 3, column 'measure': measure = 'table' is no longer used: compute the odds ratio from the counts, n11 * n00 / (n10 * n01), and enter it with measure OR.
```

For a study’s 2x2 table, compute the odds ratio yourself,
`n11 * n00 / (n10 * n01)` (both diseases, only the first, only the
second, neither), and enter it with measure `OR`. If a cell is zero,
decide how to handle it (e.g. add 0.5 to every cell) before computing
the odds ratio.

### Adjusted odds ratios

An odds ratio from a logistic regression that included other diseases of
the model is entered with those diseases in `adjusted_for`; the model
matches it as the regression coefficient it is, not as a crude odds
ratio. Covariates that are not diseases of the model (parity, herd, …)
cannot be matched: they are noted, and the row is used as a crude odds
ratio, which is an approximation:

``` r

d12 <- cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20))
a_adj <- cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(1.8, 1.2, 2.5),
                         adjusted_for = c("d3; parity; herd", NA, NA))
r_adj <- deconflate(cm_model(d12, cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5)),
                             associations = a_adj), n_draws = 0)
r_adj$notes
#> [1] "Association row 1: the adjustment for factors outside the model (parity, herd) is ignored; the adjustment for the model's diseases is used."
```

### Several rows for one item

Two studies of the same odds ratio are two rows; they are pooled (random
effects by default, see
[`?deconflate`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
Several impact rows of one disease (e.g. a crude and an adjusted
estimate) are fitted together by weighted least squares. The
five-disease example below has both; `$heterogeneity` reports how well
the rows agree.

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
#> [1] "associations.csv"             "culling.csv"                 
#> [3] "diseases.csv"                 "impacts_calving_interval.csv"
#> [5] "impacts_welfare.csv"          "impacts_yield.csv"           
#> [7] "README.md"                    "run_all_features.R"
outcomes <- c(yield = "impacts_yield.csv", calving_interval = "impacts_calving_interval.csv")
results <- lapply(outcomes, function(f) {
  m <- cm_read_inputs(diseases = file.path(five, "diseases.csv"),
                      associations = file.path(five, "associations.csv"),
                      impacts = file.path(five, f))
  deconflate(m, n_draws = 0)
})
results$calving_interval
#> <cm_result> full model; calving interval increase [days]
#> 
#>  disease raw adjusted     change rows
#>      LAM  12   10.206 -1.495e-01    1
#>      MAS  NA    4.768         NA    2
#>      MET  18   14.774 -1.792e-01    1
#>      SCK   4    4.000  2.220e-16    1
#>       RP  10   10.002  2.042e-04    1
#> 
#> Raw sum: NA (several rows per disease); adjusted total: 7.265
#> Diagnostics: residual 4.44e-16, condition number 1.67, sign changes 0
#> Pooled or over-determined items: 2 (see $heterogeneity and $estimates)
#> 
#> Notes:
#> * Several association rows were pooled for 1 item (random effects; see
#>   $heterogeneity).
#> * The impact rows over-determine the impacts (1 more rows than unknowns):
#>   fitted by weighted least squares (Q = 0.312 on 1 df, p = 0.576; see
#>   $estimates). Rows without uncertainty are exact.
```

The results can then be combined in base R, e.g. into one long table:

``` r

do.call(rbind, lapply(names(results), function(o) {
  data.frame(outcome = o, units = results[[o]]$units,
             results[[o]]$adjusted[, c("disease", "rows", "adjusted")])
}))
#>             outcome      units disease rows  adjusted
#> 1             yield % of yield     LAM    1  4.103601
#> 2             yield % of yield     MAS    1  2.284185
#> 3             yield % of yield     MET    1  4.520276
#> 4             yield % of yield     SCK    1  1.380768
#> 5             yield % of yield      RP    1  2.668537
#> 6  calving_interval       days     LAM    1 10.206306
#> 7  calving_interval       days     MAS    2  4.768205
#> 8  calving_interval       days     MET    1 14.774125
#> 9  calving_interval       days     SCK    1  4.000000
#> 10 calving_interval       days      RP    1 10.002042
```

## Event impacts

Impacts on an event such as culling or death are event impacts: hazard
ratios, rate ratios, risk ratios, odds ratios or risk differences. Their
table has a `measure` column (`HR`, `rate_ratio`, `RR`, `OR` or `RD`,
which can be mixed); rows are crude or adjusted (`adjusted_for`), about
diseases or combinations, as for additive impacts. They are adjusted
with `event_model = TRUE`, which needs the overall risk of the event in
the population over the period; risk ratios, odds ratios and risk
differences must refer to the same period. The template’s `culling.csv`:

``` r

read.csv(file.path(dir, "culling.csv"))[, c("disease", "adjusted_for", "measure", "value")]
#>   disease adjusted_for measure value
#> 1     LAM           NA      HR  1.74
#> 2     SCK           NA      HR  1.92
#> 3     MET           NA      RR  1.45
culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = file.path(dir, "associations.csv"),
                          impacts = file.path(dir, "culling.csv"))
deconflate(culling, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
#> <cm_event_result> culling; snapshot
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
#> Diagnostics: residual 3.89e-16, condition number 2.2, sign changes 0
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

read.csv(file.path(dir, "diseases.csv"))[, c("disease", "measure", "value", "dist", "p1", "p2")]
#>   disease        measure value dist    p1     p2
#> 1     LAM     prevalence  0.25 beta 78.29 227.42
#> 2     SCK incidence_rate  0.48         NA     NA
#> 3     MET     prevalence  0.10         NA     NA
read.csv(file.path(dir, "associations.csv"))[, c("disease", "associated_with", "value", "dist", "p1", "p2", "p3")]
#>   disease associated_with value   dist   p1  p2 p3
#> 1     LAM             SCK  2.01 normal 2.01 0.2  0
#> 2     MET             SCK  1.94          NA  NA NA
#> 3     MET             LAM  6.10          NA  NA NA
```

[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
keeps them in the model, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
draws them (1000 draws by default) and adds 95% intervals to the
results:

``` r

deconflate(yield, n_draws = 200, seed = 1)
#> <cm_result> pairwise shortcut; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted lower  upper   change
#>      LAM 4.81    2.871 0.600  4.697 -0.40320
#>      SCK 8.40    7.819 5.475 10.065 -0.06915
#>      MET 5.61    3.170 2.302  4.010 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015 (95% interval 3.111 to 4.932)
#> Diagnostics: residual 0.00e+00, condition number 2, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
```

The central estimate always uses the point values (`value`); the draws
give the intervals. Values are drawn on the scale they were entered on:

- a disease’s `value` (e.g. a prevalence or an incidence rate);
- an association row’s `value`, on its own measure (e.g. an odds ratio);
- an impact row’s `value`, in the units of its table (for event impacts,
  on the scale of the measure, e.g. a hazard ratio).

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
`$distributions`, keyed `prob:<disease>`, `associations:<row>` and
`impacts:<row>` (rows numbered from 1 in each table), the same keys as
`cm_model(distributions = )`, which also accepts `assoc:<d1>:<d2>` and
`impact:<disease>` when they name one row:

``` r

names(yield$distributions)
#> [1] "prob:LAM"       "associations:1" "impacts:1"      "impacts:2"
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
  diseases = data.frame(disease = c("LAM", "SCK", "SCK"), value = c(0.25, 1.4, 0.3)),
  associations = data.frame(disease = "LAM", associated_with = "MET", value = 2.0),
  impacts = data.frame(disease = c("LAM", "SCK"), value = c("4.81", "8,40"))
)
#> Found 4 problem(s) in the inputs:
#>   diseases, row 3, column 'disease': Duplicate disease 'SCK': the diseases table has one row per disease (rows about several diseases or subgroups belong in the associations table).
#>   diseases, row 2, column 'value': Gives a probability of 1.4; it must be strictly between 0 and 1.
#>   associations, row 1, column 'associated_with': unknown disease in `associated_with`: MET.
#>   impacts, row 2, column 'value': '8,40' is not a number.
```

[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
stops with the same list (an error of class `deconflate_input_problems`)
when there are errors, and prints notes otherwise.

The installed folder `example_with_errors` has a disease table, an
association table, an additive impact table (`impacts_yield.csv`) and an
event impact table (`culling.csv`), each with mistakes explained in its
`note` column: a missing distribution parameter, a prevalence above 1, a
duplicate disease, a combination in the disease table, a missing
association value, an unknown disease, an odds ratio adjusted only for
factors outside the model, measures from earlier versions, rows that are
not well formed (a subgroup naming the row’s own disease, an adjusted
risk ratio, an `OR_ratio` without a subgroup), a decimal comma, two
units in one table, an unsupported distribution, and a point value
outside its distribution:

``` r

err <- system.file("extdata", "example_with_errors", package = "deconflate")
list.files(err)
#> [1] "associations.csv"  "culling.csv"       "diseases.csv"     
#> [4] "impacts_yield.csv"
cm_check_inputs(diseases = file.path(err, "diseases.csv"),
                associations = file.path(err, "associations.csv"),
                impacts = file.path(err, "impacts_yield.csv"))
#> Found 17 problem(s) in the inputs:
#>   diseases, row 5, column 'disease': Duplicate disease 'SCK': the diseases table has one row per disease (rows about several diseases or subgroups belong in the associations table).
#>   diseases, row 6, column 'disease': Disease id 'LAM;SCK' must not contain '|', ';' or ':' (a combination of diseases belongs in the associations table).
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   diseases, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   diseases, row 4, column 'dist': Parameters p1-p4 are given without a distribution (dist).
#>   associations, row 6, column 'measure': measure = 'cond_prob' is now 'prevalence': P(disease | associated_with) is the prevalence of `disease` among animals with `associated_with`.
#>   associations, row 7, column 'measure': measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).
#>   associations, row 4, column 'value': Missing value. Leave the pair out (no row) if its association is unknown; an odds ratio of 1 states that the diseases are unrelated.
#>   associations, row 3, column 'associated_with': unknown disease in `associated_with`: CM.
#>   associations, row 8, column 'among': `among` names LAM, which the row is about.
#>   associations, row 9: an adjusted association must be an odds ratio between two single diseases (measure OR).
#>   associations, row 10, column 'among': OR_ratio needs one disease in `among` (present, e.g. 'SCK', or absent, 'no SCK').
#>   associations, row 5, column 'adjusted_for' (note): 'parity', 'herd' are not diseases of the model: treated as factors outside the model (ignored; the row is used as crude).
#>   associations, row 2, column 'dist': PERT needs min <= mode <= max and min < max.
#>   impacts, row 2, column 'value': '8,40' is not a number.
#>   impacts, row 5, column 'disease': unknown disease in `disease`: CM.
#>   impacts, column 'units': One unit per table is allowed (found: % decrease, days).
#>   impacts, row 3, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
#>   impacts, row 4, column 'dist' (note): The point value 5.61 lies outside the distribution's support [6, 8].
```

The event impact table has the `estimand` column of earlier versions
(now ignored, with a note), a row without a measure, a hazard ratio of 0
and a risk difference above 1 (the disease table’s problems, listed
above, are left out here):

``` r

p <- cm_check_inputs(diseases = file.path(err, "diseases.csv"),
                     impacts = file.path(err, "culling.csv"))
p[p$table == "impacts", ]
#> Found 3 problem(s) in the inputs:
#>   impacts, column 'estimand' (note): The estimand column is no longer used: an estimate is crude when `adjusted_for` is empty and adjusted for the diseases it names otherwise. Column ignored.
#>   impacts, row 2, column 'measure': Missing measure: event impacts need a measure in every row (HR, rate_ratio, RR, OR, RD).
#>   impacts, row 3, column 'value': Ratios must be positive.
#>   impacts, row 4, column 'value': A risk difference must lie between -1 and 1.
```

## Tables typed in R

Any table can be a data frame instead of a file, with the same columns,
including the distribution columns:

``` r

typed <- cm_read_inputs(
  diseases = data.frame(disease = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease = c("d1", "d2"), associated_with = c("d2", "d3"),
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
#>   Association rows: 2 (OR: 2)
#>   Disease pairs: 3 (2 with a crude pair row, 1 unknown)
#>   Impacts: (unlabelled) [% of yield] (additive)
#>   Impact rows: 3 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 2
deconflate(typed, n_draws = 0)$adjusted
#>   disease    kind raw adjusted      change rows
#> 1      d1 disease 2.5 1.988942 -0.20442313    1
#> 2      d2 disease 5.0 3.404334 -0.31913313    1
#> 3      d3 disease 7.5 6.906626 -0.07911658    1
```

Here the pair d1:d3 has no row, so it is unknown and filled in by the
fit. The R constructors
([`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md),
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md),
[`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md),
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md),
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md))
build the same model directly.

## Disease ids

Ids are used in keys such as `assoc:d1:d2` and in lists such as
`"d1; d2"`, so they must not contain `|`, `;` or `:`, nor leading or
trailing spaces. `"all"` and `"others"` are reserved, and an id cannot
start with `"no "` (these are used in `adjusted_for` and `among`):

``` r

tryCatch(cm_diseases(c("LAM", "SCK:early"), c(0.25, 0.40)),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "Disease ids must not contain '|', ';' or ':' or surrounding spaces (check: SCK:early)."
tryCatch(cm_diseases(c("LAM", "all"), c(0.25, 0.40)),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "'all' and 'others' are reserved and cannot be disease ids."
```

## Fuller examples

The folder `five_diseases` uses every feature on five diseases: all
probability types and association measures, two studies of one odds
ratio (pooled), an `OR_ratio` row, crude, adjusted and subgroup impact
rows, several rows for one disease, interactions (in
`impacts_welfare.csv`), event impacts with mixed measures and an
interaction (`culling.csv`), and every distribution type mixed with
point values. Its `run_all_features.R` script runs every part of the
package on it.

The folder `global_dairy_2024` holds the inputs of Rasmussen et
al. (2024) for 12 diseases: `diseases.csv`, `associations.csv` (the 38
odds ratios of Table 3), `impacts_yield.csv`, `impacts_fertility.csv`
and `culling.csv` (hazard ratios entered as crude rows), with the input
distributions of the paper’s analysis in the `dist` columns. Pairs not
in Table 3 have no row, so they are unknown:

``` r

gd <- system.file("extdata", "global_dairy_2024", package = "deconflate")
cm_read_inputs(diseases = file.path(gd, "diseases.csv"),
               associations = file.path(gd, "associations.csv"),
               impacts = file.path(gd, "impacts_yield.csv"))
#> <cm_model>
#>   Diseases: 12 (CK, CM, DA, DYS, LAM, MET, MF, OC, PTB, RP, SCK, SCM)
#>   Association rows: 38 (OR: 38)
#>   Disease pairs: 66 (38 with a crude pair row, 28 unknown)
#>   Impacts: milk yield loss [% decrease] (additive)
#>   Impact rows: 12 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 46
```

[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)
returns the same inputs as a model, and
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)
compares it with the published analysis, which treated these pairs as
independent.
