# Using your own data

``` r

library(deconflate)
```

Inputs can be kept in CSV files (one per table) or typed in R as data
frames.
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
checks every table, reports all problems at once, and returns the
population, one analysis per impact table, a hazard-ratio model if there
are hazard ratios, and a Monte Carlo sampler if there are input
distributions.

## Files

| File | Contents | Required columns |
|----|----|----|
| `diseases.csv` | one row per disease | `id`, `value` |
| `associations.csv` | one row per associated pair | `disease1`, `disease2` |
| `three_way.csv` | optional three-way association scenarios | `disease1`, `disease2`, `disease3`, `ratio` |
| `impacts_<analysis>.csv` | one impact table per analysis, e.g. `impacts_yield.csv` | `disease`, `value` |
| `interactions_<analysis>.csv` | optional impact interactions of that analysis | `disease1`, `disease2`, `value` |
| `uncertainty.csv` | optional input distributions shared by all analyses | `key`, `dist` |
| `uncertainty_<analysis>.csv` | optional input distributions of one analysis | `key`, `dist` |
| `hazard_ratios.csv` | optional culling (or mortality) hazard ratios | `disease`, `value` |

Only `diseases.csv` is required. A single `impacts.csv` is read as one
analysis named `impacts`. All analyses share the diseases, associations
and three-way terms. The optional columns are described in
[`?cm_read_inputs`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):

- `diseases.csv`: `type` (`prevalence`, `probability` or
  `incidence_rate`), `time_horizon`, `reference_population`, `source`;
- `associations.csv`: `value`, `measure` (`OR`, `RR`, `RD`, `cond_prob`,
  `phi`, `table`, `independent` or `unknown`), the counts `n11`, `n10`,
  `n01`, `n00` for `measure = table`, `adjusted`, `adjusted_for`,
  `source`;
- `impacts_<analysis>.csv`: `estimand` (`crude` or `adjusted_linear`),
  `adjusted_for`, `source`, and one `label` and one `units` per table;
- `hazard_ratios.csv`: `estimand` (`crude` or `adjusted`),
  `adjusted_for`, `source`.

Every disease needs a row in every impact table (use 0 for no impact)
and in the hazard-ratio table (use 1 for no effect). Impact values can
be in any units, the same within a table; results come back in those
units.

## Start from the template

[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set of files to a folder:

``` r

dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 8 files to /tmp/RtmpPRbRoQ/my-inputs
list.files(dir)
#> [1] "associations.csv"             "diseases.csv"                
#> [3] "hazard_ratios.csv"            "impacts_calving_interval.csv"
#> [5] "impacts_yield.csv"            "interactions_yield.csv"      
#> [7] "three_way.csv"                "uncertainty.csv"
```

It has two analyses, milk yield (in percent) and calving interval (in
days):

``` r

read.csv(file.path(dir, "impacts_yield.csv"))
#>   disease value estimand adjusted_for           label      units
#> 1     LAM  4.81    crude           NA milk yield loss % of yield
#> 2     SCK  8.40    crude           NA milk yield loss % of yield
#> 3     MET  5.61    crude           NA milk yield loss % of yield
#>                source
#> 1 Illustrative values
#> 2 Illustrative values
#> 3 Illustrative values
read.csv(file.path(dir, "impacts_calving_interval.csv"))
#>   disease value estimand adjusted_for                     label units
#> 1     LAM    12    crude           NA calving interval increase  days
#> 2     SCK     4    crude           NA calving interval increase  days
#> 3     MET    18    crude           NA calving interval increase  days
#>                source
#> 1 Illustrative values
#> 2 Illustrative values
#> 3 Illustrative values
```

## Read the files

``` r

inp <- cm_read_inputs(dir = dir)
inp
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
#>   Analyses: calving_interval, yield
#>   Hazard ratios: yes (use $hr_model with deconflate_hr())
#>   Uncertain inputs: batch sampler over 2 analyses (use $sampler with cm_monte_carlo())
```

The result holds:

- `population`: the diseases, associations and three-way terms, as a
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  object;
- `analyses`: a
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
  object with one model per impact table;
- `model`: the
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  when there is exactly one analysis (otherwise `NULL`);
- `hr_model`: the hazard-ratio model, for
  [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md);
- `sampler`: a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
  for one analysis, or a
  [`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md)
  for several, for
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md).

``` r

inp$population
#> <cm_population>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
names(inp$analyses$models)
#> [1] "calving_interval" "yield"
deconflate(inp$analyses)
#> <cm_results> 2 analyses: calving_interval, yield
#> 
#> == calving_interval ==
#> <cm_result> method: simultaneous; calving interval increase [days]
#> 
#>  disease raw adjusted  change
#>      LAM  12    8.922 -0.2565
#>      SCK   4    1.935 -0.5162
#>      MET  18   14.044 -0.2198
#> 
#> Raw sum: 6.325; adjusted total: 4.373
#> Diagnostics: residual 0.00e+00, condition number 2, sign changes 0
#> 
#> == yield ==
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0
deconflate_hr(inp$hr_model)
#> <cm_hr_result> method: snapshot
#> 
#>  disease  raw adjusted   change estimand
#>      LAM 1.74    1.532 -0.11980    crude
#>      SCK 1.92    1.786 -0.06976    crude
#>      MET 1.50    1.141 -0.23947    crude
#> 
#> Diagnostics:
#>  max_reconstruction_residual n_sign_changes condition_number
#>                      1.1e-13              0              2.1
#>                                       feasibility
#>  joint distribution fitted (max residual 6.6e-11)
```

The template has two analyses, so the sampler is a batch sampler: every
draw of the shared disease and association inputs is used by both
analyses.

``` r

inp$sampler
#> <cm_batch_sampler> 2 analyses (calving_interval, yield); 2 shared population inputs
mc <- cm_monte_carlo(inp$sampler, 200, seed = 1)
mc
#> <cm_mc_batch> 200 draws, 2 analyses
#>   calving_interval: 200 accepted, 0 rejected
#>   yield: 200 accepted, 0 rejected
s <- summary(mc, diagnose = FALSE)
s[, c("analysis", "disease", "mean", "q0.025", "q0.975")]
#>           analysis disease      mean    q0.025    q0.975
#> 1 calving_interval     LAM  8.933671 7.0997909 10.666460
#> 2 calving_interval     SCK  1.942984 1.4555402  2.520507
#> 3 calving_interval     MET 14.154484 4.2837150 22.834133
#> 4            yield     LAM  3.000444 0.7481448  4.965309
#> 5            yield     SCK  7.908887 5.7769745 10.105338
#> 6            yield     MET  3.099691 2.2065559  3.869254
```

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
when there are errors, and prints notes otherwise. The installed example
with errors shows more cases:

``` r

cm_check_inputs(dir = system.file("extdata", "example_with_errors", package = "deconflate"))
#> Found 23 problem(s) in the inputs:
#>   interactions_culling: There is no impact table for analysis 'culling'.
#>   diseases, row 5, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   associations, row 3, column 'disease2': Unknown disease 'CM' (not in the diseases table).
#>   associations, row 4: Duplicate pair SCK:LAM.
#>   associations, row 5, column 'value': Missing value for measure OR.
#>   associations, row 6, column 'adjusted': A covariate-adjusted measure is not a marginal 2x2 association. Use the crude measure, or set adjusted_associations = 'use_as_marginal' to use it as an approximation.
#>   hazard_ratios, row 3, column 'value': Hazard ratios must be positive.
#>   impacts_fertility, column 'units': One units per table is allowed (found: % increase, days).
#>   impacts_yield, row 2, column 'value': '8,40' is not a number.
#>   impacts_yield, row 3, column 'adjusted_for': adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.
#>   impacts_yield: No impact for: MET. Add a row with value 0 for no impact.
#>   uncertainty, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   uncertainty, row 3, column 'key': Key 'assoc:LAM:RP' does not match an association with a numeric measure (expected assoc:<d1>:<d2>).
#>   uncertainty, row 4, column 'key': Key 'impact:yield:CM': unknown disease 'CM'.
#>   uncertainty, row 5, column 'key': Key 'impact:culling:LAM': unknown analysis 'culling'.
#>   uncertainty, row 6, column 'key': Key 'impact:LAM': name the analysis (impact:<analysis>:<disease>), since there are several.
#>   uncertainty, row 7, column 'dist': PERT needs min <= mode <= max and min < max.
#>   uncertainty, row 9, column 'key': Key 'inter:yield:LAM:SCK' does not match an interaction of analysis 'yield'.
#>   uncertainty, row 10, column 'key': Key 'milk_price' must start with prob:, assoc:, impact: or inter:.
#>   uncertainty, row 11, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
#>   uncertainty_fertility, row 3, column 'key': Disease probabilities are shared by all analyses; put prob: keys in uncertainty.csv.
#>   uncertainty_fertility, row 2, column 'key': The same input is given more than once.
```

## The global dairy example

The folder `global_dairy_2024` holds the inputs of Rasmussen et
al. (2024): `diseases.csv` (12 diseases), `associations.csv` (38 odds
ratios), `impacts_yield.csv` and `impacts_fertility.csv` (in percent),
`hazard_ratios.csv` (culling) and `uncertainty.csv` (the input
distributions).

``` r

gd <- cm_read_inputs(dir = system.file("extdata", "global_dairy_2024", package = "deconflate"))
gd
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 12 (CK, CM, DA, DYS, LAM, MET, MF, OC, PTB, RP, SCK, SCM)
#>   Disease pairs: 66 [independent (default): 28; specified: 38]
#>   Analyses: fertility, yield
#>   Hazard ratios: yes (use $hr_model with deconflate_hr())
#>   Uncertain inputs: batch sampler over 2 analyses (use $sampler with cm_monte_carlo())
cmp <- compare_methods(gd$analyses, methods = c("published", "simultaneous"))
cmp$impacts[cmp$impacts$analysis == "yield", ]
#>    analysis disease       raw  published simultaneous
#> 13    yield      CK 0.4321944 0.02473297   -5.3709141
#> 14    yield      CM 3.2499000 1.33194993   -0.3380725
#> 15    yield      DA 2.8369300 0.79378291   -2.6316316
#> 16    yield     DYS 4.9190880 3.56881898    4.3096749
#> 17    yield     LAM 4.8061000 2.53240737    1.9792779
#> 18    yield     MET 5.6130850 2.84095146    2.7909168
#> 19    yield      MF 0.5365131 0.06897532   -1.5399737
#> 20    yield      OC 3.7478390 2.63859343    3.2174755
#> 21    yield     PTB 4.3000000 3.22917066    3.9410642
#> 22    yield      RP 4.1986640 2.25771112    2.4733053
#> 23    yield     SCK 8.3964720 7.10307185    8.2785938
#> 24    yield     SCM 6.2931840 5.58961905    6.5761493
```

The simultaneous solution changes the sign of several adjusted impacts
(see
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).
The culling hazard ratios are adjusted separately
([`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md)):

``` r

deconflate_hr(gd$hr_model, method = "published")
#> <cm_hr_result> method: published
#> 
#>  disease   raw adjusted   change estimand
#>       CK 1.500    1.178 -0.21500    crude
#>       CM 2.300    1.904 -0.17220    crude
#>       DA 2.851    2.198 -0.22912    crude
#>      DYS 1.258    1.098 -0.12699    crude
#>      LAM 1.745    1.381 -0.20877    crude
#>      MET 1.116    1.012 -0.09318    crude
#>       MF 3.000    2.648 -0.11742    crude
#>       OC 1.620    1.459 -0.09960    crude
#>      PTB 2.311    2.047 -0.11395    crude
#>       RP 1.600    1.284 -0.19715    crude
#>      SCK 1.920    1.675 -0.12747    crude
#>      SCM 1.450    1.255 -0.13453    crude
#> 
#> Diagnostics:
#>  max_reconstruction_residual n_sign_changes condition_number
#>                         0.49              0             3.25
#>                                      feasibility
#>  triple screen passed (necessary condition only)
```

## Tables typed in R

Any table can be a data frame instead of a file. Tables that come one
per analysis (`impacts`, `interactions`, `uncertainty`) are named lists,
named after the analyses; `uncertainty` can also have an element named
`shared` for the population inputs:

``` r

inp2 <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease1 = c("d1", "d2"), disease2 = c("d2", "d3"),
                            value = c(2, 3), measure = "OR"),
  impacts = list(
    yield = data.frame(disease = c("d1", "d2", "d3"), value = c(2.5, 5, 7.5),
                       units = "% of yield"),
    fertility = data.frame(disease = c("d1", "d2", "d3"), value = c(1, 2, 0),
                           units = "% of calving interval")
  ),
  uncertainty = list(
    shared = data.frame(key = "assoc:d1:d2", dist = "lognormal_ci", p1 = 2, p2 = 1.4, p3 = 2.9),
    yield = data.frame(key = "impact:d1", dist = "normal", p1 = 2.5, p2 = 0.5)
  )
)
inp2
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
#>   Analyses: yield, fertility
#>   Uncertain inputs: batch sampler over 2 analyses (use $sampler with cm_monte_carlo())
deconflate(inp2$analyses)$yield$adjusted
#>   disease raw adjusted      change estimand adjusted_for
#> 1      d1 2.5 2.143250 -0.14269981    crude         <NA>
#> 2      d2 5.0 3.387080 -0.32258408    crude         <NA>
#> 3      d3 7.5 6.934209 -0.07543874    crude         <NA>
```

With a single analysis (e.g. `impacts = list(yield = ...)`, or one data
frame, which becomes the analysis `impacts`), `$model` is set and
`$sampler` is a single-analysis sampler.

## Associations from other designs

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
#>   Disease pairs: 1 [specified: 1]
#>   Covariate-adjusted association measures used as marginal: 1
```

In files, set the column `adjusted` to `TRUE` and call
`cm_read_inputs(..., adjusted_associations = "use_as_marginal")`; each
such row is then reported as a note.

### Contingency tables and zero cells

A study table can be given as counts (`measure = "table"`). Its odds
ratio is applied to the modelled population’s own probabilities. A table
with a zero cell has no finite odds ratio, so by default 0.5 is added to
every cell (Haldane-Anscombe correction), and the row is flagged in
column `corrected`:

``` r

tab <- cm_associations("d1", "d2", measure = "table", n11 = 0, n10 = 30, n01 = 40, n00 = 930)
tab[, c("disease1", "disease2", "measure", "value", "corrected")]
#>   disease1 disease2 measure     value corrected
#> 1       d1       d2   table 0.3766444      TRUE
tryCatch(cm_associations("d1", "d2", measure = "table", n11 = 0, n10 = 30, n01 = 40,
                         n00 = 930, zero_cell = "error"),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "A contingency table has a zero cell; use zero_cell = 'haldane' to add 0.5 to every cell."
```

When reading files, the correction is applied and reported as a note:

``` r

cm_check_inputs(
  diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
  associations = data.frame(disease1 = "d1", disease2 = "d2", measure = "table",
                            n11 = 0, n10 = 30, n01 = 40, n00 = 930)
)
#> Notes on the inputs:
#>   associations, row 1 (note): The table has a zero cell; 0.5 was added to every cell (Haldane correction).
```

Tables in which every count is zero are rejected.

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

## Uncertainty

Each row of an uncertainty table gives one input a distribution. The key
names the input:

- `prob:<disease>`: the disease’s `value` (e.g. a prevalence or
  incidence rate);
- `assoc:<d1>:<d2>`: the association’s `value` (e.g. an odds ratio);
- `impact:<analysis>:<disease>`: an impact, in the units of its table;
- `inter:<analysis>:<d1>:<d2>`: an interaction.

In `uncertainty_<analysis>.csv`, or when there is only one analysis, the
analysis can be left out (`impact:<disease>`, `inter:<d1>:<d2>`).
Disease and association keys belong in the shared `uncertainty.csv`.

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

``` r

read.csv(file.path(dir, "uncertainty.csv"))
#>                           key   dist    p1     p2 p3 p4
#> 1                    prob:LAM   beta 78.29 227.42 NA NA
#> 2               assoc:LAM:SCK normal  2.01   0.20  0 NA
#> 3            impact:yield:LAM normal  4.81   0.87 NA NA
#> 4            impact:yield:SCK normal  8.40   1.19 NA NA
#> 5 impact:calving_interval:MET   pert  6.00  18.00 30 NA
#>                              note
#> 1            beta(shape1, shape2)
#> 2 normal(mean, sd) truncated at 0
#> 3       normal(mean, sd), percent
#> 4       normal(mean, sd), percent
#> 5      pert(min, mode, max), days
```

[`cm_dist_table()`](https://rasmussenphilip.github.io/deconflate/reference/cm_dist_table.md)
builds the distributions alone, if you want to pass them to
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
yourself. See
[`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md)
for the Monte Carlo analysis.
