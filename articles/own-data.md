# Using your own data

``` r

library(deconflate)
```

Inputs can be kept in CSV files (one per table) or typed in R as data
frames.
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
checks every table, reports all problems at once, and returns the
population, the model (one analysis per impact table), a hazard-ratio
model if there are hazard ratios, and a Monte Carlo sampler if any value
has a distribution.

## Files

| File | Contents | Required columns |
|----|----|----|
| `diseases.csv` | one row per disease | `id`, `value` |
| `associations.csv` | one row per associated pair | `disease1`, `disease2` |
| `three_way.csv` | optional three-way association scenarios | `disease1`, `disease2`, `disease3`, `ratio` |
| `impacts.csv` | the impact table of a single analysis | `disease`, `value` |
| `impacts_<analysis>.csv` | for several analyses, one impact table each, e.g. `impacts_yield.csv` | `disease`, `value` |
| `interactions.csv`, `interactions_<analysis>.csv` | optional impact interactions | `disease1`, `disease2`, `value` |
| `hazard_ratios.csv` | optional culling (or mortality) hazard ratios | `disease`, `value`, `estimand` |

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
- impact tables: `estimand` (`crude` or `adjusted_linear`),
  `adjusted_for`, `source`, and one `label` and one `units` per table;
- `hazard_ratios.csv`: `adjusted_for`, `source`. The `estimand` column
  is required (`snapshot_crude` or `snapshot_stratified`; see below).

Every disease needs a row in every impact table (use 0 for no impact)
and in the hazard-ratio table (use 1 for no effect). Impact values can
be in any units, the same within a table; results come back in those
units.

There is no separate uncertainty file: an uncertain value gets its
distribution in its own row, in the columns `dist` and `p1`-`p4` (see
Uncertainty below), so point values and uncertain values sit side by
side. Every table can also have a free-text `note` column.

## Start from the template

[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set of files to a folder. By default it is a single
analysis, the standard workflow:

``` r

dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 5 files to /tmp/RtmpuauJId/my-inputs
list.files(dir)
#> [1] "associations.csv" "diseases.csv"     "impacts.csv"      "interactions.csv"
#> [5] "three_way.csv"
read.csv(file.path(dir, "impacts.csv"))
#>   disease value estimand adjusted_for           label      units
#> 1     LAM  4.81    crude           NA milk yield loss % of yield
#> 2     SCK  8.40    crude           NA milk yield loss % of yield
#> 3     MET  5.61    crude           NA milk yield loss % of yield
#>                source   dist   p1   p2 p3 p4             note
#> 1 Illustrative values normal 4.81 0.87 NA NA normal(mean, sd)
#> 2 Illustrative values normal 8.40 1.19 NA NA normal(mean, sd)
#> 3 Illustrative values          NA   NA NA NA
```

`interactions.csv` and `three_way.csv` are empty (header only); fill
them in if needed.

## Read the files

``` r

inp <- cm_read_inputs(dir = dir)
inp
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
#>   Analyses: impacts
#>   Uncertain inputs: 4 (use $sampler with cm_monte_carlo())
```

The result holds:

- `population`: the diseases, associations and three-way terms, as a
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  object;
- `model`: the
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  when there is exactly one analysis (otherwise `NULL`);
- `analyses`: a
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
  object with one model per impact table (here one, named `impacts`);
- `hr_model`: the hazard-ratio model, for
  [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
  (`NULL` without `hazard_ratios.csv`);
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
deconflate(inp$model)
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

## From the files to a Monte Carlo run

In the template, four values have a distribution in their row: the
lameness prevalence (`diseases.csv`), the lameness:subclinical ketosis
odds ratio (`associations.csv`) and two yield impacts (`impacts.csv`).
The other values are point values (empty `dist`). With one analysis, the
sampler is a
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md):

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
inp$sampler
#> <cm_sampler> 4 uncertain inputs
#>   prob:LAM: beta
#>   assoc:LAM:SCK: normal
#>   impact:LAM: normal
#>   impact:SCK: normal
mc <- cm_monte_carlo(inp$sampler, 200, seed = 1)
mc
#> <cm_mc> method: simultaneous
#>   Analysis: milk yield loss [% of yield]
#>   Draws: 200, rejected: 0 (0.0%)
#>   Sampling: random
#>   Effective sample size: 200.0
summary(mc, diagnose = FALSE)[, c("disease", "mean", "mcse", "q0.025", "q0.975")]
#>   disease     mean       mcse   q0.025   q0.975
#> 1     LAM 2.941653 0.06803879 1.109592 4.926224
#> 2     SCK 7.793789 0.08327428 5.743658 9.896458
#> 3     MET 3.142470 0.02738527 2.336560 3.849221
```

A single-analysis sampler supports Latin hypercube sampling
(`sampling = "lhs"`) and importance sampling. For importance sampling,
each proposal must cover the whole support of the input’s own
distribution. A defensive mixture, part the input’s own distribution and
part a wider or shifted one, does so by construction. Here the lameness
yield impact is sampled more often around 6%:

``` r

specs <- attr(inp$sampler, "specs")
prop <- list("impact:LAM" = dist_mixture(specs[["impact:LAM"]], dist_normal(6, 1.5),
                                         weights = c(0.5, 0.5)))
mc_is <- cm_monte_carlo(inp$sampler, 200, proposal = prop, seed = 1)
mc_is
#> <cm_mc> method: simultaneous
#>   Analysis: milk yield loss [% of yield]
#>   Draws: 200, rejected: 0 (0.0%)
#>   Sampling: random, importance sampling of impact:LAM
#>   Effective sample size: 167.9
summary(mc_is, diagnose = FALSE)[, c("disease", "mean", "mcse", "q0.025", "q0.975")]
#>   disease     mean       mcse    q0.025    q0.975
#> 1     LAM 2.780003 0.07196501 0.7373599  4.709681
#> 2     SCK 7.966686 0.09387597 5.7170650 10.075637
#> 3     MET 3.180425 0.02769793 2.4727201  3.881797
```

The draws are weighted by the ratio of the densities, so the estimates
still refer to the input’s own distribution; compare the standard errors
(`mcse`) to see whether the proposal helped. See
[`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md).

## Several analyses

`cm_template(type = "analyses")` writes two analyses on one population,
milk yield (in percent) and calving interval (in days):

``` r

dir2 <- file.path(tempdir(), "my-analyses")
cm_template(dir2, type = "analyses", overwrite = TRUE)
#> Wrote 6 files to /tmp/RtmpuauJId/my-analyses
list.files(dir2)
#> [1] "associations.csv"             "diseases.csv"                
#> [3] "impacts_calving_interval.csv" "impacts_yield.csv"           
#> [5] "interactions_yield.csv"       "three_way.csv"
read.csv(file.path(dir2, "impacts_calving_interval.csv"))
#>   disease value estimand adjusted_for                     label units
#> 1     LAM    12    crude           NA calving interval increase  days
#> 2     SCK     4    crude           NA calving interval increase  days
#> 3     MET    18    crude           NA calving interval increase  days
#>                source dist p1 p2 p3 p4                       note
#> 1 Illustrative values      NA NA NA NA                           
#> 2 Illustrative values      NA NA NA NA                           
#> 3 Illustrative values pert  6 18 30 NA pert(min, mode, max), days
inp2 <- cm_read_inputs(dir = dir2)
inp2
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
#>   Analyses: calving_interval, yield
#>   Uncertain inputs: batch sampler over 2 analyses (use $sampler with cm_monte_carlo())
names(inp2$analyses$models)
#> [1] "calving_interval" "yield"
deconflate(inp2$analyses)
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
```

With several analyses, `model` is `NULL` and the sampler is a batch
sampler: every draw of the shared disease and association inputs is used
by both analyses. Batch runs use simple random sampling.

``` r

inp2$sampler
#> <cm_batch_sampler> 2 analyses (calving_interval, yield); 2 shared population inputs
mcb <- cm_monte_carlo(inp2$sampler, 200, seed = 1)
mcb
#> <cm_mc_batch> 200 draws, 2 analyses
#>   calving_interval: 200 accepted, 0 rejected
#>   yield: 200 accepted, 0 rejected
s <- summary(mcb, diagnose = FALSE)
s[, c("analysis", "disease", "mean", "q0.025", "q0.975")]
#>           analysis disease      mean    q0.025    q0.975
#> 1 calving_interval     LAM  8.933671 7.0997909 10.666460
#> 2 calving_interval     SCK  1.942984 1.4555402  2.520507
#> 3 calving_interval     MET 14.154484 4.2837150 22.834133
#> 4            yield     LAM  3.000444 0.7481448  4.965309
#> 5            yield     SCK  7.908887 5.7769745 10.105338
#> 6            yield     MET  3.099691 2.2065559  3.869254
```

## Hazard ratios

Hazard ratios go in `hazard_ratios.csv`. They are adjusted by a separate
model
([`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md)),
and the `estimand` column states how each value is to be read:
`snapshot_crude`, or `snapshot_stratified` with `adjusted_for`. Entering
a published Cox hazard ratio as either is an approximation and your
assumption, so there is no default:

``` r

write.csv(data.frame(disease = c("LAM", "SCK", "MET"), value = c(1.74, 1.5, 1.9),
                     estimand = "snapshot_crude", source = "Illustrative"),
          file.path(dir, "hazard_ratios.csv"), row.names = FALSE)
inp_hr <- cm_read_inputs(dir = dir)
deconflate_hr(inp_hr$hr_model)
#> <cm_hr_result> method: snapshot
#> 
#>  disease  raw adjusted   change       estimand
#>      LAM 1.74    1.496 -0.13995 snapshot_crude
#>      SCK 1.50    1.366 -0.08945 snapshot_crude
#>      MET 1.90    1.524 -0.19768 snapshot_crude
#> 
#> Diagnostics:
#>  max_reconstruction_residual n_sign_changes condition_number
#>                     2.69e-13              0             2.21
#>                                       feasibility
#>  joint distribution fitted (max residual 6.6e-11)
```

A table without the `estimand` column, or with the names used before
version 0.3.0 (`crude`, `adjusted`), is reported:

``` r

d3 <- data.frame(id = c("LAM", "SCK", "MET"), value = c(0.25, 0.38, 0.10))
cm_check_inputs(diseases = d3,
                hazard_ratios = data.frame(disease = c("LAM", "SCK", "MET"),
                                           value = c(1.74, 1.5, 1.9)))
#> Found 1 problem(s) in the inputs:
#>   hazard_ratios, column 'estimand': Required column is missing.
cm_check_inputs(diseases = d3,
                hazard_ratios = data.frame(disease = c("LAM", "SCK", "MET"),
                                           value = c(1.74, 1.5, 1.9), estimand = "crude"))
#> Found 3 problem(s) in the inputs:
#>   hazard_ratios, row 1, column 'estimand': Unknown estimand 'crude' (use snapshot_crude or snapshot_stratified; see ?cm_hazard_ratios).
#>   hazard_ratios, row 2, column 'estimand': Unknown estimand 'crude' (use snapshot_crude or snapshot_stratified; see ?cm_hazard_ratios).
#>   hazard_ratios, row 3, column 'estimand': Unknown estimand 'crude' (use snapshot_crude or snapshot_stratified; see ?cm_hazard_ratios).
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
#> Found 19 problem(s) in the inputs:
#>   uncertainty: Uncertainty files are no longer read: give each uncertain value a distribution in its own table, in the columns dist and p1-p4 (see ?cm_read_inputs).
#>   interactions_culling: There is no impact table for analysis 'culling'.
#>   diseases, row 5, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   diseases, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   diseases, row 4, column 'dist': Parameters p1-p4 are given without a distribution (dist).
#>   associations, row 3, column 'disease2': Unknown disease 'CM' (not in the diseases table).
#>   associations, row 4: Duplicate pair SCK:LAM.
#>   associations, row 5, column 'value': Missing value for measure OR.
#>   associations, row 6, column 'adjusted': A covariate-adjusted measure is not a marginal 2x2 association. Use the crude measure, or set adjusted_associations = 'use_as_marginal' to use it as an approximation.
#>   associations, row 2, column 'dist': PERT needs min <= mode <= max and min < max.
#>   associations, row 7, column 'dist': A distribution needs a numeric measure (OR, RR, RD, cond_prob or phi); give the measure and its point value.
#>   hazard_ratios, column 'dist' (note): Uncertainty is not used for hazard ratios (Monte Carlo runs adjust additive impacts); the dist and p1-p4 columns are ignored.
#>   hazard_ratios, row 3, column 'value': Hazard ratios must be positive.
#>   impacts_fertility, column 'units': One units per table is allowed (found: % increase, days).
#>   impacts_fertility, row 3, column 'p2': '8,5' is not a number.
#>   impacts_fertility, row 2, column 'dist' (note): The point value 1.12 lies outside the distribution's support [1.5, 3].
#>   impacts_yield, row 2, column 'value': '8,40' is not a number.
#>   impacts_yield, row 3, column 'adjusted_for': adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.
#>   impacts_yield: No impact for: MET. Add a row with value 0 for no impact.
#>   impacts_yield, row 3, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
```

## The global dairy example

The folder `global_dairy_2024` holds the inputs of Rasmussen et
al. (2024): `diseases.csv` (12 diseases), `associations.csv` (38 odds
ratios), `impacts_yield.csv` and `impacts_fertility.csv` (in percent)
and `hazard_ratios.csv` (culling, entered as `snapshot_crude`). The
input distributions of the paper’s Monte Carlo analysis are in the
`dist` and `p1`-`p4` columns of the association and impact tables.

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
#>  disease   raw adjusted   change       estimand
#>       CK 1.500    1.178 -0.21500 snapshot_crude
#>       CM 2.300    1.904 -0.17220 snapshot_crude
#>       DA 2.851    2.198 -0.22912 snapshot_crude
#>      DYS 1.258    1.098 -0.12699 snapshot_crude
#>      LAM 1.745    1.381 -0.20877 snapshot_crude
#>      MET 1.116    1.012 -0.09318 snapshot_crude
#>       MF 3.000    2.648 -0.11742 snapshot_crude
#>       OC 1.620    1.459 -0.09960 snapshot_crude
#>      PTB 2.311    2.047 -0.11395 snapshot_crude
#>       RP 1.600    1.284 -0.19715 snapshot_crude
#>      SCK 1.920    1.675 -0.12747 snapshot_crude
#>      SCM 1.450    1.255 -0.13453 snapshot_crude
#> 
#> Diagnostics:
#>  max_reconstruction_residual n_sign_changes condition_number
#>                         0.49              0             3.25
#>                                      feasibility
#>  triple screen passed (necessary condition only)
```

## Tables typed in R

Any table can be a data frame instead of a file. Tables that come one
per analysis (`impacts`, `interactions`) are named lists, named after
the analyses. Distributions are columns of the tables, as in the files:

``` r

inp3 <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease1 = c("d1", "d2"), disease2 = c("d2", "d3"),
                            value = c(2, 3), measure = "OR",
                            dist = c("lognormal_ci", NA), p1 = c(2, NA), p2 = c(1.4, NA),
                            p3 = c(2.9, NA)),
  impacts = list(
    yield = data.frame(disease = c("d1", "d2", "d3"), value = c(2.5, 5, 7.5),
                       units = "% of yield", dist = c("normal", NA, NA),
                       p1 = c(2.5, NA, NA), p2 = c(0.5, NA, NA)),
    fertility = data.frame(disease = c("d1", "d2", "d3"), value = c(1, 2, 0),
                           units = "% of calving interval")
  )
)
inp3
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
#>   Analyses: yield, fertility
#>   Uncertain inputs: batch sampler over 2 analyses (use $sampler with cm_monte_carlo())
inp3$sampler
#> <cm_batch_sampler> 2 analyses (yield, fertility); 1 shared population inputs
deconflate(inp3$analyses)$yield$adjusted
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

The tables `diseases`, `associations`, `three_way`, `impacts` and
`interactions` can have the columns `dist` and `p1`-`p4`. A row with a
`dist` has an uncertain value; a row with an empty `dist` is a point
value. The deterministic methods use the point value (`value`, or
`ratio` for a three-way term); Monte Carlo runs draw from the
distribution, on the scale the value was entered on:

- a disease’s `value` (e.g. a prevalence or an incidence rate);
- an association’s `value` (e.g. an odds ratio); the measure must be
  numeric (`OR`, `RR`, `RD`, `cond_prob` or `phi`), not `table`,
  `independent` or `unknown`;
- a three-way `ratio`;
- an impact or an interaction, in the units of its table.

Hazard ratios cannot have distributions (the Monte Carlo tools adjust
additive impacts). `dist` and the parameters `p1`-`p4`:

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

The checks report a distribution with missing or invalid parameters,
parameters without a `dist`, and a distribution on a non-numeric
association measure; a point value outside its distribution’s range is
noted. In the sampler, inputs are keyed `prob:<disease>`,
`assoc:<d1>:<d2>`, `three:<d1>:<d2>:<d3>`, `impact:<disease>` and
`inter:<d1>:<d2>` (the keys of `mc$params`, of proposals in importance
sampling and of
[`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)):

``` r

names(attr(inp2$sampler$samplers$yield, "specs"))
#> [1] "prob:LAM"      "assoc:LAM:SCK" "impact:LAM"    "impact:SCK"
```

Files named `uncertainty.csv` or `uncertainty_<analysis>.csv`, used by
version 0.2, are reported as an error: move each distribution into the
row of its value.

[`cm_dist_table()`](https://rasmussenphilip.github.io/deconflate/reference/cm_dist_table.md)
builds distributions from a table with columns `key`, `dist` and
`p1`-`p4`, if you want to pass them to
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
yourself. See
[`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md)
for the Monte Carlo analysis.

## A fuller example

The folder `five_diseases` uses every feature on five diseases: all
probability types and association measures, a three-way term, crude and
adjusted impacts, interactions, hazard ratios, and every distribution
type mixed with point values. Its `run_all_features.R` script runs every
part of the package on it:

``` r

five <- system.file("extdata", "five_diseases", package = "deconflate")
list.files(five)
#>  [1] "associations.csv"             "diseases.csv"                
#>  [3] "hazard_ratios.csv"            "impacts_calving_interval.csv"
#>  [5] "impacts_welfare.csv"          "impacts_yield.csv"           
#>  [7] "interactions_welfare.csv"     "README.md"                   
#>  [9] "run_all_features.R"           "three_way.csv"
inp5 <- cm_read_inputs(dir = five)
inp5
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 5 (LAM, MAS, MET, SCK, RP)
#>   Disease pairs: 10 [specified: 10]
#>   Three-way terms: 1
#>   Analyses: calving_interval, welfare, yield
#>   Hazard ratios: yes (use $hr_model with deconflate_hr())
#>   Uncertain inputs: batch sampler over 3 analyses (use $sampler with cm_monte_carlo())
read.csv(file.path(five, "impacts_yield.csv"))[, c("disease", "value", "dist", "p1", "p2", "p3")]
#>   disease value      dist  p1  p2 p3
#> 1     LAM   4.8    normal 4.8 0.9 NA
#> 2     MAS   3.3      pert 1.5 3.3  6
#> 3     MET   5.6 pert_mean 3.0 5.6  8
#> 4     SCK   2.5     fixed 2.5  NA NA
#> 5      RP   4.2            NA  NA NA
```
