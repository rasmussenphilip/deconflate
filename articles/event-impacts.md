# Event impacts (culling and mortality)

``` r

library(deconflate)
```

Some impacts are not amounts but comparisons of the risk of an event,
such as culling or death, between animals with and without a disease.
They are reported as hazard ratios (HR) from survival models, rate
ratios, risk ratios (RR), odds ratios (OR) or risk differences (RD).
Such an estimate for one disease is conflated with the effects of
associated diseases, as an additive impact is, but the effects combine
multiplicatively, not additively. These are **event impacts**, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
adjusts them with a separate model, the snapshot hazard model, when
called with `event_model = TRUE`.

All five measures are mapped onto the same model, so they can be mixed
in one impact table (one row per disease): an HR for one disease, a risk
difference for another. The result is the same for every measure: an
adjusted hazard ratio per disease (its own hazard multiplier), and the
risk of the event attributable to disease, allocated to diseases.

## Describing event impacts

An impact table holds event impacts when it has a `measure` column:
`"HR"`, `"rate_ratio"`, `"RR"`, `"OR"` or `"RD"`. Every disease needs a
row (1 for no effect on a ratio scale, 0 for a risk difference). The
estimand must be stated; there is no default, because it is an
assumption about what the source estimate measures (see below):

``` r

pop <- example_supplement()   # its yield impacts are not used here
cull <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), measure = "HR",
                                 estimand = "snapshot_crude", label = "culling"))
cull
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (3 with an association, 0 unknown)
#>   Impacts: culling (event impacts: use event_model = TRUE)
#>   Estimands: snapshot_crude: 3
#>   Measures: HR: 3
tryCatch(cm_impacts("d1", 1.5, measure = "HR"),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "State the estimand of event impacts: estimand = \"snapshot_crude\" or \"snapshot_stratified\" (see ?cm_impacts)."
```

In a CSV file read by `cm_read_inputs(impacts = )`, the same information
is in the columns `disease`, `value`, `measure`, `estimand` and
`adjusted_for` (see the five-disease example at the end).

[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
needs `event_model = TRUE` and the overall risk of the event:

``` r

res <- deconflate(cull, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
res
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure raw adjusted_hr
#>       d1      HR 1.5       1.383
#>       d2      HR 2.0       1.891
#>       d3      HR 1.3       1.144
#> 
#> Overall risk 0.25; disease-free risk 0.2134; attributable to disease 0.03665 (14.7% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable  share
#>       d1     0.007371 0.2011
#>       d2     0.023488 0.6409
#>       d3     0.005791 0.1580
#> 
#> Diagnostics: residual 2.22e-16, condition number 1.68, sign changes 0
```

The kind of impacts and the call must agree. Event impacts without
`event_model = TRUE`, additive impacts with it, or a missing
`overall_risk` are errors:

``` r

tryCatch(deconflate(cull), deconflate_error = function(e) conditionMessage(e))
#> [1] "These are event impacts (the impact table has a measure column): use event_model = TRUE, with overall_risk."
tryCatch(deconflate(pop, event_model = TRUE, overall_risk = 0.25),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "event_model = TRUE needs event impacts: an impact table with a measure column (HR, rate_ratio, RR, OR or RD) and snapshot estimands (see ?cm_impacts)."
tryCatch(deconflate(cull, event_model = TRUE),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "event_model = TRUE needs `overall_risk`: the overall risk of the event in the population over the period of the estimates (a proportion, e.g. 0.25, or a distribution such as dist_beta(250, 750))."
```

## The snapshot model

Within the period of the overall risk, an animal with disease
combination `d` has a constant hazard `h0 * exp(sum_i beta_i * d_i)`, so
its risk of the event over the period is
`R(d) = 1 - exp(-h0 * exp(sum_i beta_i * d_i))`. The baseline `h0` is
set so that the population risk, averaged over the distribution of
disease combinations
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)),
equals `overall_risk`. The `beta`s are solved so that the model
reproduces each raw estimate:

- a hazard ratio or rate ratio by the ratio of the average hazard
  multipliers among animals with and without the disease;
- a risk ratio, odds ratio or risk difference by the ratio, odds ratio
  or difference of the average risks `R` among animals with and without
  the disease.

The adjusted hazard ratios are `exp(beta)`. The comparison is made at
the start of the period (the “snapshot”); its estimands are defined by
the model:

- `"snapshot_crude"`: the comparison of animals with and without the
  disease at the start of the period, in the population described by the
  probabilities and associations;
- `"snapshot_stratified"`: the same comparison within strata of the
  diseases in `adjusted_for`, combined across strata (see below).

Because the model uses the whole distribution of disease combinations,
it depends on the maximum-entropy assumption of
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
and on any three-way terms
([`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)).
Pairs without an association are unknown, and the fit fills them in, as
in the global method for additive impacts; the result lists them with
their fitted odds ratios (`$unknown_pairs`).

### What the snapshot model is not

A Cox HR estimated over follow-up is not, in general, the snapshot
ratio: animals with high hazards leave first, so the mixture of disease
combinations among survivors changes over time, and the marginal HR
changes with it. The snapshot model is exact for its own estimands.
Entering a published HR as `"snapshot_crude"` is an approximation that
you state, not a conversion the package makes. It is reasonable when
follow-up is short relative to the hazards or the diseases are rare, and
it should be reported and checked.

## The overall risk

`overall_risk` is the proportion of animals with the event over the
period (e.g. the annual culling rate). It is required. For hazard and
rate ratios it fixes the baseline hazard, and so the attributable risk,
but not the adjusted hazard ratios:

``` r

risks <- c(0.10, 0.25, 0.40)
sapply(risks, function(r) {
  deconflate(cull, event_model = TRUE, overall_risk = r, n_draws = 0)$adjusted$adjusted
})
#>          [,1]     [,2]     [,3]
#> [1,] 1.383153 1.383153 1.383153
#> [2,] 1.890932 1.890932 1.890932
#> [3,] 1.143987 1.143987 1.143987
```

Risk ratios, odds ratios and risk differences compare risks over a
period, so they depend on the risk level, and the adjusted hazard ratios
change with `overall_risk`. They must refer to the same period as
`overall_risk` (a risk ratio of culling within the lactation needs the
culling risk over the lactation):

``` r

cull_rr <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(1.4, 1.8, 1.2), measure = "RR",
                                    estimand = "snapshot_crude"))
sapply(risks, function(r) {
  deconflate(cull_rr, event_model = TRUE, overall_risk = r, n_draws = 0)$adjusted$adjusted
})
#>          [,1]     [,2]     [,3]
#> [1,] 1.332659 1.386210 1.480956
#> [2,] 1.800381 1.948599 2.225537
#> [3,] 1.082457 1.093058 1.110192
```

Not every set of estimates fits every overall risk. Risk ratios of 4 for
every disease cannot hold when half of the animals have the event:
animals with d2 would need a risk above 1.

``` r

cull_rr4 <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(4, 4, 4), measure = "RR",
                                     estimand = "snapshot_crude"))
tryCatch(deconflate(cull_rr4, event_model = TRUE, overall_risk = 0.5, n_draws = 0),
         deconflate_error = function(e) conditionMessage(e))
#> [1] "Could not solve the snapshot model (max residual 5.75e-01): the estimates may be incompatible with the overall risk or with each other."
```

The overall risk can be uncertain: give it as a distribution, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
draws it with the other uncertain inputs
([`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md)).
Here only the overall risk is uncertain; the adjusted hazard ratios do
not depend on it, so their intervals have zero width, while the
attributable risk gets an interval. The central estimate uses the mean
of the distribution:

``` r

deconflate(cull, event_model = TRUE, overall_risk = dist_beta(250, 750),
           n_draws = 200, seed = 1)
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure raw adjusted_hr lower upper
#>       d1      HR 1.5       1.383 1.383 1.383
#>       d2      HR 2.0       1.891 1.891 1.891
#>       d3      HR 1.3       1.144 1.144 1.144
#> 
#> Overall risk 0.25; disease-free risk 0.2134; attributable to disease 0.03665 (14.7% of the overall risk); 95% interval [0.03421, 0.03929]
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable  share    lower    upper
#>       d1     0.007371 0.2011 0.006870 0.007917
#>       d2     0.023488 0.6409 0.021960 0.025136
#>       d3     0.005791 0.1580 0.005385 0.006237
#> 
#> Diagnostics: residual 2.22e-16, condition number 1.68, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
```

## Estimates adjusted for other diseases

If a source estimate comes from a model that already included other
diseases, the closest snapshot estimand is `"snapshot_stratified"`, with
those diseases in `adjusted_for` (ids separated by `";"`), or `"all"`
for every other disease. The model value of such an estimate is computed
within strata of its adjustment set and combined across strata with
Mantel-Haenszel-type weights. With `adjusted_for = "all"`, a hazard
ratio is the disease’s own hazard multiplier and is used as it is. A Cox
coefficient adjusted for the same diseases matches this only under the
snapshot model’s assumptions:

``` r

cull_adj <- cm_model(pop, cm_impacts(
  c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), measure = "HR",
  estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_stratified"),
  adjusted_for = c("d2", NA, "all")
))
deconflate(cull_adj, event_model = TRUE, overall_risk = 0.25, n_draws = 0)$adjusted
#>   disease measure raw adjusted       change            estimand adjusted_for
#> 1      d1      HR 1.5 1.509402  0.006268033 snapshot_stratified           d2
#> 2      d2      HR 2.0 1.821704 -0.089148021      snapshot_crude         <NA>
#> 3      d3      HR 1.3 1.300000  0.000000000 snapshot_stratified          all
```

The estimate of d1 is a ratio within strata of d2 that still averages
over d3, which it is not adjusted for. Hazard ratios are not
collapsible, so it differs slightly from d1’s own hazard multiplier, the
adjusted HR (1.509 against the estimate of 1.5).

## Risk attributable to disease

The risk attributable to disease is the overall risk minus the risk of a
disease-free animal, `1 - exp(-h0)`. An animal’s risk cannot exceed 1,
so an animal with several diseases is counted once, and the attributable
risk is smaller than the sum of per-disease excess risks when diseases
co-occur. It is allocated to diseases by Shapley values over all disease
combinations
([`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md));
the allocation adds up to the attributable risk (`unallocated` is zero
up to rounding):

``` r

res$attributable$summary
#>   overall_risk disease_free_risk attributable attributable_fraction
#> 1         0.25          0.213351   0.03664899             0.1465959
#>     unallocated
#> 1 -3.469447e-17
res$attributable$by_disease
#>   disease attributable     share
#> 1      d1  0.007370744 0.2011173
#> 2      d2  0.023487728 0.6408834
#> 3      d3  0.005790514 0.1579993
```

A disease whose adjusted hazard ratio is below 1 gets a negative share.
The risks are proportions of animals with the event in the period.
Converting them into numbers of animals or a value (e.g. times the
replacement price less salvage value) is left to the user.

## Comparing methods

`compare_methods(event_model = TRUE)` sets the snapshot model beside two
approximations, each with the attributable risk at `overall_risk`:

- `"snapshot"`: the model of
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md);
- `"first_order"`: the log-linear approximation `log(HR_raw) = A beta`,
  with the conflation matrix `A` of the additive model. It needs only
  the pairwise tables, so every pair needs an association;
- `"published"`: the approach of Rasmussen et al. (2024), HR - 1
  adjusted with the proportional approximation (eq. 16 of Rasmussen et
  al. 2022) and 1 added back. It is not a method of
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md);
  it is kept for comparison and reproduction only.

The first-order and published methods need hazard (or rate) ratios, and
the published approach crude estimates. Methods that cannot be run are
listed with the reason:

``` r

compare_methods(cull, event_model = TRUE, overall_risk = 0.25)
#> <cm_comparison> methods: published, first_order, snapshot
#> Impacts: culling
#> 
#> Adjusted hazard ratios:
#>  disease measure raw published first_order snapshot
#>       d1      HR 1.5      1.41        1.40     1.38
#>       d2      HR 2.0      1.91        1.89     1.89
#>       d3      HR 1.3      1.19        1.17     1.14
#> 
#> Totals:
#>       method overall_risk disease_free_risk attributable
#>    published         0.25            0.2107      0.03934
#>  first_order         0.25            0.2122      0.03785
#>     snapshot         0.25            0.2134      0.03665
compare_methods(cull_adj, event_model = TRUE, overall_risk = 0.25)$failed
#>                                                                published 
#> "The published (2024) approach is defined for crude hazard ratios only."
```

## Sensitivity tools

[`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
[`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md),
[`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md)
and
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
take `event_model = TRUE` and `overall_risk` as well. Their total is
then the attributable risk, and disease rankings use its Shapley
allocation. They report point estimates (a distribution of the overall
risk is used at its mean).
[`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md)
also varies the overall risk:

``` r

screen_associations(cull, event_model = TRUE, overall_risk = 0.25)
#> <cm_screen> 6 scenarios; baseline total 0.03665
#>   pair    status scenario  total   change rel_change max_rank_shift rank_corr
#>  d2:d3 specified OR x 0.5 0.0404  0.00375     0.1022              1       0.5
#>  d2:d3 specified   OR x 2 0.0334 -0.00327    -0.0893              0       1.0
#>  d1:d2 specified OR x 0.5 0.0392  0.00252     0.0687              0       1.0
#>  d1:d2 specified   OR x 2 0.0343 -0.00236    -0.0644              1       0.5
#>  d1:d3 specified   OR x 2 0.0354 -0.00122    -0.0332              0       1.0
#>  d1:d3 specified OR x 0.5 0.0378  0.00113     0.0309              0       1.0
#>  failed
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
sensitivity_oat(cull, event_model = TRUE, overall_risk = 0.25)
#>          input value  total_low total_high        swing  rel_swing
#> 1    impact:d3  1.30 0.02889187 0.04351877 0.0146268990 0.39910789
#> 2    impact:d2  2.00 0.03045764 0.04211735 0.0116597041 0.31814535
#> 3         risk  0.25 0.03068666 0.04186459 0.0111779319 0.30499977
#> 4    impact:d1  1.50 0.03224435 0.04056558 0.0083212307 0.22705215
#> 5      prob:d2  0.15 0.03313608 0.04027235 0.0071362672 0.19471937
#> 6      prob:d1  0.10 0.03539108 0.03791441 0.0025233336 0.06885139
#> 7      prob:d3  0.20 0.03558047 0.03778608 0.0022056091 0.06018200
#> 8  assoc:d2:d3  3.00 0.03781031 0.03573794 0.0020723758 0.05654661
#> 9  assoc:d1:d2  2.00 0.03746031 0.03599871 0.0014616022 0.03988111
#> 10 assoc:d1:d3  1.00 0.03703142 0.03632958 0.0007018458 0.01915048
```

In
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md),
the overall risk is the input `"risk"`. At what overall risk does the
attributable risk reach 0.04, and what raw hazard ratio of d1 would
change the ranking of diseases?

``` r

cm_threshold(cull, "risk", c(0.05, 0.6), conclusion = "total", target = 0.04,
             event_model = TRUE, overall_risk = 0.25)$thresholds[, c("item", "status", "threshold")]
#>    item    status threshold
#> 1 total threshold 0.2812647
cm_threshold(cull, "impact:d1", c(0.5, 3), conclusion = "rank",
             event_model = TRUE, overall_risk = 0.25)$thresholds[, c("item", "status", "threshold")]
#>       item    status threshold
#> 1 d1 vs d3 threshold  1.404700
#> 2 d1 vs d2 threshold  2.428348
```

See
[`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md)
for how threshold searches work.

## Five diseases, mixed measures

The folder `five_diseases` contains an illustrative population of five
diseases (made-up values) and a culling table, `culling.csv`, that mixes
measures and estimands:

``` r

dir <- system.file("extdata", "five_diseases", package = "deconflate")
read.csv(file.path(dir, "culling.csv"))[, c("disease", "value", "measure", "estimand", "adjusted_for")]
#>   disease value measure            estimand adjusted_for
#> 1     LAM  1.74      HR      snapshot_crude             
#> 2     MAS  1.60      HR snapshot_stratified      LAM;SCK
#> 3     MET  1.23      RR      snapshot_crude             
#> 4     SCK  1.45      OR      snapshot_crude             
#> 5      RP  0.10      RD snapshot_stratified          all
```

LAM has a crude hazard ratio, MAS a hazard ratio stratified by LAM and
SCK, MET a risk ratio and SCK an odds ratio of culling within the
lactation, and RP a risk difference stratified by all other diseases.
The risk-based estimates need the culling risk over the lactation, here
taken to be 0.25. The population includes a three-way term, which the
snapshot model uses:

``` r

five <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                       associations = file.path(dir, "associations.csv"),
                       three_way = file.path(dir, "three_way.csv"),
                       impacts = file.path(dir, "culling.csv"))
five
#> <cm_model>
#>   Diseases: 5 (LAM, MAS, MET, SCK, RP)
#>   Disease pairs: 10 (10 with an association, 0 unknown)
#>   Three-way terms: 1
#>   Impacts: culling (event impacts: use event_model = TRUE)
#>   Estimands: snapshot_crude: 3; snapshot_stratified: 2
#>   Measures: HR: 2; OR: 1; RD: 1; RR: 1
#>   Uncertain inputs (with a distribution): 15
deconflate(five, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure  raw adjusted_hr            estimand
#>      LAM      HR 1.74       1.583      snapshot_crude
#>      MAS      HR 1.60       1.579 snapshot_stratified
#>      MET      RR 1.23       1.081      snapshot_crude
#>      SCK      OR 1.45       1.191      snapshot_crude
#>       RP      RD 0.10       1.502 snapshot_stratified
#> 
#> Overall risk 0.25; disease-free risk 0.1802; attributable to disease 0.06975 (27.9% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable   share
#>      LAM     0.024872 0.35658
#>      MAS     0.025619 0.36730
#>      MET     0.001570 0.02251
#>      SCK     0.012331 0.17679
#>       RP     0.005358 0.07682
#> 
#> Diagnostics: residual 5.90e-15, condition number 5.07, sign changes 0
```

The tables have distributions (the `dist` columns), so the default call,
without `n_draws = 0`, would also report 95% intervals from 1000 draws.
Only the snapshot model applies to these measures:

``` r

compare_methods(five, event_model = TRUE, overall_risk = 0.25)$failed
#>                                                                                                                    published 
#>   "The published method needs hazard or rate ratios; risk ratios, odds ratios and risk differences need the snapshot model." 
#>                                                                                                                  first_order 
#> "The first-order method needs hazard or rate ratios; risk ratios, odds ratios and risk differences need the snapshot model."
```

## The global dairy example

`example_global_dairy("culling")` contains the culling HRs of Rasmussen
et al. (2024) on the 12-disease population of that paper. They are
published Cox HRs, entered as `"snapshot_crude"`; that is an assumption,
and the results below hold under it. The global average replacement rate
is 23.66% (Rasmussen et al. 2024, Table 1). Pairs not in the paper’s
Table 3 are unknown, and the joint distribution fills them in. It is
fitted once and passed on with `joint`, so that it can be reused below:

``` r

gd_cull <- example_global_dairy("culling")
j_gd <- fit_joint(gd_cull)
gd <- deconflate(gd_cull, event_model = TRUE, overall_risk = 0.2366, n_draws = 0, joint = j_gd)
#> Warning: Adjusted hazard ratios are on the other side of 1 from the raw
#> estimates for CK, DYS, MET. The raw estimates are weaker than the associated
#> diseases alone would produce; check the estimands and source populations.
gd
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure   raw adjusted_hr
#>       CK      HR 1.500      0.9771
#>       CM      HR 2.300      1.8182
#>       DA      HR 2.851      1.7870
#>      DYS      HR 1.258      0.9540
#>      LAM      HR 1.745      1.2401
#>      MET      HR 1.116      0.6928
#>       MF      HR 3.000      2.4247
#>       OC      HR 1.620      1.4022
#>      PTB      HR 2.311      1.9208
#>       RP      HR 1.600      1.1832
#>      SCK      HR 1.920      1.6452
#>      SCM      HR 1.450      1.1547
#> 
#> Overall risk 0.2366; disease-free risk 0.128; attributable to disease 0.1086 (45.9% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable     share
#>       CK   -0.0001265 -0.001164
#>       CM    0.0304920  0.280670
#>       DA    0.0027374  0.025197
#>      DYS   -0.0004590 -0.004225
#>      LAM    0.0087303  0.080360
#>      MET   -0.0053919 -0.049630
#>       MF    0.0047215  0.043460
#>       OC    0.0065598  0.060381
#>      PTB    0.0138246  0.127252
#>       RP    0.0035096  0.032305
#>      SCK    0.0341130  0.313999
#>      SCM    0.0099293  0.091396
#> 
#> Diagnostics: residual 7.33e-15, condition number 3.82, sign changes 3 (CK, DYS, MET)
#> Unknown pairs: 28 without an association; the global fit gave them odds ratios from 1.08 to 2.37 (see $unknown_pairs).
```

For clinical ketosis, dystocia and metritis, the adjusted HRs are below
1: under the multiplicative model, their raw HRs are smaller than the
associated diseases alone would produce, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
warns about this. Their shares of the attributable risk are negative.

The 2024 paper instead converted each adjusted HR to an excess risk
relative to the overall rate, `HR * r / (HR * r + 1 - r) - r`, and
summed these excess risks weighted by prevalence. Applied to the
adjusted HRs above, that sum is larger than the attributable risk,
because it counts animals with several diseases more than once:

``` r

p <- gd_cull$diseases$prob
h <- gd$adjusted$adjusted
c(excess_risk_sum = sum(p * (h * 0.2366 / (h * 0.2366 + 1 - 0.2366) - 0.2366)),
  attributable = gd$attributable$summary$attributable)
#> excess_risk_sum    attributable 
#>       0.1173110       0.1086402
```

The paper treated unlisted pairs as independent and used the published
approximation;
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
recomputes its tables
([`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).

## Event impacts alongside additive impacts

Event impacts and additive impacts are adjusted in separate calls: the
additive analyses return impacts in their own units, and the event model
returns hazard ratios and a proportion of animals. To set them side by
side, convert both to a common unit, here a value per cow and year with
the global average yield of 5013 kg, a milk price of 0.5981 per kg, and
a replacement price less culled-cow price of 1299.33 - 785.86 (Rasmussen
et al. 2024, Table 1). The yield analysis has the same population, so
the joint distribution is reused. The yield impacts are percent losses,
so the gap is computed as in
[`vignette("deconflate")`](https://rasmussenphilip.github.io/deconflate/articles/deconflate.md):

``` r

yres <- deconflate(example_global_dairy("yield"), n_draws = 0, joint = j_gd)
#> Warning: Adjusted impacts change sign for CK, CM, DA, MF. The raw impacts are
#> smaller than the associated diseases alone would produce under the additive
#> model; check the estimands, and whether the estimates come from populations
#> with different comorbidity patterns.
yield_gap <- 5013 / (1 - yres$totals$adjusted_total / 100) - 5013
c(yield = yield_gap * 0.5981,
  culling = gd$attributable$summary$attributable * (1299.33 - 785.86))
#>     yield   culling 
#> 222.73500  55.78346
```

The warning says that, with the unknown pairs filled in by the global
fit, the raw yield impacts of four diseases are smaller than their
associated diseases alone would produce under the additive model.

## Methods that are not available for new analyses

The published analyses converted hazard ratios in ways that the package
does not offer for new analyses:

- Rasmussen et al. (2022) treated each HR as an odds ratio of a 2x2
  table of disease by culling, to obtain an excess annual culling risk,
  and recovered adjusted HRs by rescaling (eq. 23);
- Rasmussen et al. (2024) adjusted HR - 1 as if it were an additive
  impact, and valued culling with the excess risk relative to the
  overall rate given above.

These conversions are used only inside
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
and
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md),
so that the published tables can be recomputed (see
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).

## Caveats

- The model assumes a constant hazard within the period and that the
  estimates apply for the whole period. Most source estimates come from
  models in which disease status changes during the period.
- Hazard ratios and odds ratios are not collapsible: even without
  confounding, an average ratio differs from the ratio within subgroups.
  The snapshot model accounts for this given the distribution of disease
  combinations; the first-order method and the published approach do
  not.
- The snapshot model and the attributable risk use the whole
  distribution of disease combinations, so they depend on the
  maximum-entropy assumption, on how unknown pairs are filled in, and on
  any three-way terms.
- Culling and death compete. Make sure the source estimates are of the
  same kind (cause-specific or subdistribution), and that `overall_risk`
  is the risk of the same event.
