# Getting started with deconflate

``` r

library(deconflate)
```

## The problem

Single-disease impact estimates compare animals with and without a
disease. Animals with one disease are more (or less) likely to have
associated diseases, so each estimate partly reflects the impacts of
those other diseases. Adding such estimates double counts. `deconflate`
adjusts (“de-conflates”) the estimates so that they can be aggregated.

The package assumes no species, outcome or unit. The examples use dairy
cattle and milk yield, but any population of animals (or people) and any
additive outcome works the same way.

## A population and an impact table

The population holds what every outcome shares: the disease
probabilities and the pairwise associations. This is the three-disease
example from the Supplementary File of Rasmussen et al. (2022):

``` r

pop <- cm_population(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20), type = "prevalence"),
  cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3), measure = "OR")
)
pop
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (3 with an association, 0 unknown)
```

Associations are odds ratios by default (risk ratios, risk differences,
conditional probabilities and phi coefficients are also accepted; see
[`?cm_associations`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)).
An odds ratio of 1, as for d1 and d3, states that two diseases are
unrelated. A pair with no association at all is *unknown*:
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
then fills it in from the other associations (below).

A model adds one impact table: one raw estimate per disease, all in the
same units. The units can be anything additive (percent of yield, kg of
milk, days, euros, a welfare score); the engine does not convert them,
and the results come back in the units supplied.

``` r

yield <- cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5),
                    label = "milk yield loss", units = "% of yield")
m <- cm_model(pop, yield)
m
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (3 with an association, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Estimands: crude: 3
```

The same impacts in kg give the same adjustment, in kg:

``` r

kg <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(250, 500, 750),
                               units = "kg/animal/year"))
deconflate(kg)$adjusted
#>   disease raw adjusted      change estimand adjusted_for
#> 1      d1 250 214.3250 -0.14269981    crude         <NA>
#> 2      d2 500 338.7080 -0.32258408    crude         <NA>
#> 3      d3 750 693.4209 -0.07543874    crude         <NA>
```

Each
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
call adjusts one impact table. For another outcome, make another model
on the same population and call
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
again:

``` r

ci <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(4, 10, 2),
                               label = "calving interval increase", units = "days"))
deconflate(ci)$adjusted
#>   disease raw  adjusted      change estimand adjusted_for
#> 1      d1   4 2.9784447 -0.25538884    crude         <NA>
#> 2      d2  10 9.6989318 -0.03010682    crude         <NA>
#> 3      d3   2 0.3798539 -0.81007307    crude         <NA>
```

The tables can also be read from CSV files (with any names) or data
frames with
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md),
which checks them and reports every problem at once (see
[`vignette("own-data")`](https://rasmussenphilip.github.io/deconflate/articles/own-data.md)):

``` r

ex <- system.file("extdata", "five_diseases", package = "deconflate")
five <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
                       associations = file.path(ex, "associations.csv"),
                       impacts = file.path(ex, "impacts_yield.csv"))
five
#> <cm_model>
#>   Diseases: 5 (LAM, MAS, MET, SCK, RP)
#>   Disease pairs: 10 (10 with an association, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Estimands: crude: 5
#>   Uncertain inputs (with a distribution): 15
```

## Adjusting impacts

``` r

res <- deconflate(m)
res
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    2.143 -0.14270
#>       d2 5.0    3.387 -0.32258
#>       d3 7.5    6.934 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> 
#> Notes:
#> * No input has a distribution, so no draws were run: the results are point
#>   estimates.
```

[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
has two exact methods, chosen automatically (`method = "auto"`, the
default):

- the **simultaneous** method solves the system of impact equations from
  the pairwise 2x2 tables. It works for any number of diseases;
- the **global** method fits the maximum-entropy distribution of disease
  combinations and solves the same equations from it. It also handles
  unknown pairs, impact interactions and three-way terms
  ([`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)).
  Up to 20 diseases it enumerates all combinations; beyond that it uses
  the sampled backend of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  ([`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md)).

When every pair has an association and there are no interactions or
three-way terms, the two give the same answer and the simultaneous
method is used. Otherwise the global method is used, and a note says
why. Leaving out the d1:d3 odds ratio makes that pair unknown:

``` r

pop_unknown <- cm_population(pop$diseases,
                             cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)))
res_unknown <- deconflate(cm_model(pop_unknown, yield))
res_unknown$notes
#> [1] "The global method was used because of 1 pair without an association (unknown)."     
#> [2] "No input has a distribution, so no draws were run: the results are point estimates."
res_unknown$unknown_pairs
#>   disease1 disease2 fitted_or
#> 1       d1       d3  1.143088
rbind(odds_ratio_1 = setNames(res$adjusted$adjusted, res$adjusted$disease),
      unknown = res_unknown$adjusted$adjusted)
#>                    d1       d2       d3
#> odds_ratio_1 2.143250 3.387080 6.934209
#> unknown      1.988942 3.404334 6.906626
```

The global fit gives d1 and d3 a positive association (an odds ratio of
about 1.14), implied by their associations with d2, so the adjusted
impacts differ from those with an odds ratio of 1. Unknown is not the
same as independent; see
[`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md).

[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
sets the methods side by side, without switching:

``` r

compare_methods(m)
#> <cm_comparison> methods: published, simultaneous, global
#> Impacts: milk yield loss
#> Units: % of yield
#> 
#> Adjusted values:
#>  disease raw published simultaneous global
#>       d1 2.5      2.07         2.14   2.14
#>       d2 5.0      3.70         3.39   3.39
#>       d3 7.5      6.75         6.93   6.93
#> 
#> Totals:
#>        method raw_sum adjusted_total
#>     published     2.5          2.111
#>  simultaneous     2.5          2.109
#>        global     2.5          2.109
```

The `"published"` column is the proportional approximation of Rasmussen
et al. (2022, 2024), eq. 16. It is not a method of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md):
it is kept for comparison
([`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md))
and for reproducing the published analyses
([`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md))
only. It is defined for crude estimates only.

The diagnostics report how well the adjusted impacts reconstruct the raw
ones, the condition number of the conflation matrix, and adjusted
impacts whose sign differs from the raw impact. A sign change means that
the raw impact is smaller than the associated diseases alone would
produce, so the inputs are inconsistent under the additive model (or the
estimands are not what they were assumed to be):

``` r

low_d2 <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(2.5, 0.2, 7.5)))
tryCatch(deconflate(low_d2), deconflate_sign_change = function(w) conditionMessage(w))
#> [1] "Adjusted impacts change sign for d2. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns."
```

## Contributions

The expected aggregate impact per animal is `sum_i p_i b_i` (plus any
interaction terms). Each disease’s contribution is its own term plus
half of each interaction term it is involved in, which is its Shapley
value in closed form. Contributions add up to the aggregate:

``` r

res$contributions
#>   disease      main interaction     total     share
#> 1      d1 0.2143250           0 0.2143250 0.1016130
#> 2      d2 0.5080619           0 0.5080619 0.2408757
#> 3      d3 1.3868419           0 1.3868419 0.6575113
res$totals
#>   raw_sum adjusted_total interaction_total
#> 1     2.5       2.109229                 0
```

`raw_sum` is the naive aggregate of the unadjusted impacts. The shares
are `NA` when the aggregate is zero, since a share of zero is not
defined.
[`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) combine the
impacts and contributions in one table.

How robust is a ranking or a sign to one uncertain input?
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
finds the input values at which the conclusion changes, for example the
odds ratio of d2 and d3 at which d1 would contribute more than d2:

``` r

cm_threshold(m, "assoc:d2:d3", c(0.1, 100), conclusion = "rank")$thresholds[
  , c("item", "status", "threshold", "description")]
#>       item    status threshold                         description
#> 1 d1 vs d2 threshold  10.20365 d1 moves above d2 (by contribution)
```

See
[`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md).

## Estimands

Each raw impact must be one of two estimands, set explicitly in
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md):

- `"crude"` (default): the difference in the outcome between animals
  with and without the disease, not adjusted for other diseases.
- `"adjusted_linear"`: the coefficient of the disease in an additive
  (linear) regression of the outcome on the disease and the diseases
  listed in `adjusted_for` (ids separated by `";"`, or `"all"` for every
  other disease).

`adjusted_for` is never used to infer the estimand. Giving it for a
crude estimate is an error:

``` r

tryCatch(cm_impacts("d1", 2, adjusted_for = "d2"),
         deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "`adjusted_for` is given for crude estimates (d1). Set estimand = 'adjusted_linear' if these estimates are adjusted for those diseases; other adjusted estimands are not supported."
```

Other adjusted estimands (matched or propensity-score estimates, for
example) are not supported. For an adjusted estimate, the package uses
the population linear projection of the omitted diseases, which depends
only on the pairwise tables. This requires that the probabilities and
associations describe the source population of the estimates.

A worked example: in the population above, the true additive impacts are
0.02, 0.04 and 0.06. A study of d1 that adjusted for d2, and crude
studies of d2 and d3, would report:

``` r

truth <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)
raw <- simulate_raw_impacts(pop, truth,
                            estimand = c("adjusted_linear", "crude", "crude"),
                            adjusted_for = c("d2", NA, NA))
round(raw$value, 4)
#> [1] 0.0187 0.0541 0.0667
```

With the estimands recorded, the adjustment recovers the truth:

``` r

imp_adj <- cm_impacts(c("d1", "d2", "d3"), raw$value,
                      estimand = c("adjusted_linear", "crude", "crude"),
                      adjusted_for = c("d2", NA, NA))
m_adj <- cm_model(pop, imp_adj)
deconflate(m_adj)$adjusted
#>   disease        raw adjusted      change        estimand adjusted_for
#> 1      d1 0.01866481     0.02  0.07153525 adjusted_linear           d2
#> 2      d2 0.05406438     0.04 -0.26014129           crude         <NA>
#> 3      d3 0.06668175     0.06 -0.10020359           crude         <NA>
```

Treating the d1 estimate as crude adjusts it a second time for d2 and
gives the wrong answer:

``` r

m_wrong <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), imp_adj$value))
deconflate(m_wrong)$adjusted$adjusted
#> [1] 0.01440635 0.04043097 0.05992801
```

The published approximation is not defined for adjusted estimands;
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
records the reason:

``` r

compare_methods(m_adj)$failed
#>                                                                                                                          published 
#> "The published approximation is defined for crude estimates only; use method = 'simultaneous' or 'global' for adjusted estimands."
```

## Uncertainty

When inputs have distributions (the `dist` columns of the input tables,
or `distributions` in
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)),
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
also draws them (`n_draws = 1000` by default) and adds 95% intervals.
The central estimate always uses the point values. Without distributions
no draws are run, and a note says so (as in the results above). The
five-disease example has distributions for most inputs; 200 draws keep
this example fast:

``` r

res_five <- deconflate(five, n_draws = 200, seed = 1)
res_five
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease raw adjusted  lower upper  change
#>      LAM 4.8    4.139 2.1878 5.863 -0.1377
#>      MAS 3.3    2.356 0.8542 4.515 -0.2861
#>      MET 5.6    4.514 2.6746 6.125 -0.1939
#>      SCK 2.5    1.366 0.6589 1.930 -0.4535
#>       RP 4.2    2.665 1.6610 3.299 -0.3654
#> 
#> Raw sum: 3.742; adjusted total: 2.735 (95% interval 2.208 to 3.592)
#> Diagnostics: residual 0.00e+00, condition number 1.78, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
```

Pass `n_draws = 0` for point estimates only. See
[`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md)
for the draws, their stability checks and the sensitivity tools.

## Event impacts

Impacts on the risk of an event, such as culling or death, are hazard
ratios, rate ratios, risk ratios, odds ratios or risk differences rather
than additive amounts. They are given with a `measure` (one per disease;
measures can be mixed) and a snapshot estimand, and adjusted with
`event_model = TRUE` and the overall risk of the event over the period:

``` r

cull <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.2),
                                 measure = c("HR", "HR", "RR"),
                                 estimand = "snapshot_crude", label = "culling"))
res_cull <- deconflate(cull, event_model = TRUE, overall_risk = 0.25)
res_cull$adjusted
#>   disease measure raw adjusted      change       estimand adjusted_for
#> 1      d1      HR 1.5 1.381707 -0.07886199 snapshot_crude         <NA>
#> 2      d2      HR 2.0 1.908294 -0.04585287 snapshot_crude         <NA>
#> 3      d3      RR 1.2 1.097741          NA snapshot_crude         <NA>
res_cull$attributable$summary
#>   overall_risk disease_free_risk attributable attributable_fraction unallocated
#> 1         0.25         0.2147648   0.03523521             0.1409408           0
```

The adjusted values are each disease’s own hazard ratio, and
`$attributable` gives the risk attributable to disease and its
allocation to diseases. See
[`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md).

## From adjusted impacts to gaps and values

[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
returns impacts and contributions in the units of the impacts; it does
not convert them into other quantities. Turning them into a gap between
the observed and the disease-free level of an outcome, or into a value,
needs to know what the impacts mean, and takes a few lines of base R.

When the impacts are percentage losses relative to the disease-free
level (as in `m`, with an observed mean of 10,000 units), the
disease-free level is `observed / (1 - aggregate / 100)`, the gap is the
disease-free level minus the observed level, and each disease’s part of
the gap is its share of the aggregate:

``` r

observed <- 10000
aggregate <- res$totals$adjusted_total            # in percent
disease_free <- observed / (1 - aggregate / 100)
gap <- disease_free - observed
ct <- res$contributions
by_disease <- data.frame(disease = ct$disease, gap = gap * ct$share)
by_disease
#>   disease       gap
#> 1      d1  21.89431
#> 2      d2  51.90090
#> 3      d3 141.67238
c(disease_free = disease_free, gap = gap)
#> disease_free          gap 
#>   10215.4676     215.4676
```

A value per unit of outcome turns the gaps into values; any costs that
are not part of the adjustment are added separately:

``` r

unit_value <- 0.30
by_disease$value <- by_disease$gap * unit_value
sum(by_disease$value)
#> [1] 64.64028
```

`ct$share` is `NA` when the aggregate is zero;
`disease_free * ct$total / 100` gives the same gaps and stays defined.
For an outcome that disease increases, the disease-free level is
`observed / (1 + aggregate / 100)` and the gap is
`observed - disease_free`. When the impacts are absolute effects (in the
outcome’s own units, such as the days of the calving-interval model `ci`
above), no conversion is needed: each disease’s contribution
(`ct$total`) is its part of the gap, and the aggregate is the gap. With
draws, the gap increases with the aggregate, so the same formula applied
to the interval bounds of the aggregate (`adjusted_total_lower` and
`adjusted_total_upper` in `$totals`) gives an interval for the total
gap.

## Feasibility

Each pairwise table can be valid on its own while no population has all
of them. By default the simultaneous method screens every triple of
diseases (a necessary condition), and `feasibility = "lp"` runs an exact
check (it needs the `lpSolve` package):

``` r

res$diagnostics$feasibility
#> [1] "triple screen passed (necessary condition only)"
```

Three diseases with probability 0.5, two strong positive associations
and one strong negative one cannot coexist:

``` r

bad <- cm_model(
  cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
  cm_impacts(c("a", "b", "c"), c(1, 2, 3)),
  associations = cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05))
)
check_feasibility(bad, method = "triples")
#> <cm_feasibility> NOT jointly feasible [triples]
#> 
#> Conflicting triples:
#>  disease1 disease2 disease3   gap
#>         a        b        c 0.226
```

## Conditions

Errors and warnings have classes, so that callers (and the draws and
sensitivity tools of the package) can handle them separately:

| Class | Meaning |
|----|----|
| `deconflate_infeasible` | inputs that admit no valid probabilities (e.g. the example above) |
| `deconflate_singular` | a singular or non-identifiable system |
| `deconflate_nonconvergence` | a numerical procedure did not converge |
| `deconflate_unsupported` | an estimand and method combination that is not supported |
| `deconflate_nonfinite` | non-finite results |
| `deconflate_input_problems` | [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md) found errors in the input tables |

All errors also have class `deconflate_error`, and warnings
`deconflate_warning` (e.g. `deconflate_sign_change`).

``` r

e <- tryCatch(deconflate(bad), deconflate_error = function(e) e)
class(e)
#> [1] "deconflate_infeasible" "deconflate_error"      "error"                
#> [4] "condition"
conditionMessage(e)
#> [1] "The pairwise associations are jointly infeasible: no population has all of them (conflicting triples: a-b-c). See check_feasibility()."
```

[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
also needs at least one association; without any, there is nothing to
de-conflate (the diseases would be treated as independent), and it
points to the tools that show how much associations could matter:

``` r

no_assoc <- cm_model(pop$diseases, yield)
tryCatch(deconflate(no_assoc), deconflate_error = function(e) conditionMessage(e))
#> [1] "No association estimates were given, so there is nothing to de-conflate (the diseases would be treated as independent). To see how much associations could change the results, use screen_associations() or cm_threshold()."
```

## Checking the method with a simulated herd

[`simulate_raw_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/simulate_raw_impacts.md)
computes the raw impacts that single-disease studies would report in a
herd with known true impacts. Adjusting them should recover the truth;
the published approximation (shown for comparison) does so only
approximately:

``` r

m_sim <- cm_model(pop, simulate_raw_impacts(pop, c(d1 = 2, d2 = 4, d3 = 6)))
compare_methods(m_sim, methods = c("published", "simultaneous"))
#> <cm_comparison> methods: published, simultaneous
#> 
#> Adjusted values:
#>  disease  raw published simultaneous
#>       d1 2.42      1.96            2
#>       d2 5.41      4.19            4
#>       d3 6.67      5.87            6
#> 
#> Totals:
#>        method raw_sum adjusted_total
#>     published   2.387          1.998
#>  simultaneous   2.387          2.000
```

## Where next

- [`vignette("own-data")`](https://rasmussenphilip.github.io/deconflate/articles/own-data.md):
  reading your own inputs from CSV files or data frames.
- [`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md):
  culling and other event impacts, and the risk attributable to disease.
- [`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md):
  the 2022 and 2024 analyses.
- [`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md):
  the joint distribution, unknown associations, impact interactions,
  three-way scenarios and attribution.
- [`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md):
  draws and intervals, stability checks and sensitivity screening.
- [`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md):
  threshold searches with
  [`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
  (at what input value does a ranking, a sign or the aggregate change?),
  and the sampled backend of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  for the global method with more than about 20 diseases.
