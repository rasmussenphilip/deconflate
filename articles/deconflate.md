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

## A population and an impact vector

The population holds what all analyses share: the disease probabilities
and the pairwise associations. This is the three-disease example from
the Supplementary File of Rasmussen et al. (2022):

``` r

pop <- cm_population(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20), type = "prevalence"),
  cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3), measure = "OR")
)
pop
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
```

Pairs without an association are independent by default
(`missing_associations = "independent"`). Use `"unknown"` to leave them
unconstrained instead; this is a different assumption (see
[`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)).

An analysis adds one impact vector: one value per disease, all in the
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
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
#>   Impacts: milk yield loss [% of yield]
#>   Estimands: crude: 3
```

The same impacts in kg give the same adjustment, in kg:

``` r

kg <- cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(250, 500, 750),
                               units = "kg/cow/year"))
deconflate(kg)$adjusted
#>   disease raw adjusted      change estimand adjusted_for
#> 1      d1 250 214.3250 -0.14269981    crude         <NA>
#> 2      d2 500 338.7080 -0.32258408    crude         <NA>
#> 3      d3 750 693.4209 -0.07543874    crude         <NA>
```

Several types of impact are separate analyses on one population:

``` r

a <- cm_analyses(pop,
  yield = yield,
  calving_interval = cm_impacts(c("d1", "d2", "d3"), c(4, 10, 2), units = "days")
)
a
#> <cm_analyses> 2 analyses on one population
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
#>   - yield [% of yield]
#>   - calving_interval [days]
```

`deconflate(a)` adjusts each analysis and returns a named list of
results.

## Adjusting impacts

``` r

compare_methods(m)
#> <cm_comparison> methods: published, simultaneous, global
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

[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
has two methods:

- `"simultaneous"` (the default) solves the system of impact equations
  exactly, using the pairwise 2x2 tables. Use it when impacts are
  additive and every pair of diseases has an association estimate or a
  defensible independence assumption.
- `"global"` fits the maximum-entropy distribution of disease
  combinations and can include impact interactions
  ([`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)).
  Use it otherwise: when impacts interact, or when some pairs are
  unknown. Without interactions and unknown pairs, it equals
  `"simultaneous"`. The simultaneous method works for any number of
  diseases; the global method enumerates combinations up to about 20
  diseases, and beyond that uses the sampled backend of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  ([`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md)).

The `"published"` row above is the proportional approximation of
Rasmussen et al. (2022, 2024), eq. 16. It is not a method of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md):
it is kept for comparison
([`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md),
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md))
and for reproducing the published analyses
([`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md))
only. It is defined for crude estimates only. (Hazard ratios have their
own adapter, with the snapshot model as the default; see
[`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md).)

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
```

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
#> [1] "`adjusted_for` is given for crude estimates (d1). Set estimand = 'adjusted_linear' if these are coefficients from an additive regression adjusted for those diseases; other adjusted estimands are not supported."
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
raw$value
#> [1] 0.01866481 0.05406438 0.06668175
```

that is, 0.018664808346, 0.054064376914 and 0.066681750858. With the
estimands recorded, the simultaneous method recovers the truth:

``` r

imp_adj <- cm_impacts(c("d1", "d2", "d3"),
                      c(0.018664808346, 0.054064376914, 0.066681750858),
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
outcome’s own units, such as the days of the calving-interval analysis
above), no conversion is needed: each disease’s contribution
(`ct$total`) is its part of the gap, and the aggregate is the gap.

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

Errors and warnings have classes, so that callers (and the Monte Carlo
and screening functions) can handle them separately:

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
  reading your own inputs from CSV files.
- [`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md):
  culling hazard ratios and the culling attributable to disease.
- [`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md):
  the 2022 and 2024 analyses.
- [`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md):
  the joint distribution, unknown associations, impact interactions,
  three-way scenarios and attribution.
- [`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md):
  Monte Carlo analysis, stability checks, scenarios and sensitivity
  screening.
- [`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md):
  threshold searches with
  [`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
  (at what input value does a ranking, a sign or the aggregate change?),
  and the sampled backend of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  for the global method with more than about 20 diseases.
