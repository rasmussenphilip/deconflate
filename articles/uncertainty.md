# Uncertainty and sensitivity

``` r

library(deconflate)
```

Any input can be uncertain: a disease probability, an association, a
three-way term, an impact, an interaction, and, for event impacts, the
overall risk.
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
draws the uncertain inputs from their distributions, adjusts each draw
in the same way as the central estimate, and reports intervals. Monte
Carlo runs are part of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
(`n_draws`).

## Input distributions

`dist_*()` functions describe uncertain inputs. Draws are on the scale
the inputs were entered on (e.g. percent impacts, odds ratios, incidence
rates).

``` r

dist_pert(1.19, 3.30, 10.71)        # min, mode, max
#> <cm_dist> pert(min = 1.19, mode = 3.3, max = 10.71, lambda = 4), mean 4.183
dist_lognormal_ci(2.7, 1.5, 4.9)    # odds ratio and 95% CI
#> <cm_dist> lognormal(meanlog = 0.9933, sdlog = 0.302), mean 2.826
dist_normal(2.63, 1.44, lower = 0)  # truncated at zero
#> <cm_dist> normal(mean = 2.63, sd = 1.44, lower = 0, upper = Inf), mean 2.742
```

The others are
[`dist_fixed()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
[`dist_lognormal()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
[`dist_beta()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
[`dist_pert_mean()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md),
[`dist_uniform()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
and
[`dist_mixture()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
(e.g. two pooled sources of evidence); see
[`?distributions`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md).

## Where distributions are given

In input tables, an uncertain value gets its distribution in its own
row, in the columns `dist` and `p1`-`p4`, and
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
keeps it in the model (see
[`vignette("own-data")`](https://rasmussenphilip.github.io/deconflate/articles/own-data.md)).
In R, give them to
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
as `distributions`, a named list whose names are keys: `prob:<disease>`,
`assoc:<d1>:<d2>`, `three:<d1>:<d2>:<d3>`, `impact:<disease>` and
`inter:<d1>:<d2>`. Inputs without a distribution keep their point value
in every draw.

``` r

m <- example_supplement()
mu <- cm_model(m, m$impacts, distributions = list(
  "prob:d1" = dist_beta(20, 180),
  "assoc:d1:d2" = dist_lognormal_ci(2, 1.4, 2.9),
  "impact:d1" = dist_normal(2.5, 0.5),
  "impact:d2" = dist_normal(5, 1),
  "impact:d3" = dist_pert(5, 7.5, 9)
))
mu
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (3 with an association, 0 unknown)
#>   Impacts: yield [%] (additive)
#>   Estimands: crude: 3
#>   Uncertain inputs (with a distribution): 5
```

## Draws and intervals

When any input has a distribution,
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
runs `n_draws` draws (1000 by default; here 200 to keep this vignette
fast). Each draw is a complete set of inputs, so the uncertainty of
probabilities, associations and impacts is propagated jointly:

``` r

r <- deconflate(mu, n_draws = 200, seed = 1)
r
#> <cm_result> method: simultaneous; yield [%]
#> 
#>  disease raw adjusted  lower upper   change
#>       d1 2.5    2.143 0.9455 3.134 -0.14270
#>       d2 5.0    3.387 1.5783 5.325 -0.32258
#>       d3 7.5    6.934 5.0366 8.176 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109 (95% interval 1.668 to 2.454)
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
```

The central estimate (`adjusted`) always uses the point values; the
draws give the 95% intervals (`lower`, `upper`: the 2.5% and 97.5%
quantiles over the accepted draws). The intervals are added to the
adjusted impacts, the contributions and the adjusted total:

``` r

r$contributions
#>   disease      main interaction     total     share      lower     upper
#> 1      d1 0.2143250           0 0.2143250 0.1016130 0.08916943 0.3660317
#> 2      d2 0.5080619           0 0.5080619 0.2408757 0.23674792 0.7987588
#> 3      d3 1.3868419           0 1.3868419 0.6575113 1.00731482 1.6352118
r$totals
#>   raw_sum adjusted_total interaction_total adjusted_total_lower
#> 1     2.5       2.109229                 0             1.667907
#>   adjusted_total_upper
#> 1             2.453791
```

`$draws` holds the details: `summary` (mean, standard deviation, Monte
Carlo standard error `mcse`, interval, median and a stability check for
every quantity), the accepted draws (`values`, and the drawn inputs in
`params`), the rejections and the seed:

``` r

r$draws$summary
#>          quantity      mean         sd        mcse      lower    median
#> 1     adjusted:d1 2.1675643 0.54205131 0.038328816 0.94554635 2.1903500
#> 2     adjusted:d2 3.4613776 1.01702150 0.071914280 1.57831944 3.4670986
#> 3     adjusted:d3 6.7080184 0.86100899 0.060882529 5.03657409 6.7445200
#> 4 contribution:d1 0.2195163 0.06943075 0.004909496 0.08916943 0.2154583
#> 5 contribution:d2 0.5192066 0.15255323 0.010787142 0.23674792 0.5200648
#> 6 contribution:d3 1.3416037 0.17220180 0.012176506 1.00731482 1.3489040
#> 7           total 2.0803266 0.20281949 0.014341504 1.66790746 2.0800081
#> 8         raw_sum 2.4669067 0.24419083 0.017266899 1.97123335 2.4596112
#>       upper stability
#> 1 3.1335943        ok
#> 2 5.3250588        ok
#> 3 8.1760588        ok
#> 4 0.3660317        ok
#> 5 0.7987588        ok
#> 6 1.6352118        ok
#> 7 2.4537906        ok
#> 8 2.9201849        ok
```

The stability check is `ok`, `imprecise` (the Monte Carlo standard error
is more than 5% of the mean: use more draws, or Latin hypercube
sampling), `heavy_tail` (the most extreme 1% of draws contribute more
than 60% of the variance) or `insufficient_info` (the precision cannot
be assessed, e.g. because the mean is zero; `mcse` still gives the
absolute error). A note lists the quantities that are not `ok`.

Without any distribution, no draws are run, whatever `n_draws`, and a
note says so; `n_draws = 0` skips the draws (and the note) for a model
that has distributions:

``` r

deconflate(example_supplement())$notes
#> [1] "No input has a distribution, so no draws were run: the results are point estimates."
```

### Seed

Draws are reproducible with `seed`. Without it, a seed is chosen and
stored in the result, so a run can always be repeated:

``` r

r_again <- deconflate(mu, n_draws = 200, seed = 1)
identical(r$adjusted, r_again$adjusted)
#> [1] TRUE
r_unseeded <- deconflate(mu, n_draws = 200)
r_unseeded$draws$seed
#> [1] 132362068
```

### Rejected draws

Some draws are impossible: a probability outside (0, 1), a non-positive
odds ratio, or associations that no population can have together. Such
draws are rejected and not replaced, so the intervals describe the
accepted draws only. The number rejected is printed with the result,
each rejection is listed by type with its reason in `$draws$rejections`,
and a note appears when more than 10% are rejected. Here a prevalence is
given a normal distribution that reaches below zero:

``` r

m_rej <- cm_model(m, m$impacts, distributions = list("prob:d1" = dist_normal(0.05, 0.05)))
r_rej <- deconflate(m_rej, n_draws = 200, seed = 1)
r_rej
#> <cm_result> method: simultaneous; yield [%]
#> 
#>  disease raw adjusted lower upper   change
#>       d1 2.5    2.143 2.108 2.159 -0.14270
#>       d2 5.0    3.387 3.324 3.542 -0.32258
#>       d3 7.5    6.934 6.908 6.945 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109 (95% interval 1.924 to 2.209)
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (22 rejected; random sampling; seed 1).
#> 
#> Notes:
#> * 11% of the draws were rejected (infeasible 22): the input distributions
#>   often combine values that cannot hold together. The intervals describe the
#>   accepted draws only; see $draws$rejections.
head(r_rej$draws$rejections, 3)
#>   draw       type
#> 1   10 infeasible
#> 2   24 infeasible
#> 3   27 infeasible
#>                                                                                        reason
#> 1  prob:d1 = -0.0269975 is not a valid disease value (it gives a probability outside (0, 1)).
#> 2 prob:d1 = -0.00738285 is not a valid disease value (it gives a probability outside (0, 1)).
#> 3   prob:d1 = -0.060735 is not a valid disease value (it gives a probability outside (0, 1)).
```

Draws can also be rejected because the values cannot hold together, even
when each is valid alone. Three diseases with strong odds ratios are
feasible at their point values, but not when one of them is much weaker
(`warn = FALSE` hides a warning that these made-up impacts change sign):

``` r

strong <- cm_model(
  cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
  cm_impacts(c("a", "b", "c"), c(4, 5, 6)),
  associations = cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 20)),
  distributions = list("assoc:a:b" = dist_lognormal_ci(20, 1, 400))
)
r_strong <- deconflate(strong, n_draws = 200, seed = 1, warn = FALSE)
r_strong$draws$n_rejected
#> [1] 15
table(r_strong$draws$rejections$type)
#> 
#> infeasible 
#>         15
r_strong$draws$rejections$reason[1]
#> [1] "The pairwise associations are jointly infeasible: no population has all of them (conflicting triples: a-b-c). See check_feasibility()."
```

Report the rejection share: conditioning on acceptance changes the
effective input distribution. A high share means the input distributions
often combine values that contradict each other; narrower or correlated
evidence may be needed.

### Latin hypercube sampling

`sampling = "lhs"` stratifies each input’s distribution, which usually
reduces the Monte Carlo error of means. Latin hypercube draws are not
independent, so they are split into `lhs_replicates` independent blocks
(10 by default), and the Monte Carlo standard error is estimated from
the spread of the block means. Compare `mcse` with the random draws
above:

``` r

r_lhs <- deconflate(mu, n_draws = 200, seed = 1, sampling = "lhs")
r_lhs
#> <cm_result> method: simultaneous; yield [%]
#> 
#>  disease raw adjusted lower upper   change
#>       d1 2.5    2.143 1.168 3.137 -0.14270
#>       d2 5.0    3.387 1.434 5.538 -0.32258
#>       d3 7.5    6.934 5.168 8.360 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109 (95% interval 1.658 to 2.39)
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; Latin hypercube sampling; seed 1).
data.frame(quantity = r$draws$summary$quantity,
           mcse_random = r$draws$summary$mcse,
           mcse_lhs = r_lhs$draws$summary$mcse)
#>          quantity mcse_random     mcse_lhs
#> 1     adjusted:d1 0.038328816 0.0033480607
#> 2     adjusted:d2 0.071914280 0.0059934476
#> 3     adjusted:d3 0.060882529 0.0041949602
#> 4 contribution:d1 0.004909496 0.0014896424
#> 5 contribution:d2 0.010787142 0.0008990171
#> 6 contribution:d3 0.012176506 0.0008389920
#> 7           total 0.014341504 0.0020168647
#> 8         raw_sum 0.017266899 0.0022140969
```

## Event impacts: the overall risk

For event impacts
([`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md)),
the impacts can have distributions like any other input, and
`overall_risk` can be a distribution too. The central estimate uses the
mean of that distribution:

``` r

cull <- cm_model(m, cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.2),
                               measure = c("HR", "HR", "RR"),
                               estimand = "snapshot_crude", label = "culling"))
e <- deconflate(cull, event_model = TRUE, overall_risk = dist_beta(250, 750),
                n_draws = 200, seed = 1)
e
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure raw adjusted_hr lower upper
#>       d1      HR 1.5       1.382 1.382 1.382
#>       d2      HR 2.0       1.908 1.906 1.910
#>       d3      RR 1.2       1.098 1.093 1.104
#> 
#> Overall risk 0.25; disease-free risk 0.2148; attributable to disease 0.03524 (14.1% of the overall risk); 95% interval [0.03276, 0.03796]
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable  share    lower    upper
#>       d1     0.007366 0.2091 0.006866 0.007912
#>       d2     0.023888 0.6780 0.022378 0.025504
#>       d3     0.003981 0.1130 0.003517 0.004549
#> 
#> Diagnostics: residual 5.27e-16, condition number 1.72, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
e$attributable$summary
#>   overall_risk disease_free_risk attributable attributable_fraction unallocated
#> 1         0.25         0.2147648   0.03523521             0.1409408           0
#>   attributable_lower attributable_upper
#> 1         0.03276159         0.03796451
```

Here only the overall risk is uncertain. The attributable risk and its
allocation to diseases vary with it, while the adjusted hazard ratios
hardly change: raw hazard ratios do not depend on the overall risk, and
only the risk ratio of d3 does.

## Comparing methods on the same draws

[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
runs several methods on the same inputs, including the approximation of
Rasmussen et al. (2022, eq. 16) (`"published"`), kept for comparison
with earlier analyses only. With `n_draws` (default 0), every method is
applied to the same draws; a draw is rejected if any method fails on it.
There is no automatic switch of method here: methods that cannot be run
are listed in `$failed`.

``` r

mu2 <- cm_model(m, m$impacts, distributions = list(
  "assoc:d1:d2" = dist_lognormal_ci(2, 1.4, 2.9),
  "impact:d2" = dist_normal(5, 1)
))
cmp <- compare_methods(mu2, methods = c("published", "simultaneous"), n_draws = 200, seed = 1)
cmp
#> <cm_comparison> methods: published, simultaneous
#> Impacts: yield
#> Units: %
#> 
#> Adjusted values:
#>  disease raw published simultaneous
#>       d1 2.5      2.07         2.14
#>       d2 5.0      3.70         3.39
#>       d3 7.5      6.75         6.93
#> 
#> Totals:
#>        method raw_sum adjusted_total
#>     published     2.5          2.111
#>  simultaneous     2.5          2.109
#> 
#> Uncertainty (200 draws, 0 rejected): mean and 95% interval of the total
#>        method  mean lower upper stability
#>     published 2.107 1.945 2.333        ok
#>  simultaneous 2.102 1.928 2.333        ok
s <- cmp$draws$summary
s[s$quantity == "adjusted:d2", c("method", "mean", "lower", "upper")]
#>          method     mean    lower    upper
#> 2     published 3.658327 2.303267 5.664957
#> 10 simultaneous 3.332220 1.759484 5.555188
```

## Gaps and values over the draws

The package reports impacts in their own units. A gap or a value per
draw is computed in base R from the accepted draws in `$draws$values`
(one row per draw, one column per quantity). Here the impacts are
percent losses of an observed level of 10,000 units, and the value of a
unit is drawn per draw:

``` r

total <- r$draws$values[, "total"]
gap <- 10000 / (1 - total / 100) - 10000
set.seed(5)
value <- gap * runif(length(gap), 0.25, 0.35)
sapply(list(gap = gap, value = value),
       function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975))))
#>            gap    value
#> mean  212.4962 63.99437
#> 2.5%  173.6852 50.78312
#> 97.5% 251.5701 81.64235
```

The central gap uses the central total:

``` r

10000 / (1 - r$totals$adjusted_total / 100) - 10000
#> [1] 215.4676
```

## Sensitivity

The sensitivity tools give point estimates (they do not run draws), and
choose the method as
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
does. They also accept `event_model = TRUE` with `overall_risk`.
One-at-a-time sensitivity (as in Rasmussen et al. 2024, Fig. 7):

``` r

sensitivity_oat(example_supplement(), variation = 0.2)
#>         input value total_low total_high      swing  rel_swing
#> 1   impact:d3  7.50  1.844070   2.374387 0.53031715 0.25142703
#> 2     prob:d3  0.20  1.851628   2.374637 0.52300893 0.24796215
#> 3   impact:d2  5.00  1.998423   2.220035 0.22161264 0.10506809
#> 4     prob:d2  0.15  2.026335   2.197276 0.17094138 0.08104449
#> 5 assoc:d2:d3  3.00  2.170224   2.062015 0.10820854 0.05130242
#> 6   impact:d1  2.50  2.063348   2.155110 0.09176175 0.04350488
#> 7     prob:d1  0.10  2.069218   2.149697 0.08047878 0.03815555
#> 8 assoc:d1:d3  1.00  2.138876   2.084124 0.05475254 0.02595856
#> 9 assoc:d1:d2  2.00  2.127728   2.093846 0.03388228 0.01606382
```

Screening associations, including pairs with no estimate, finds those
worth estimating. Specified associations are multiplied by each of
`multipliers`; pairs without an estimate (here d1 and d3) are set to
each of `or_values`:

``` r

m_pairs <- cm_model(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5)),
  associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
)
screen_associations(m_pairs, or_values = c(0.5, 2), multipliers = c(0.5, 2))
#> <cm_screen> 6 scenarios; baseline total 2.091
#>   pair    status scenario total  change rel_change max_rank_shift rank_corr
#>  d2:d3 specified OR x 0.5  2.30  0.2094     0.1001              0         1
#>  d2:d3 specified   OR x 2  1.91 -0.1769    -0.0846              0         1
#>  d1:d3   unknown OR = 0.5  2.20  0.1043     0.0499              0         1
#>  d1:d3   unknown   OR = 2  2.01 -0.0780    -0.0373              0         1
#>  d1:d2 specified   OR x 2  2.02 -0.0755    -0.0361              0         1
#>  d1:d2 specified OR x 0.5  2.16  0.0733     0.0351              0         1
#>  failed
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
```

The screens report scenarios that cannot be run, with the reason in
column `failed`, rather than skipping them. Here, weakening any of the
three strong associations above makes the pairs jointly infeasible:

``` r

screen_associations(strong, multipliers = c(0.05, 2))
#> <cm_screen> 6 scenarios; baseline total 3.305
#>  pair    status  scenario total  change rel_change max_rank_shift rank_corr
#>   b:c specified    OR x 2  3.17 -0.1309   -0.03960              0         1
#>   a:c specified    OR x 2  3.22 -0.0808   -0.02443              0         1
#>   a:b specified    OR x 2  3.27 -0.0306   -0.00927              0         1
#>   a:b specified OR x 0.05    NA      NA         NA             NA        NA
#>   a:c specified OR x 0.05    NA      NA         NA             NA        NA
#>   b:c specified OR x 0.05    NA      NA         NA             NA        NA
#>      failed
#>        <NA>
#>        <NA>
#>        <NA>
#>  infeasible
#>  infeasible
#>  infeasible
```

The screens also work without any association estimates: the baseline
then has every disease independent, and the screen shows which pairs
would matter if they were associated.
[`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
and
[`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md)
(see
[`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md))
work in the same way, and
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
finds the association at which a result changes (see
[`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md)).
