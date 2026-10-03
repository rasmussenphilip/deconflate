# Uncertainty, scenarios and sensitivity

``` r

library(deconflate)
```

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

## Monte Carlo analysis

[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
attaches distributions to a model. Keys name the inputs:
`prob:<disease>`, `assoc:<d1>:<d2>`, `three:<d1>:<d2>:<d3>`,
`impact:<disease>` and `inter:<d1>:<d2>`. (When inputs are read from
files with
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md),
distributions are given in the rows of the tables, in the columns `dist`
and `p1`-`p4`, and the sampler is built for you; see
[`vignette("own-data")`](https://rasmussenphilip.github.io/deconflate/articles/own-data.md).)
Each draw is a complete model, so the uncertainty of probabilities,
associations and impacts is propagated jointly. Here the association of
d2 and d3 is given a defensive mixture (a wider component with weight
0.2), which is used for a scenario below:

``` r

m <- example_supplement()
s <- cm_sampler(
  m,
  diseases = list(d1 = dist_beta(20, 180)),
  associations = list(
    "d1:d2" = dist_lognormal_ci(2, 1.4, 2.9),
    "d2:d3" = dist_mixture(dist_lognormal_ci(3, 2, 4.5), dist_lognormal(log(3), 0.8),
                           weights = c(0.8, 0.2))
  ),
  impacts = list(d1 = dist_normal(2.5, 0.5), d2 = dist_normal(5, 1),
                 d3 = dist_pert(5, 7.5, 9))
)
s
#> <cm_sampler> 6 uncertain inputs
#>   prob:d1: beta
#>   assoc:d1:d2: lognormal
#>   assoc:d2:d3: mixture
#>   impact:d1: normal
#>   impact:d2: normal
#>   impact:d3: pert
mc <- cm_monte_carlo(s, 400, seed = 1)
mc
#> <cm_mc> method: simultaneous
#>   Analysis: yield [%]
#>   Draws: 400, rejected: 0 (0.0%)
#>   Sampling: random
#>   Effective sample size: 400.0
summary(mc)
#>   disease       method     mean        sd       mcse    q0.025     q0.5
#> 1      d1 simultaneous 2.076567 0.5145080 0.02572540 1.1335619 2.097841
#> 2      d2 simultaneous 3.428872 1.2330534 0.06165267 0.9580435 3.407767
#> 3      d3 simultaneous 6.758133 0.8088373 0.04044186 5.1267784 6.763320
#>     q0.975 trimmed_mean    rel_mcse tail_share stability
#> 1 3.021695     2.074861 0.012388427 0.07431351        ok
#> 2 5.930501     3.411216 0.017980452 0.09468703        ok
#> 3 8.217512     6.758384 0.005984177 0.05810142        ok
summary(mc, what = "total")
#>         quantity       method     mean        sd       mcse   q0.025     q0.5
#> 1 adjusted_total simultaneous 2.071742 0.2101324 0.01050662 1.678115 2.059572
#> 2        raw_sum          raw 2.448621 0.2247769 0.01123884 2.004242 2.449640
#>     q0.975 trimmed_mean    rel_mcse tail_share stability
#> 1 2.491225     2.068209 0.005071395  0.1216242        ok
#> 2 2.864391     2.448254 0.004589868  0.1003012        ok
```

Draws with impossible inputs (e.g. a probability outside (0, 1), or
associations that no population can have together) are rejected and
counted by type with `summary(mc, what = "rejections")`. A draw whose
results are not finite is rejected as a whole. Report the rejection
rate: conditioning on acceptance changes the effective input
distribution.

## Several analyses on shared draws

[`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md)
builds one sampler per analysis of a
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
object. Each draw of the disease and association inputs is shared by all
analyses, with the same draw identifiers; each analysis draws its own
impacts:

``` r

pop <- cm_population(m$diseases, m$associations)
a <- cm_analyses(pop,
  yield = cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
  calving_interval = cm_impacts(c("d1", "d2", "d3"), c(4, 10, 2), units = "days")
)
bs <- cm_batch_sampler(a,
  diseases = list(d1 = dist_beta(20, 180)),
  associations = list("d1:d2" = dist_lognormal_ci(2, 1.4, 2.9)),
  impacts = list(yield = list(d1 = dist_normal(2.5, 0.5)),
                 calving_interval = list(d2 = dist_normal(10, 3)))
)
bs
#> <cm_batch_sampler> 2 analyses (yield, calving_interval); 2 shared population inputs
mcb <- cm_monte_carlo(bs, 200, seed = 1)
mcb
#> <cm_mc_batch> 200 draws, 2 analyses
#>   yield: 200 accepted, 0 rejected
#>   calving_interval: 200 accepted, 0 rejected
summary(mcb, diagnose = FALSE)[, c("analysis", "disease", "mean", "q0.025", "q0.975")]
#>           analysis disease      mean     q0.025    q0.975
#> 1            yield      d1 2.1669988  1.2546517  3.193909
#> 2            yield      d2 3.3854750  3.2277223  3.478043
#> 3            yield      d3 6.9344775  6.9188907  6.960504
#> 4 calving_interval      d1 2.9659706  1.7987533  3.602805
#> 5 calving_interval      d2 9.6606875  4.2743121 15.179476
#> 6 calving_interval      d3 0.3862423 -0.6117922  1.271536
```

Each analysis’s run is in `mcb$analyses`, and can be used with the
functions below.

## Comparing methods on the same draws

Give several methods to
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
and every accepted draw is adjusted with each of them (a draw is
rejected if any method fails on it).
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
then tabulates the results:

``` r

mc2 <- cm_monte_carlo(s, 400, method = c("published", "simultaneous"), seed = 1)
compare_methods(mc2)
#> <cm_comparison> methods: published, simultaneous
#> Units: %
#> Monte Carlo: 400 draws (0 rejected); statistic: mean
#> 
#> Adjusted values:
#>  disease raw_mean published simultaneous
#>       d1     2.45      2.01         2.08
#>       d2     4.98      3.75         3.43
#>       d3     7.30      6.57         6.76
#> 
#> Totals:
#>        quantity       method  mean  q0.5 trimmed_mean    mcse stability
#>  adjusted_total    published 2.076 2.069        2.072 0.01047        ok
#>  adjusted_total simultaneous 2.072 2.060        2.068 0.01051        ok
#>         raw_sum          raw 2.449 2.450        2.448 0.01124        ok
```

## Unstable estimates

Monte Carlo means can be unstable.
[`summary()`](https://rdrr.io/r/base/summary.html) checks every estimate
and, when one looks unstable, prints why and what to do. Here the raw
impact of `d2` is uncertain enough to be near zero or negative:

``` r

s_d2 <- cm_sampler(example_supplement(), impacts = list(d2 = dist_normal(0.5, 1.5)))
mc_d2 <- cm_monte_carlo(s_d2, 400, method = c("published", "simultaneous"), seed = 2)
sm <- summary(mc_d2)
#> 2 Monte Carlo estimate(s) may be unstable:
#> * d2 (published): the published approximation divides by m + c, which changes sign within the sampled inputs (in 34.2% of draws its sign differs from m), so the estimate has a pole inside the input distribution and its mean may not exist.
#>     Suggestion: report quantiles (e.g. the median), or use method = "simultaneous" (exact, no division). A trimmed mean is a different estimand. If the mean does not exist, more draws or importance sampling will not make it converge.
#> * d2 (simultaneous): the Monte Carlo standard error is 5.8% of the mean.
#>     Suggestion: increase n_draws to about 3,400 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> (See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)
sm[, c("disease", "method", "mean", "q0.5", "trimmed_mean", "mcse", "stability")]
#>   disease       method      mean       q0.5 trimmed_mean        mcse
#> 1      d1    published  2.466043  2.4622674    2.4630193 0.007941041
#> 2      d2    published -4.036952  0.4811604    0.4388501 3.764532112
#> 3      d3    published  7.438133  7.4397288    7.4357484 0.012726664
#> 4      d1 simultaneous  2.646521  2.6534324    2.6467940 0.008543511
#> 5      d2 simultaneous -1.391114 -1.4658329   -1.4088810 0.081114483
#> 6      d3 simultaneous  7.732377  7.7433377    7.7328094 0.013549669
#>       stability
#> 1            ok
#> 2 possible_pole
#> 3            ok
#> 4            ok
#> 5     imprecise
#> 6            ok
```

Each estimate gets one of five statuses (see
[`?cm_diagnose`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)):

- **`ok`**: no problem detected.
- **`imprecise`**: the Monte Carlo standard error is more than 5% of the
  mean. Use more draws, or Latin hypercube sampling.
- **`insufficient_info`**: the precision cannot be assessed, because
  fewer than two Latin hypercube blocks have draws with positive weight,
  or because the mean is (nearly) zero, so that a relative error is
  undefined. In the second case `mcse` still gives the absolute error,
  which can be judged against the size of effect that matters. A missing
  precision is never reported as `ok`.
- **`heavy_tail`**: the most extreme 1% of draws contribute more than
  60% of the variance. Importance sampling can help if they come from
  one region of one input (below).
- **`possible_pole`**: the published approximation, `m^2 / (m + c)`,
  divides by `m + c`, which takes both signs within the sampled inputs
  and in some draws has the opposite sign to `m`. The estimate then has
  a pole inside the input distribution, and its mean may not exist; no
  number of draws, and no sampling scheme, makes such a mean converge.
  Report quantiles (e.g. the median), or use the exact method. A trimmed
  mean is a different estimand. Removable cases are not flagged: when
  `c = 0` (e.g. a disease with no associated impacts) the formula
  reduces to `m`, and a sign change of `m` alone does not create a pole.

Draws that give non-finite results are rejected, and are noted as
`non_finite`.
[`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)
lists the flagged estimates with the suggestions:

``` r

cm_diagnose(mc_d2)
#> 2 Monte Carlo estimate(s) may be unstable:
#> * d2 (published): the published approximation divides by m + c, which changes sign within the sampled inputs (in 34.2% of draws its sign differs from m), so the estimate has a pole inside the input distribution and its mean may not exist.
#>     Suggestion: report quantiles (e.g. the median), or use method = "simultaneous" (exact, no division). A trimmed mean is a different estimand. If the mean does not exist, more draws or importance sampling will not make it converge.
#> * d2 (simultaneous): the Monte Carlo standard error is 5.8% of the mean.
#>     Suggestion: increase n_draws to about 3,400 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> (See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)
```

### Latin hypercube sampling

`sampling = "lhs"` stratifies each input’s distribution, which usually
reduces the Monte Carlo error of means. LHS draws are not independent,
so the draws are split into `lhs_replicates` independent blocks (10 by
default), and the standard error is estimated from the spread of the
block means. The pooled mean is a ratio estimator over the blocks (each
block’s weighted sum over the total weight), and its standard error uses
all blocks as sampled (`mc$n_blocks`), including any whose draws were
all rejected; with fewer than two blocks with positive weight the status
is `insufficient_info`:

``` r

mc_lhs <- cm_monte_carlo(s_d2, 400, method = "simultaneous", sampling = "lhs", seed = 2)
mc_lhs
#> <cm_mc> method: simultaneous
#>   Analysis: yield [%]
#>   Draws: 400, rejected: 0 (0.0%)
#>   Sampling: Latin hypercube, 10 replicate blocks
#>   Effective sample size: 400.0
summary(mc_lhs, diagnose = FALSE)[, c("disease", "mean", "mcse")]
#>   disease      mean         mcse
#> 1      d1  2.638968 0.0004894368
#> 2      d2 -1.319403 0.0046468501
#> 3      d3  7.720398 0.0007762274
sm[sm$method == "simultaneous", c("disease", "mean", "mcse")]
#>   disease      mean        mcse
#> 4      d1  2.646521 0.008543511
#> 5      d2 -1.391114 0.081114483
#> 6      d3  7.732377 0.013549669
```

### Importance sampling

When extreme values come from one region of one input,
[`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md)
builds a defensive mixture that samples that region more often. The
draws are then weighted by the ratio of the densities, so the estimates
still refer to the original input distributions:

``` r

prop <- cm_suggest_proposal(mc_d2, "d2", method = "simultaneous")
#> Proposal for impact:d2: 50% its own distribution, 50% uniform on [-3.972, 5.56], the range of impact:d2 in the 8 most extreme draws of d2 (simultaneous).
mc_is <- cm_monte_carlo(s_d2, 400, method = "simultaneous", proposal = prop, seed = 3)
mc_is
#> <cm_mc> method: simultaneous
#>   Analysis: yield [%]
#>   Draws: 400, rejected: 0 (0.0%)
#>   Sampling: random, importance sampling of impact:d2
#>   Effective sample size: 329.8
summary(mc_is, diagnose = FALSE)[, c("disease", "mean", "mcse")]
#>   disease      mean       mcse
#> 1      d1  2.638857 0.00727967
#> 2      d2 -1.318345 0.06911522
#> 3      d3  7.720221 0.01154527
```

Each proposal must cover the whole support of the input’s own
distribution, and point masses (fixed values) cannot be
importance-sampled. Both are checked before any draw is made. A proposal
that is too narrow is rejected:

``` r

tryCatch(cm_monte_carlo(s_d2, 10, proposal = list("impact:d2" = dist_uniform(-5, 5))),
         deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "The proposal for 'impact:d2' has support [-5, 5], which does not cover the input's support [-Inf, Inf]; use a defensive mixture that includes the input's own distribution (see cm_suggest_proposal())."
```

So is a mixture whose range spans the support but leaves a gap: draws
would never fall in the gap, and the weighted estimates would converge
to the mean over the covered part only. Here the raw impact of d3 has a
PERT distribution on \[5, 9\], and the proposal has no mass between 6.5
and 7.5:

``` r

gap <- dist_mixture(dist_uniform(4, 6.5), dist_uniform(7.5, 10))
tryCatch(cm_monte_carlo(s, 10, proposal = list("impact:d3" = gap)),
         deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "The proposal for 'impact:d3' has support [4, 6.5] u [7.5, 10], which does not cover the input's support [5, 9]; use a defensive mixture that includes the input's own distribution (see cm_suggest_proposal())."
```

Including the input’s own distribution as a mixture component, as
[`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md)
does, always covers its support.

A defensive mixture guarantees support, but not a finite variance or
better precision: compare the standard errors. Importance sampling
reduces the error of a mean that exists. It cannot fix a `possible_pole`
estimate.

## Scenarios by reweighting

[`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
replaces input distributions and reweights the existing draws
(importance sampling), so no re-run is needed. The sampler above used a
defensive mixture for the d2:d3 odds ratio, so a scenario with a
stronger association still has support:

``` r

sc <- cm_scenario(mc, list("assoc:d2:d3" = dist_lognormal_ci(4, 2.5, 6.4)))
sc$ess
#> [1] 135.8104
summary(sc, what = "total", diagnose = FALSE)
#>         quantity       method     mean        sd       mcse   q0.025     q0.5
#> 1 adjusted_total simultaneous 1.990148 0.1805131 0.01443515 1.665442 1.994400
#> 2        raw_sum          raw 2.447616 0.2184128 0.01792253 2.065892 2.436139
#>     q0.975 trimmed_mean    rel_mcse tail_share stability
#> 1 2.340430     1.991178 0.007253305  0.1512718        ok
#> 2 2.862396     2.449011 0.007322444  0.1687684        ok
```

Check the effective sample size: a small value means the scenario is
poorly covered by the original draws. A scenario distribution must lie
within the range the input was sampled from, which is checked.
Reweighting cannot recover rejected draws.

[`cm_reweight()`](https://rasmussenphilip.github.io/deconflate/reference/cm_reweight.md)
takes any log density ratio as a function of the sampled inputs
(`mc$params`). For example, conditioning on a raw impact of d1 above 2%:

``` r

cond <- cm_reweight(mc, function(p) ifelse(p[["impact:d1"]] > 2, 0, -Inf))
cond$ess
#> [1] 319
summary(cond, diagnose = FALSE)[, c("disease", "mean", "q0.5")]
#>   disease     mean     q0.5
#> 1      d1 2.246454 2.227384
#> 2      d2 3.471261 3.508868
#> 3      d3 6.765979 6.793630
```

## Productivity gaps over the draws

[`cm_mc_gap()`](https://rasmussenphilip.github.io/deconflate/reference/cm_mc_gap.md)
applies
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
to every accepted draw. The observed mean and the unit value can be
fixed or drawn:

``` r

g <- cm_mc_gap(mc, observed = 10000, direction = "decrease", effect = "percent",
               unit_value = dist_uniform(0.25, 0.35), seed = 5)
g$summary
#>         method gap_mean gap_q0.025 gap_q0.975 value_mean value_q0.025
#> 1 simultaneous 211.6042   170.6756   255.4873   63.43205     46.30161
#>   value_q0.975
#> 1     82.77978
```

## Sensitivity

One-at-a-time sensitivity (as in Rasmussen et al. 2024, Fig. 7):

``` r

oat <- sensitivity_oat(example_supplement(), variation = 0.2)
oat
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

Screening associations, including pairs with no estimate, to find those
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
#> <cm_screen> 6 scenarios; baseline total 2.109
#>   pair                status scenario total  change rel_change max_rank_shift
#>  d2:d3             specified OR x 0.5  2.31  0.1974     0.0936              0
#>  d2:d3             specified   OR x 2  1.94 -0.1643    -0.0779              0
#>  d1:d3 independent (default)   OR = 2  2.01 -0.0963    -0.0457              0
#>  d1:d3 independent (default) OR = 0.5  2.20  0.0859     0.0407              0
#>  d1:d2             specified   OR x 2  2.05 -0.0582    -0.0276              0
#>  d1:d2             specified OR x 0.5  2.16  0.0549     0.0260              0
#>  rank_corr failed
#>          1   <NA>
#>          1   <NA>
#>          1   <NA>
#>          1   <NA>
#>          1   <NA>
#>          1   <NA>
```

The screens report scenarios that cannot be run, with the reason in
column `failed`, rather than skipping them. Here, weakening any of three
strong associations makes the pairs jointly infeasible:

``` r

strong <- cm_model(
  cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
  cm_impacts(c("a", "b", "c"), c(1, 2, 3)),
  associations = cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 20))
)
screen_associations(strong, multipliers = c(0.05, 2))
#> <cm_screen> 6 scenarios; baseline total 1.322
#>  pair    status  scenario total  change rel_change max_rank_shift rank_corr
#>   b:c specified    OR x 2  1.24 -0.0824    -0.0624              0         1
#>   a:c specified    OR x 2  1.29 -0.0323    -0.0244              0         1
#>   a:b specified    OR x 2  1.34  0.0178     0.0135              0         1
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

[`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
and
[`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md)
(see
[`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md))
work in the same way. Comparing complete scenarios:

``` r

base <- example_supplement()
compare_scenarios(base = base, d1_d3_linked = set_association(base, "d1", "d3", 2))$totals
#>       scenario    total failed rel_to_first
#> 1         base 2.109229   <NA>   0.00000000
#> 2 d1_d3_linked 2.012881   <NA>  -0.04567916
```
