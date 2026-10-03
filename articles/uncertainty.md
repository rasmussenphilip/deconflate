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
attaches distributions to a model. Each draw is a complete model, so
disease probabilities and associations are shared by all outcomes within
a draw. `outcome_correlation` additionally correlates a disease’s
impacts across outcomes (Gaussian copula):

``` r

m <- example_uk_dairy_2022()
R <- matrix(c(1, 0.5, 0.5, 1), 2,
            dimnames = list(c("yield", "fertility"), c("yield", "fertility")))
s <- cm_sampler(
  m,
  diseases = list(LAM = dist_beta(30, 70), MAS = dist_beta(30, 70)),
  associations = list("LAM:PTB" = dist_mixture(dist_lognormal_ci(2.7, 1.5, 4.9),
                                               dist_lognormal(log(2.7), 0.8),
                                               weights = c(0.8, 0.2)),
                      "LAM:SCK" = dist_lognormal_ci(2.01, 1.62, 2.49)),
  impacts = list("yield:LAM" = dist_normal(5.54, 0.8),
                 "fertility:LAM" = dist_normal(12.47, 2)),
  outcome_correlation = R
)
s
#> <cm_sampler> 6 uncertain inputs (2 drawn with outcome correlation)
#>   prob:LAM: beta
#>   prob:MAS: beta
#>   assoc:LAM:PTB: mixture
#>   assoc:LAM:SCK: lognormal
#>   impact:yield:LAM: normal
#>   impact:fertility:LAM: normal
mc <- cm_monte_carlo(s, 200, method = "simultaneous",
                     economics = uk_dairy_2022_economics(), seed = 1)
mc
#> <cm_mc> method: simultaneous
#>   Draws: 200, rejected as infeasible: 0 (0.0%)
#>   Sampling: random
#>   Effective sample size: 200.0
summary(mc, what = "total")
#>     outcome disease       method      mean       sd      mcse    q0.025
#> 1   culling     all simultaneous  58.78375  3.45915 0.2445989  51.54509
#> 2 fertility     all simultaneous  95.44686 11.02336 0.7794689  76.09317
#> 3     yield     all simultaneous 168.20789 11.23947 0.7947506 146.12134
#> 4     total     all simultaneous 322.43849 23.08920 1.6326528 275.35343
#>        q0.5    q0.975 trimmed_mean    rel_mcse tail_share stability
#> 1  58.74946  64.82462     58.76218 0.004160995 0.09153850        ok
#> 2  94.12680 120.57118     94.98735 0.008166523 0.07612709        ok
#> 3 168.49768 189.95906    168.03140 0.004724811 0.08047083        ok
#> 4 322.75810 372.44177    321.92942 0.005063455 0.06821893        ok
```

Draws with impossible inputs (e.g. a probability outside (0, 1), or an
association incompatible with the sampled probabilities) are rejected
and counted. Report the rejection rate: conditioning on feasibility
changes the effective input distribution.

## Comparing methods on the same draws

Give several methods to
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
and every accepted draw is adjusted with each of them.
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
then tabulates the results:

``` r

mc2 <- cm_monte_carlo(s, 200, method = c("published", "simultaneous"), seed = 1)
compare_methods(mc2)
#> <cm_comparison> methods: published, simultaneous
#> Monte Carlo: 200 draws (0 rejected); statistic: mean
#> 
#> Adjusted impacts:
#>    outcome disease unit raw_mean published simultaneous
#>      yield      CO    %     0.00     0.000    -6.04e-01
#>      yield      DA    %     4.04     2.370     2.20e+00
#>      yield     DYS    %     4.05     2.940     2.95e+00
#>      yield     FAS    %     7.33     7.330     7.33e+00
#>      yield     GIN    %     3.28     3.280     3.28e+00
#>      yield     LAM    %     5.49     4.800     5.00e+00
#>      yield     MAS    %     4.57     3.760     4.08e+00
#>      yield     MET    %     3.95     2.410     1.87e+00
#>      yield      MF    %     0.41     0.100    -5.13e-01
#>      yield     NEO    %     4.20     4.200     4.20e+00
#>      yield     PTB    %     5.90     4.490     4.20e+00
#>      yield      RP    %     7.38     6.100     6.58e+00
#>      yield     SCK    %     3.05     1.930     1.67e+00
#>  fertility      CO    %    11.30    11.100     1.14e+01
#>  fertility      DA    %     0.00     0.000    -4.08e-01
#>  fertility     DYS    %     6.96     6.020     6.08e+00
#>  fertility     FAS    %     0.00     0.000    -1.23e-16
#>  fertility     GIN    %     1.20     1.200     1.20e+00
#>  fertility     LAM    %    12.50    11.900     1.24e+01
#>  fertility     MAS    %     0.00     0.000    -9.48e-01
#>  fertility     MET    %     4.74     4.100     4.61e+00
#>  fertility      MF    %     0.00     0.000     2.88e-01
#>  fertility     NEO    %     7.21     7.210     7.21e+00
#>  fertility     PTB    %     5.79     4.000     3.14e+00
#>  fertility      RP    %     2.74     1.680     1.17e+00
#>  fertility     SCK    %     1.50     0.563    -9.18e-01
#>    culling      CO    %     0.00     0.000    -2.67e+00
#>    culling      DA    %    31.40    21.900     2.25e+01
#>    culling     DYS    %    14.20    11.900     1.28e+01
#>    culling     FAS    %     0.00     0.000    -7.58e-16
#>    culling     GIN    %     0.00     0.000    -5.90e-16
#>    culling     LAM    %    25.60    22.600     2.38e+01
#>    culling     MAS    %    21.30    16.900     1.78e+01
#>    culling     MET    %    17.40    12.300     1.26e+01
#>    culling      MF    %    20.60    15.600     1.60e+01
#>    culling     NEO    %     9.89     9.890     9.89e+00
#>    culling     PTB    %    19.60    13.600     1.17e+01
#>    culling      RP    %     0.00     0.000    -5.18e+00
#>    culling     SCK    %    15.70     9.710     7.70e+00
```

## Unstable estimates

Monte Carlo means can be unstable.
[`summary()`](https://rdrr.io/r/base/summary.html) checks every estimate
and, when one looks unstable, prints why and what to do. Here the raw
yield impact of `d2` is uncertain enough to be near zero or negative:

``` r

s_d2 <- cm_sampler(example_supplement(),
                   impacts = list("yield:d2" = dist_normal(0.5, 1.5)))
mc_d2 <- cm_monte_carlo(s_d2, 400, method = c("published", "simultaneous"), seed = 2)
sm <- summary(mc_d2)
#> 2 Monte Carlo estimate(s) may be unstable:
#> * yield / d2 (published): the published approximation divides by m + c, which changes sign in 8.0% of draws, so the mean does not exist.
#>     Suggestion: report the median or trimmed_mean, or use method = "simultaneous" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge.
#> * yield / d2 (simultaneous): the Monte Carlo standard error is 5.8% of the mean.
#>     Suggestion: increase n_draws to about 3,400 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> (See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)
sm[, c("disease", "method", "mean", "q0.5", "trimmed_mean", "stability")]
#>   disease       method        mean         q0.5 trimmed_mean stability
#> 1      d1    published  0.02466043  0.024622674  0.024630193        ok
#> 2      d2    published -0.04036952  0.004811604  0.004388501   no_mean
#> 3      d3    published  0.07438133  0.074397288  0.074357484        ok
#> 4      d1 simultaneous  0.02646521  0.026534324  0.026467940        ok
#> 5      d2 simultaneous -0.01391114 -0.014658329 -0.014088810 imprecise
#> 6      d3 simultaneous  0.07732377  0.077433377  0.077328094        ok
```

There are three kinds of instability (see
[`?cm_diagnose`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)):

- **`no_mean`**: the published approximation, `m^2 / (m + c)`, divides
  by a quantity that changes sign across draws. Its distribution then
  has no mean, so no number of draws, and no sampling scheme, makes the
  mean converge. Report the median or trimmed mean, or use the exact
  method.
- **`heavy_tail`**: a few draws dominate the variance. Importance
  sampling can help (below).
- **`imprecise`**: the Monte Carlo error is large relative to the mean.
  Use more draws, or Latin hypercube sampling.

[`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)
lists the flagged estimates with the suggestions:

``` r

cm_diagnose(mc_d2)
#> 2 Monte Carlo estimate(s) may be unstable:
#> * yield / d2 (published): the published approximation divides by m + c, which changes sign in 8.0% of draws, so the mean does not exist.
#>     Suggestion: report the median or trimmed_mean, or use method = "simultaneous" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge.
#> * yield / d2 (simultaneous): the Monte Carlo standard error is 5.8% of the mean.
#>     Suggestion: increase n_draws to about 3,400 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> (See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)
```

### Latin hypercube sampling

`sampling = "lhs"` stratifies each input’s distribution, which usually
reduces the Monte Carlo error of means:

``` r

mc_lhs <- cm_monte_carlo(s_d2, 400, method = "simultaneous", sampling = "lhs", seed = 2)
summary(mc_lhs, diagnose = FALSE)[, c("disease", "mean", "mcse")]
#>   disease        mean         mcse
#> 1      d1  0.02638348 8.241601e-05
#> 2      d2 -0.01313510 7.824806e-04
#> 3      d3  0.07719414 1.307085e-04
```

### Importance sampling

When extreme values come from one region of one input,
[`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md)
builds a defensive mixture that samples that region more often. The
draws are then weighted by the ratio of the densities, so the estimates
still refer to the original input distributions:

``` r

prop <- cm_suggest_proposal(mc_d2, "yield", "d2", method = "simultaneous")
#> Proposal for impact:yield:d2: 50% its own distribution, 50% uniform on [-3.972, 5.56], the range of impact:yield:d2 in the 8 most extreme draws of yield / d2 (simultaneous).
mc_is <- cm_monte_carlo(s_d2, 400, method = "simultaneous", proposal = prop, seed = 3)
mc_is
#> <cm_mc> method: simultaneous
#>   Draws: 400, rejected as infeasible: 0 (0.0%)
#>   Sampling: random, importance sampling of impact:yield:d2
#>   Effective sample size: 329.8
summary(mc_is, diagnose = FALSE)[, c("disease", "mean", "mcse")]
#>   disease        mean         mcse
#> 1      d1  0.02638857 8.872334e-05
#> 2      d2 -0.01318345 8.423642e-04
#> 3      d3  0.07720221 1.407117e-04
```

Importance sampling reduces the error of a mean that exists. It cannot
fix a `no_mean` estimate.

## Scenarios by reweighting

[`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
replaces input distributions and reweights the existing draws
(importance sampling), so no re-run is needed. Here the sampler used a
defensive mixture for the lameness-paratuberculosis odds ratio, so a
scenario with a stronger association still has support:

``` r

sc <- cm_scenario(mc, list("assoc:LAM:PTB" = dist_lognormal_ci(4, 2.5, 6.4)))
sc$ess
#> [1] 71.16948
summary(sc, what = "total")
#>     outcome disease       method      mean        sd     mcse    q0.025
#> 1   culling     all simultaneous  58.51338  3.639799 0.431450  50.74542
#> 2 fertility     all simultaneous  94.67543 11.854718 1.405220  76.78881
#> 3     yield     all simultaneous 167.47339 10.797443 1.279894 148.26144
#> 4     total     all simultaneous 320.66220 23.084688 2.736384 272.15182
#>        q0.5    q0.975 trimmed_mean    rel_mcse tail_share stability
#> 1  58.33510  65.74893     58.56282 0.007373527  0.1790683        ok
#> 2  92.33746 124.24201     93.98920 0.014842500  0.2615265        ok
#> 3 167.60112 190.77552    167.57575 0.007642373  0.1952306        ok
#> 4 321.06635 379.84216    320.24870 0.008533542  0.2348487        ok
```

Check the effective sample size: a small value means the scenario is
poorly covered by the original draws. Inputs drawn with an outcome
correlation cannot be reweighted individually.

## Sensitivity

One-at-a-time sensitivity (as in Rasmussen et al. 2024, Fig. 7):

``` r

oat <- sensitivity_oat(example_supplement(), variation = 0.2)
oat
#>             input value  total_low total_high        swing  rel_swing
#> 1 impact:yield:d3  7.50 0.01844070 0.02374387 0.0053031715 0.25142703
#> 2         prob:d3  0.20 0.01851628 0.02374637 0.0052300893 0.24796215
#> 3 impact:yield:d2  5.00 0.01998423 0.02220035 0.0022161264 0.10506809
#> 4         prob:d2  0.15 0.02026335 0.02197276 0.0017094138 0.08104449
#> 5     assoc:d2:d3  3.00 0.02170224 0.02062015 0.0010820854 0.05130242
#> 6 impact:yield:d1  2.50 0.02063348 0.02155110 0.0009176175 0.04350488
#> 7         prob:d1  0.10 0.02069218 0.02149697 0.0008047878 0.03815555
#> 8     assoc:d1:d3  1.00 0.02138876 0.02084124 0.0005475254 0.02595856
#> 9     assoc:d1:d2  2.00 0.02127728 0.02093846 0.0003388228 0.01606382
```

Screening associations, including pairs with no estimate, to find those
worth estimating:

``` r

screen_associations(example_supplement(), or_values = c(0.5, 2))
#> <cm_screen> 6 scenarios; baseline total 0.02109
#>   pair    status scenario  total    change rel_change max_rank_shift rank_corr
#>  d2:d3 specified OR x 0.5 0.0231  0.001974     0.0936              0         1
#>  d2:d3 specified   OR x 2 0.0194 -0.001643    -0.0779              0         1
#>  d1:d3 specified   OR x 2 0.0201 -0.000963    -0.0457              0         1
#>  d1:d3 specified OR x 0.5 0.0220  0.000859     0.0407              0         1
#>  d1:d2 specified   OR x 2 0.0205 -0.000582    -0.0276              0         1
#>  d1:d2 specified OR x 0.5 0.0216  0.000549     0.0260              0         1
```

Comparing complete scenarios:

``` r

base <- example_supplement()
compare_scenarios(base = base, d1_d3_linked = set_association(base, "d1", "d3", 2))$totals
#>       scenario      total rel_to_first
#> 1         base 0.02109229   0.00000000
#> 2 d1_d3_linked 0.02012881  -0.04567916
```
