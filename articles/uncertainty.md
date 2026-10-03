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
summary(mc, what = "total")
#>     outcome disease      mean       sd      mcse    q0.025      q0.5    q0.975
#> 1   culling     all  58.78375  3.45915 0.2445989  51.54509  58.74946  64.82462
#> 2 fertility     all  95.44686 11.02336 0.7794689  76.09317  94.12680 120.57118
#> 3     total     all 322.43849 23.08920 1.6326528 275.35343 322.75810 372.44177
#> 4     yield     all 168.20789 11.23947 0.7947506 146.12134 168.49768 189.95906
```

Draws with impossible inputs (e.g. a probability outside (0, 1), or an
association incompatible with the sampled probabilities) are rejected
and counted. Report the rejection rate: conditioning on feasibility
changes the effective input distribution.

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
#>     outcome disease      mean        sd     mcse    q0.025      q0.5    q0.975
#> 1   culling     all  58.51338  3.639799 0.431450  50.74542  58.33510  65.74893
#> 2 fertility     all  94.67543 11.854718 1.405220  76.78881  92.33746 124.24201
#> 3     total     all 320.66220 23.084688 2.736384 272.15182 321.06635 379.84216
#> 4     yield     all 167.47339 10.797443 1.279894 148.26144 167.60112 190.77552
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
