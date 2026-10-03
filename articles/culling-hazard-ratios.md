# Culling and hazard ratios

``` r

library(deconflate)
```

Culling (and mortality) impacts are usually reported as hazard ratios
(HRs) from survival models. A raw HR for one disease is conflated with
the HRs of associated diseases, just like a yield impact, but HRs
combine multiplicatively, and the culling they cause has to be counted
without counting a cow twice. This vignette compares the ways the
package can handle them.

## Options

| Approach | What is adjusted | How diseases combine | How to use it |
|----|----|----|----|
| Multiplicative, exact | log HR, over the distribution of disease combinations | HRs multiply | `scale = "hazard_ratio"`, `method = "global"` |
| Multiplicative, first-order | log HR, pairwise | HRs multiply (approximately) | `scale = "hazard_ratio"`, `method = "simultaneous"` |
| HR - 1 (Rasmussen et al. 2024) | HR - 1, with eq. 16 | excess HRs add | `scale = "hazard_ratio"`, `method = "published"` |
| Excess risk (Rasmussen et al. 2022) | excess culling risk | excess risks add | [`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md), [`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md) |

Under the multiplicative model, a cow’s hazard is
`h0 * exp(sum_i beta_i * D_i)`. The raw HR of disease `i` is the ratio
of the average hazard among cows with and without `i`. The global method
solves for the `beta`s exactly, so that the adjusted HRs reproduce the
raw HRs over the fitted distribution of disease combinations. The
simultaneous method uses the first-order version, which needs only the
pairwise tables.

If a source HR comes from a model that already included the associated
diseases, set `adjusted_for` in
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
so that it is not adjusted again (with the simultaneous method).

## The global dairy example

`example_global_dairy(culling_scale = "hazard_ratio")` enters the
culling HRs of Rasmussen et al. (2024) as hazard ratios:

``` r

m <- example_global_dairy(culling_scale = "hazard_ratio")
cmp <- compare_methods(m)
cu <- cmp$impacts[cmp$impacts$outcome == "culling", ]
cu[, c("disease", "raw", "published", "simultaneous", "global")]
#>    disease      raw published simultaneous    global
#> 25      CK 1.500100  1.177580    0.9619070 0.9821306
#> 26      CM 2.300000  1.903941    1.8350864 1.8193536
#> 27      DA 2.851179  2.197930    1.8250542 1.8476664
#> 28     DYS 1.258143  1.098377    1.1165198 1.0505432
#> 29     LAM 1.744976  1.380683    1.2864923 1.2230210
#> 30     MET 1.116444  1.012411    0.7369879 0.7198427
#> 31      MF 2.999886  2.647637    2.6233669 2.6601486
#> 32      OC 1.620000  1.458644    1.5487477 1.5579660
#> 33     PTB 2.310508  2.047235    2.0203082 2.0165455
#> 34      RP 1.599928  1.284496    1.2305974 1.1816995
#> 35     SCK 1.920000  1.675253    1.7412069 1.7107120
#> 36     SCM 1.449996  1.254928    1.2520842 1.2109005
```

The published column reproduces the paper’s adjustment at the central
values (Table 5 reports Monte Carlo means). The methods agree within
about 0.1 for most diseases, and differ most where comorbidity is
heaviest (clinical ketosis, lameness, displaced abomasum, metritis).

## Culling attributable to disease

[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
turns adjusted HRs into the part of the overall culling risk that
disease causes. It finds the baseline hazard for which the population
culling risk equals the observed rate, and compares it with the risk of
a disease-free cow. A cow’s risk cannot exceed 1, so a cow with several
diseases is counted once. The attributable risk is then allocated to
diseases with Shapley values.

With the global average replacement rate of 23.66% (Rasmussen et
al. 2024, Table 1):

``` r

glob <- cmp$results$global
ar <- attributable_risk(glob, overall_risk = 0.2366)
ar
#> <cm_attributable> outcome: culling
#>   Overall risk 0.2366; disease-free risk 0.12; attributable 0.1166 (49.3% of the overall risk)
#> 
#>  disease hr_adjusted attributable      share
#>       CK      0.9821   -9.665e-05 -0.0008292
#>       CM      1.8194    2.981e-02  0.2557165
#>       DA      1.8477    2.842e-03  0.0243839
#>      DYS      1.0505    4.714e-04  0.0040443
#>      LAM      1.2230    7.989e-03  0.0685350
#>      MET      0.7198   -4.730e-03 -0.0405729
#>       MF      2.6601    5.123e-03  0.0439501
#>       OC      1.5580    8.424e-03  0.0722667
#>      PTB      2.0165    1.455e-02  0.1247869
#>       RP      1.1817    3.403e-03  0.0291945
#>      SCK      1.7107    3.590e-02  0.3079776
#>      SCM      1.2109    1.289e-02  0.1105465
#> 
#> Combinations skipped by `max_present`: probability 2.76e-06
```

The 2024 paper instead converted each adjusted HR to an excess risk
relative to the overall rate, `HR * r / (HR * r + 1 - r) - r`, and
summed the excess risks weighted by prevalence. That sum counts cows
with several diseases more than once:

``` r

r <- 0.2366
P <- m$diseases$prob
hr_2024 <- cmp$impacts$published[cmp$impacts$outcome == "culling"]
sum(P * (hr_2024 * r / (hr_2024 * r + 1 - r) - r))   # paper's approach
#> [1] 0.1475373
ar$summary$attributable                               # model-based
#> [1] 0.1165722
```

The paper’s approach attributes about 14.8 percentage points of the
23.7% culling rate to disease; the multiplicative model attributes about
11.7.

## Valuing culling with the other outcomes

In [`summary()`](https://rdrr.io/r/base/summary.html),
[`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)
and
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md),
a hazard-ratio outcome is valued with
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md):
its `observed` value is the overall risk (a proportion) and its
`unit_value` the value of a cow removed.

``` r

eco <- list(observed = c(yield = 5013, culling = 0.2366),
            unit_value = c(yield = 0.5981, culling = 1299.33 - 785.86))
summary(glob, economics = eco)$totals
#>     outcome   raw_loss adjusted_loss reduction  observed disease_free
#> 1     yield 0.09931174    0.07494192 0.2453871 5013.0000 5419.1192214
#> 2   culling         NA            NA        NA    0.2366    0.1200278
#> 3 fertility 0.07471713    0.03933738 0.4735159        NA           NA
#>           gap     value
#> 1 406.1192214 242.89991
#> 2   0.1165722  59.85632
#> 3          NA        NA
```

(Global averages from Rasmussen et al. 2024, Table 1: milk yield, milk
price per kg, and replacement price less culled-cow price.)

## Caveats

- The models assume a constant hazard within the period and that the HRs
  apply for the whole period. Most source HRs come from models in which
  disease status changes during the lactation.
- Hazard ratios are not collapsible: even without confounding, an
  average HR differs from the HR within subgroups. The exact
  multiplicative model accounts for this given the joint distribution;
  the first-order method and the HR - 1 approach do not.
- Culling and death compete. Make sure the source HRs are of the same
  kind (cause-specific or subdistribution).
