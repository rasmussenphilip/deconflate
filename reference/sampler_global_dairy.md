# Monte Carlo sampler for the global dairy inputs (Rasmussen et al. 2024)

Input distributions for
[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md).
PERT distributions use the reported central value as the mode (shape 4),
and normal distributions of odds ratios are truncated at zero.

## Usage

``` r
sampler_global_dairy(inputs = c("analysis", "tables"), culling = TRUE)
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`; see the section above.

- culling:

  Include the culling outcome?

## Value

A
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

## Two versions of the inputs

- `inputs = "analysis"` (default) follows the published analysis code
  (1st revision): disease probabilities are fixed (not drawn), and the
  impact distributions use the unrounded parameters. For culling, the
  analysis entered HR - 1 and scaled the standard deviations of normal
  distributions by (HR - 1) / HR, which narrows them (e.g. clinical
  ketosis: SD 0.10 instead of 0.30). This is kept for reproduction.

- `inputs = "tables"` draws incidence from the global distributions of
  Table 2 and uses Tables 3-4 as printed. Culling hazard ratios are
  shifted by 1 with their standard deviations unchanged.

## Reproducing Table 5

With `inputs = "analysis"` and `method = "published"`, Monte Carlo means
of the adjusted impacts reproduce Table 5 for yield (within about 0.02
points, except ovarian cyst 2.51 vs 2.59 and paratuberculosis 3.32 vs
3.37) and culling (within about 0.04, after adding 1). The analysis code
used negative odds-ratio draws as they were; the package truncates them
at zero, which accounts for the paratuberculosis difference. Fertility
means for impacts whose raw distributions extend below zero (displaced
abomasum, metritis, subclinical ketosis, subclinical mastitis) are
unstable under the published approximation. See
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md).

## Examples

``` r
# \donttest{
mc <- cm_monte_carlo(sampler_global_dairy(), 200, method = "published", seed = 1)
s <- summary(mc)
s[s$outcome == "yield", c("disease", "mean")]
#>    disease         mean
#> 3       CK 0.0003629443
#> 6       CM 0.0139387593
#> 9       DA 0.0123009108
#> 12     DYS 0.0341056713
#> 15     LAM 0.0265759525
#> 18     MET 0.0296544668
#> 21      MF 0.0006981994
#> 24      OC 0.0247995885
#> 27     PTB 0.0319448378
#> 30      RP 0.0233288844
#> 33     SCK 0.0711448636
#> 36     SCM 0.0565964340
# Culling: adjusted hazard ratio = adjusted (HR - 1) + 1
cull <- s[s$outcome == "culling", ]
data.frame(disease = cull$disease, hr_adjusted = 1 + cull$mean)
#>    disease hr_adjusted
#> 1       CK    1.182933
#> 2       CM    1.887017
#> 3       DA    2.874214
#> 4      DYS    1.184845
#> 5      LAM    1.391928
#> 6      MET    1.028400
#> 7       MF    2.638197
#> 8       OC    1.474289
#> 9      PTB    2.062104
#> 10      RP    1.283543
#> 11     SCK    1.659755
#> 12     SCM    1.262059
# }
```
