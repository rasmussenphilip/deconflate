# Monte Carlo sampler for the global dairy inputs (Rasmussen et al. 2024)

Input distributions for
[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md).
PERT distributions use the reported central value as the mode (shape 4),
and normal distributions of odds ratios are truncated at zero.

## Usage

``` r
sampler_global_dairy(
  inputs = c("analysis", "tables"),
  culling = TRUE,
  culling_scale = c("excess_hr", "hazard_ratio")
)
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`; see the section above.

- culling:

  Include the culling outcome?

- culling_scale:

  `"excess_hr"` (HR - 1) or `"hazard_ratio"`, as in
  [`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md).
  Hazard-ratio distributions are the HR - 1 distributions shifted by 1.

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
#> 8 Monte Carlo estimate(s) may be unstable:
#> * yield / DA (published): the Monte Carlo standard error is 5.9% of the mean.
#>     Suggestion: increase n_draws to about 1,760 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> * fertility / DA (published): the published approximation divides by m + c, which changes sign in 0.5% of draws, so the mean does not exist.
#>     Suggestion: report the median or trimmed_mean, or use method = "simultaneous" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge.
#> * culling / DYS (published): the Monte Carlo standard error is 6.7% of the mean.
#>     Suggestion: increase n_draws to about 2,222 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> * fertility / MET (published): the published approximation divides by m + c, which changes sign in 0.5% of draws, so the mean does not exist.
#>     Suggestion: report the median or trimmed_mean, or use method = "simultaneous" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge.
#> * culling / MET (published): the Monte Carlo standard error is 8.1% of the mean.
#>     Suggestion: increase n_draws to about 3,317 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> * fertility / PTB (published): the published approximation divides by m + c, which changes sign in 0.5% of draws, so the mean does not exist.
#>     Suggestion: report the median or trimmed_mean, or use method = "simultaneous" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge.
#> * fertility / SCK (published): the published approximation divides by m + c, which changes sign in 1.5% of draws, so the mean does not exist.
#>     Suggestion: report the median or trimmed_mean, or use method = "simultaneous" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge.
#> * fertility / SCM (published): the Monte Carlo standard error is 8.1% of the mean.
#>     Suggestion: increase n_draws to about 3,321 for a 2% standard error, or use sampling = "lhs" (Latin hypercube).
#> (See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)
s[s$outcome == "yield", c("disease", "mean")]
#>    disease         mean
#> 1       CK 0.0003629443
#> 4       CM 0.0139387593
#> 7       DA 0.0123009108
#> 10     DYS 0.0341056713
#> 13     LAM 0.0265759525
#> 16     MET 0.0296544668
#> 19      MF 0.0006981994
#> 22      OC 0.0247995885
#> 25     PTB 0.0319448378
#> 28      RP 0.0233288844
#> 31     SCK 0.0711448636
#> 34     SCM 0.0565964340
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
