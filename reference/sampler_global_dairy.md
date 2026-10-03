# Monte Carlo sampler for the global dairy inputs (Rasmussen et al. 2024)

A
[`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md)
for
[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md):
shared draws of the odds ratios (and, with `inputs = "tables"`, the
incidences), and each analysis's impact distributions. PERT
distributions use the reported central value as the mode (shape 4), and
normal distributions of odds ratios are truncated at zero.

## Usage

``` r
sampler_global_dairy(inputs = c("analysis", "tables"))
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`.

## Value

A
[`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md).

## Details

With `inputs = "analysis"` the disease probabilities are fixed and the
impact distributions use the unrounded parameters of the analysis code.
With `inputs = "tables"`, incidences are drawn from the global
distributions of Table 2 and Tables 3-4 are used as printed. (The 2024
culling analysis, HR - 1 adjusted as an additive impact, is sampled only
inside
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).)

With `method = "published"`, Monte Carlo means reproduce Table 5 (see
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
and
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).

## Examples

``` r
# \donttest{
mc <- cm_monte_carlo(sampler_global_dairy(), 200, method = "published", seed = 1)
s <- summary(mc, diagnose = FALSE)
s[s$analysis == "yield", c("disease", "mean")]
#>    disease       mean
#> 1       CK 0.03544741
#> 2       CM 1.40230252
#> 3       DA 1.06739016
#> 4      DYS 3.49063187
#> 5      LAM 2.51461430
#> 6      MET 2.88009853
#> 7       MF 0.07078056
#> 8       OC 2.54091312
#> 9      PTB 3.23721070
#> 10      RP 2.30415060
#> 11     SCK 7.09380108
#> 12     SCM 5.61896630
# }
```
