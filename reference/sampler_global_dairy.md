# Monte Carlo sampler for the global dairy inputs (Rasmussen et al. 2024)

Input distributions from Tables 2-4 of Rasmussen et al. (2024) for
[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md):
global beta distributions for incidence, PERT distributions (the
reported central value is used as the mode), normal distributions for
odds ratios (truncated at zero) and impacts.

## Usage

``` r
sampler_global_dairy()
```

## Value

A
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

## Details

Running
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
with `method = "published"` and averaging the adjusted yield impacts
over draws reproduces Table 5 closely. Fertility impacts whose raw
distributions extend below zero (e.g. metritis, subclinical ketosis,
subclinical mastitis) are unstable under the published approximation,
and their means are sensitive to how such draws are handled.

## Examples

``` r
# \donttest{
mc <- cm_monte_carlo(sampler_global_dairy(), 200, method = "published", seed = 1)
s <- summary(mc)
s[s$outcome == "yield", c("disease", "mean")]
#>    disease         mean
#> 2       CK 0.0003354865
#> 4       CM 0.0134195678
#> 6       DA 0.0124002061
#> 8      DYS 0.0347643531
#> 10     LAM 0.0274527187
#> 12     MET 0.0288749662
#> 14      MF 0.0006987672
#> 16      OC 0.0255543076
#> 18     PTB 0.0330164657
#> 20      RP 0.0230020159
#> 22     SCK 0.0706897214
#> 24     SCM 0.0558738798
# }
```
