# Global dairy culling hazard ratios (Rasmussen et al. 2024)

The culling hazard ratios of Rasmussen et al. (2024), Table 4 (or the
analysis inputs; see
[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md))
as a hazard-ratio model for
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
and
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).
With `inputs = "analysis"`, metritis uses the analysis value (PERT mode
1.116) rather than Table 4 (normal, mean 1.05).

## Usage

``` r
example_global_dairy_hr(inputs = c("analysis", "tables"))
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`.

## Value

A
[`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md).

## Examples

``` r
hr <- example_global_dairy_hr()
deconflate_hr(hr, method = "published")$adjusted
#>    disease      raw adjusted      change estimand adjusted_for
#> 1       CK 1.500100 1.177580 -0.21499875    crude         <NA>
#> 2       CM 2.300000 1.903941 -0.17219973    crude         <NA>
#> 3       DA 2.851179 2.197930 -0.22911542    crude         <NA>
#> 4      DYS 1.258143 1.098377 -0.12698583    crude         <NA>
#> 5      LAM 1.744976 1.380683 -0.20876694    crude         <NA>
#> 6      MET 1.116444 1.012411 -0.09318284    crude         <NA>
#> 7       MF 2.999886 2.647637 -0.11742065    crude         <NA>
#> 8       OC 1.620000 1.458644 -0.09960240    crude         <NA>
#> 9      PTB 2.310508 2.047235 -0.11394597    crude         <NA>
#> 10      RP 1.599928 1.284496 -0.19715366    crude         <NA>
#> 11     SCK 1.920000 1.675253 -0.12747225    crude         <NA>
#> 12     SCM 1.449996 1.254928 -0.13453033    crude         <NA>
```
