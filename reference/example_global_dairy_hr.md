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
deconflate_hr(hr)$adjusted
#> Warning: Adjusted hazard ratios cross 1 for CK, MET. The raw hazard ratios are smaller than the associated diseases alone would produce; check the estimands and source populations.
#>    disease      raw  adjusted      change       estimand adjusted_for
#> 1       CK 1.500100 0.9821306 -0.34528988 snapshot_crude         <NA>
#> 2       CM 2.300000 1.8193536 -0.20897670 snapshot_crude         <NA>
#> 3       DA 2.851179 1.8476664 -0.35196408 snapshot_crude         <NA>
#> 4      DYS 1.258143 1.0505432 -0.16500491 snapshot_crude         <NA>
#> 5      LAM 1.744976 1.2230210 -0.29911873 snapshot_crude         <NA>
#> 6      MET 1.116444 0.7198427 -0.35523620 snapshot_crude         <NA>
#> 7       MF 2.999886 2.6601486 -0.11325012 snapshot_crude         <NA>
#> 8       OC 1.620000 1.5579660 -0.03829256 snapshot_crude         <NA>
#> 9      PTB 2.310508 2.0165455 -0.12722851 snapshot_crude         <NA>
#> 10      RP 1.599928 1.1816995 -0.26140459 snapshot_crude         <NA>
#> 11     SCK 1.920000 1.7107120 -0.10900416 snapshot_crude         <NA>
#> 12     SCM 1.449996 1.2109005 -0.16489389 snapshot_crude         <NA>
```
