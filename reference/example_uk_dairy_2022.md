# UK dairy example from Rasmussen et al. (2022)

Thirteen endemic diseases and conditions of UK dairy cattle, with
cow-level prevalence (Table 2), inter-disease odds ratios (Table 3, all
other pairs independent), and three analyses: milk yield (Table 4, %
decrease), calving interval (Table 5, % increase) and culling (Table 6).
Use
[`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md)
for the observed means and unit values (Table 1), and
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
for the published tables.

## Usage

``` r
example_uk_dairy_2022(yield_sck = 3.05)
```

## Arguments

- yield_sck:

  Raw yield impact of subclinical ketosis, in percent (the printed 3.05,
  or 100 \* 340 / 8737 = 3.89 from note i of Table 4).

## Value

A
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
object with analyses `yield`, `fertility` and `culling`, and attribute
`"hazard_ratios"` (a
[`cm_hazard_ratios()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hazard_ratios.md)).

## Details

The culling analysis reproduces the paper's approach, which converted
hazard ratios to excess annual culling risks by treating them as odds
ratios. That conversion is not supported for new analyses (use
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md));
the hazard ratios are attached as attribute `"hazard_ratios"` for that
purpose.

Disease ids: CO cystic ovary, DA displaced abomasum, DYS dystocia, FAS
fasciolosis, GIN gastrointestinal nematodes, LAM lameness, MAS mastitis,
MET metritis, MF milk fever, NEO neosporosis, PTB paratuberculosis, RP
retained placenta, SCK subclinical ketosis.

## Examples

``` r
uk <- example_uk_dairy_2022()
deconflate(uk, method = "published")
#> <cm_results> 3 analyses: yield, fertility, culling
#> 
#> == yield ==
#> <cm_result> method: published; milk yield loss [% decrease]
#> 
#>  disease  raw adjusted  change
#>       CO 0.00  0.00000      NA
#>       DA 4.04  2.36694 -0.4141
#>      DYS 4.05  2.93975 -0.2741
#>      FAS 7.33  7.33000  0.0000
#>      GIN 3.28  3.28000  0.0000
#>      LAM 5.54  4.83063 -0.1280
#>      MAS 4.57  3.76265 -0.1767
#>      MET 3.95  2.40914 -0.3901
#>       MF 0.41  0.09989 -0.7564
#>      NEO 4.20  4.20000  0.0000
#>      PTB 5.90  4.43066 -0.2490
#>       RP 7.38  6.10015 -0.1734
#>      SCK 3.05  1.91867 -0.3709
#> 
#> Raw sum: 7.168; adjusted total: 6.046
#> Diagnostics: residual 6.14e-01, condition number 2.08, sign changes 0
#> 
#> == fertility ==
#> <cm_result> method: published; calving interval increase [% increase]
#> 
#>  disease   raw adjusted   change
#>       CO 11.26  11.1282 -0.01171
#>       DA  0.00   0.0000       NA
#>      DYS  6.96   6.0198 -0.13509
#>      FAS  0.00   0.0000       NA
#>      GIN  1.20   1.2000  0.00000
#>      LAM 12.47  11.8896 -0.04655
#>      MAS  0.00   0.0000       NA
#>      MET  4.74   4.0989 -0.13526
#>       MF  0.00   0.0000       NA
#>      NEO  7.21   7.2100  0.00000
#>      PTB  5.79   3.8550 -0.33419
#>       RP  2.74   1.6765 -0.38815
#>      SCK  1.50   0.5495 -0.63364
#> 
#> Raw sum: 7.573; adjusted total: 6.907
#> Diagnostics: residual 1.52e+00, condition number 2.08, sign changes 0
#> 
#> == culling ==
#> <cm_result> method: published; excess annual culling risk (hazard ratio treated as an odds ratio; legacy) [proportional increase in the culling rate]
#> 
#>  disease     raw adjusted  change
#>       CO 0.00000  0.00000      NA
#>       DA 0.31384  0.21891 -0.3025
#>      DYS 0.14205  0.11896 -0.1626
#>      FAS 0.00000  0.00000      NA
#>      GIN 0.00000  0.00000      NA
#>      LAM 0.25565  0.22552 -0.1179
#>      MAS 0.21307  0.16920 -0.2059
#>      MET 0.17386  0.12261 -0.2948
#>       MF 0.20567  0.15552 -0.2438
#>      NEO 0.09889  0.09889  0.0000
#>      PTB 0.19637  0.13437 -0.3157
#>       RP 0.00000  0.00000      NA
#>      SCK 0.15726  0.09665 -0.3854
#> 
#> Raw sum: 0.2499; adjusted total: 0.1976
#> Diagnostics: residual 5.08e-02, condition number 2.08, sign changes 0
#> 
```
