# UK dairy example from Rasmussen et al. (2022)

Thirteen endemic diseases and conditions of UK dairy cattle, with
cow-level prevalence (Table 2), inter-disease odds ratios (Table 3, all
other pairs independent), and two analyses: milk yield (Table 4, %
decrease) and calving interval (Table 5, % increase). Use
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
for the published tables, including the economic values.

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
object with analyses `yield` and `fertility`, and attribute
`"hazard_ratios"` (a
[`cm_hazard_ratios()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hazard_ratios.md)).

## Details

The culling hazard ratios of Table 6 are attached as attribute
`"hazard_ratios"` (with estimand `"snapshot_crude"`, an assumption) for
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md).
The paper converted them to excess annual culling risks by treating them
as odds ratios; that historical calculation is only available inside
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).

Disease ids: CO cystic ovary, DA displaced abomasum, DYS dystocia, FAS
fasciolosis, GIN gastrointestinal nematodes, LAM lameness, MAS mastitis,
MET metritis, MF milk fever, NEO neosporosis, PTB paratuberculosis, RP
retained placenta, SCK subclinical ketosis.

## Examples

``` r
uk <- example_uk_dairy_2022()
deconflate(uk)
#> Warning: Adjusted impacts change sign for MF. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
#> Warning: Adjusted impacts change sign for SCK. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
#> <cm_results> 2 analyses: yield, fertility
#> 
#> == yield ==
#> <cm_result> method: simultaneous; milk yield loss [% decrease]
#> 
#>  disease  raw adjusted   change
#>       CO 0.00  -0.6059       NA
#>       DA 4.04   2.2006 -0.45530
#>      DYS 4.05   2.9495 -0.27174
#>      FAS 7.33   7.3300  0.00000
#>      GIN 3.28   3.2800  0.00000
#>      LAM 5.54   5.0329 -0.09154
#>      MAS 4.57   4.0846 -0.10622
#>      MET 3.95   1.8610 -0.52887
#>       MF 0.41  -0.5145 -2.25479
#>      NEO 4.20   4.2000  0.00000
#>      PTB 5.90   4.1323 -0.29962
#>       RP 7.38   6.5817 -0.10817
#>      SCK 3.05   1.6425 -0.46148
#> 
#> Raw sum: 7.168; adjusted total: 5.982
#> Diagnostics: residual 4.44e-16, condition number 2.08, sign changes 1 (MF)
#> 
#> == fertility ==
#> <cm_result> method: simultaneous; calving interval increase [% increase]
#> 
#>  disease   raw adjusted    change
#>       CO 11.26  11.4517  0.017022
#>       DA  0.00  -0.3935        NA
#>      DYS  6.96   6.0841 -0.125855
#>      FAS  0.00   0.0000        NA
#>      GIN  1.20   1.2000  0.000000
#>      LAM 12.47  12.3725 -0.007815
#>      MAS  0.00  -0.9409        NA
#>      MET  4.74   4.6170 -0.025941
#>       MF  0.00   0.2951        NA
#>      NEO  7.21   7.2100  0.000000
#>      PTB  5.79   3.0435 -0.474354
#>       RP  2.74   1.1739 -0.571558
#>      SCK  1.50  -0.9642 -1.642802
#> 
#> Raw sum: 7.573; adjusted total: 6.448
#> Diagnostics: residual 8.88e-16, condition number 2.08, sign changes 1 (SCK)
#> 
```
