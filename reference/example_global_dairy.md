# Global dairy inputs from Rasmussen et al. (2024)

The 12 diseases of Rasmussen et al. (2024): global (herd-size weighted)
lactational incidence (prevalence for PTB), pooled inter-disease odds
ratios (Table 3), and the raw impacts of one outcome: milk yield (%
decrease), fertility (% increase in calving interval) or culling (hazard
ratios, event impacts). The model includes the input distributions of
the analysis, so
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
reports intervals (use `n_draws = 0` for point estimates only).

## Usage

``` r
example_global_dairy(
  outcome = c("yield", "fertility", "culling"),
  inputs = c("analysis", "tables")
)
```

## Arguments

- outcome:

  `"yield"`, `"fertility"` or `"culling"`.

- inputs:

  `"analysis"` or `"tables"`.

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Details

Pairs not in Table 3 are unknown, as for any input in this package, so
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
uses the global method. The 2024 analysis treated them as independent;
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
does the same.

## Two versions of the inputs

- `inputs = "analysis"` (default) uses the inputs of the published
  analysis code (1st revision), which reproduce Table 5: disease
  probabilities fixed at `1 - exp(-incidence)` at the unrounded global
  means, except subclinical mastitis (entered without conversion, 0.4094
  rather than 0.3360, as in the analysis), and unrounded impacts.

- `inputs = "tables"` uses Tables 2-4 as printed, with incidence
  converted for every disease except PTB; disease probabilities are
  uncertain too.

Central values are the means of normal distributions and the modes of
PERT distributions; normal distributions of odds ratios are truncated at
zero.

## Culling

`outcome = "culling"` gives the hazard ratios as event impacts (estimand
`"snapshot_crude"`, an assumption), for
`deconflate(..., event_model = TRUE, overall_risk = 0.2366)` (the global
culling risk of the 2024 analysis). The 2024 analysis adjusted hazard
ratios minus 1 as if they were additive impacts; that historical
calculation is only available inside
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).

## Examples

``` r
gd <- example_global_dairy()
deconflate(gd, n_draws = 0)
#> Warning: Adjusted impacts change sign for CK, CM, DA, MF. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
#> <cm_result> method: global; milk yield loss [% decrease]
#> 
#>  disease    raw adjusted    change
#>       CK 0.4322  -5.0307 -12.64000
#>       CM 3.2499  -0.1396  -1.04294
#>       DA 2.8369  -2.6153  -1.92186
#>      DYS 4.9191   3.9893  -0.18901
#>      LAM 4.8061   2.0544  -0.57255
#>      MET 5.6131   2.3826  -0.57553
#>       MF 0.5365  -2.2102  -5.11949
#>       OC 3.7478   2.6804  -0.28483
#>      PTB 4.3000   3.1214  -0.27409
#>       RP 4.1987   2.1445  -0.48925
#>      SCK 8.3965   7.8863  -0.06076
#>      SCM 6.2932   5.9495  -0.05461
#> 
#> Raw sum: 9.931; adjusted total: 6.915
#> Diagnostics: residual 1.78e-15, condition number 3.24, sign changes 4 (CK, CM, DA, MF)
#> Unknown pairs: 28 without an association; the global fit gave them odds ratios from 1.08 to 2.37 (see $unknown_pairs).
#> 
#> Notes:
#> * The global method was used because of 28 pairs without an association
#>   (unknown).
# \donttest{
deconflate(example_global_dairy("culling"), event_model = TRUE,
           overall_risk = 0.2366, n_draws = 0)
#> Warning: Adjusted hazard ratios are on the other side of 1 from the raw estimates for CK, DYS, MET. The raw estimates are weaker than the associated diseases alone would produce; check the estimands and source populations.
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure   raw adjusted_hr
#>       CK      HR 1.500      0.9771
#>       CM      HR 2.300      1.8182
#>       DA      HR 2.851      1.7870
#>      DYS      HR 1.258      0.9540
#>      LAM      HR 1.745      1.2401
#>      MET      HR 1.116      0.6928
#>       MF      HR 3.000      2.4247
#>       OC      HR 1.620      1.4022
#>      PTB      HR 2.311      1.9208
#>       RP      HR 1.600      1.1832
#>      SCK      HR 1.920      1.6452
#>      SCM      HR 1.450      1.1547
#> 
#> Overall risk 0.2366; disease-free risk 0.128; attributable to disease 0.1086 (45.9% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable     share
#>       CK   -0.0001265 -0.001164
#>       CM    0.0304920  0.280670
#>       DA    0.0027374  0.025197
#>      DYS   -0.0004590 -0.004225
#>      LAM    0.0087303  0.080360
#>      MET   -0.0053919 -0.049630
#>       MF    0.0047215  0.043460
#>       OC    0.0065598  0.060381
#>      PTB    0.0138246  0.127252
#>       RP    0.0035096  0.032305
#>      SCK    0.0341130  0.313999
#>      SCM    0.0099293  0.091396
#> 
#> Diagnostics: residual 7.33e-15, condition number 3.82, sign changes 3 (CK, DYS, MET)
#> Unknown pairs: 28 without an association; the global fit gave them odds ratios from 1.08 to 2.37 (see $unknown_pairs).
# }
```
