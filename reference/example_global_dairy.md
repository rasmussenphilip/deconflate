# Global dairy inputs from Rasmussen et al. (2024)

The 12 diseases of Rasmussen et al. (2024): global (herd-size weighted)
lactational incidence (prevalence for PTB), pooled inter-disease odds
ratios (Table 3; pairs not listed are independent), and raw impacts on
yield (% decrease) and fertility (% increase in calving interval), at
their central values (means of normal distributions, modes of PERT
distributions). Each impact type is one analysis.

## Usage

``` r
example_global_dairy(inputs = c("analysis", "tables"))
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`.

## Value

A
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
object.

## Two versions of the inputs

- `inputs = "analysis"` (default) uses the inputs of the published
  analysis code (1st revision), which reproduce Table 5: disease
  probabilities fixed at `1 - exp(-incidence)` at the unrounded global
  means, except subclinical mastitis (entered without conversion, 0.4094
  rather than 0.3360, as in the analysis), and unrounded impacts.

- `inputs = "tables"` uses Tables 2-4 as printed, with incidence
  converted for every disease except PTB.

## Culling

The culling hazard ratios are available as a hazard-ratio model from
[`example_global_dairy_hr()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy_hr.md)
(see
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)).
The 2024 analysis adjusted hazard ratios minus 1 as if they were
additive impacts; that historical calculation is only available inside
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).

## Examples

``` r
gd <- example_global_dairy()
deconflate(gd, method = "published")$yield$adjusted
#>    disease       raw   adjusted     change estimand adjusted_for
#> 1       CK 0.4321944 0.02473297 -0.9427735    crude         <NA>
#> 2       CM 3.2499000 1.33194993 -0.5901566    crude         <NA>
#> 3       DA 2.8369300 0.79378291 -0.7201965    crude         <NA>
#> 4      DYS 4.9190880 3.56881898 -0.2744958    crude         <NA>
#> 5      LAM 4.8061000 2.53240737 -0.4730848    crude         <NA>
#> 6      MET 5.6130850 2.84095146 -0.4938699    crude         <NA>
#> 7       MF 0.5365131 0.06897532 -0.8714378    crude         <NA>
#> 8       OC 3.7478390 2.63859343 -0.2959694    crude         <NA>
#> 9      PTB 4.3000000 3.22917066 -0.2490301    crude         <NA>
#> 10      RP 4.1986640 2.25771112 -0.4622787    crude         <NA>
#> 11     SCK 8.3964720 7.10307185 -0.1540409    crude         <NA>
#> 12     SCM 6.2931840 5.58961905 -0.1117979    crude         <NA>
```
