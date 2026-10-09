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
deconflate(gd)$yield$adjusted
#> Warning: Adjusted impacts change sign for CK, CM, DA, MF. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
#> Warning: Adjusted impacts change sign for CK, DA, LAM, SCK, SCM. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
#>    disease       raw   adjusted       change estimand adjusted_for
#> 1       CK 0.4321944 -5.3709141 -13.42707923    crude         <NA>
#> 2       CM 3.2499000 -0.3380725  -1.10402550    crude         <NA>
#> 3       DA 2.8369300 -2.6316316  -1.92763360    crude         <NA>
#> 4      DYS 4.9190880  4.3096749  -0.12388742    crude         <NA>
#> 5      LAM 4.8061000  1.9792779  -0.58817379    crude         <NA>
#> 6      MET 5.6130850  2.7909168  -0.50278379    crude         <NA>
#> 7       MF 0.5365131 -1.5399737  -3.87033758    crude         <NA>
#> 8       OC 3.7478390  3.2174755  -0.14151181    crude         <NA>
#> 9      PTB 4.3000000  3.9410642  -0.08347344    crude         <NA>
#> 10      RP 4.1986640  2.4733053  -0.41093041    crude         <NA>
#> 11     SCK 8.3964720  8.2785938  -0.01403901    crude         <NA>
#> 12     SCM 6.2931840  6.5761493   0.04496377    crude         <NA>
compare_methods(gd$models$yield)$totals
#>         method  raw_sum adjusted_total
#> 1    published 9.931174       7.280553
#> 2 simultaneous 9.931174       7.494192
#> 3       global 9.931174       7.494192
```
