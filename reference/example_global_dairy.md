# Global dairy inputs from Rasmussen et al. (2024), at their central values

The 12 diseases of Rasmussen et al. (2024): global (herd-size weighted)
lactational incidence (prevalence for PTB), pooled inter-disease odds
ratios (Table 3, pairs not listed are independent), and raw impacts on
yield (% decrease), fertility (% increase in calving interval) and,
optionally, culling (hazard ratio minus 1). Central values are means of
normal distributions and modes of PERT distributions.

## Usage

``` r
example_global_dairy(inputs = c("analysis", "tables"), culling = TRUE)
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`; see the section above.

- culling:

  Include the culling outcome?

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
with outcomes `"yield"`, `"fertility"` and, if `culling = TRUE`,
`"culling"`.

## Two versions of the inputs

- `inputs = "analysis"` (default) uses the inputs of the published
  analysis code (1st revision), which reproduce Table 5:

  - disease probabilities are fixed at `1 - exp(-incidence)`, using the
    unrounded global mean incidence. Subclinical mastitis (SCM) was
    entered without conversion (0.4094 rather than 0.3360); this is kept
    here (`type = "probability"`);

  - impacts are the unrounded central values (e.g. clinical ketosis
    yield 0.4322% rather than 0.43%);

  - culling is entered as HR - 1, as in the analysis. Metritis uses the
    analysis value (PERT mode 1.116) rather than Table 4 (normal, mean
    1.05).

- `inputs = "tables"` uses Tables 2-4 as printed, with incidence
  converted for every disease except PTB.

## Culling

Culling impacts are hazard ratios minus 1 on the `"absolute"` scale, so
that
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
adjusts the excess hazard ratio. Convert results back with
`adjusted_hr(res, method = "excess_hr")`. These impacts are not excess
risks, so do not pass them to
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md);
the paper converts adjusted hazard ratios to excess culling risk with
`hr_to_risk(..., method = "overall_odds")`.

The paper's Table 5 reports means of adjusted impacts over Monte Carlo
draws
([`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)),
which differ from adjusting the central values (e.g. displaced abomasum
yield: 1.18 vs 0.79).

## Examples

``` r
res <- deconflate(example_global_dairy(), method = "published")
adjusted_hr(res, method = "excess_hr")
#>    disease       hr   excess excess_adjusted hr_adjusted
#> 1       CK 1.500100 0.500100      0.17758038    1.177580
#> 2       CM 2.300000 1.300000      0.90394063    1.903941
#> 3       DA 2.851179 1.851179      1.19792993    2.197930
#> 4      DYS 1.258143 0.258143      0.09837667    1.098377
#> 5      LAM 1.744976 0.744976      0.38068270    1.380683
#> 6      MET 1.116444 0.116444      0.01241058    1.012411
#> 7       MF 2.999886 1.999886      1.64763745    2.647637
#> 8       OC 1.620000 0.620000      0.45864412    1.458644
#> 9      PTB 2.310508 1.310508      1.04723492    2.047235
#> 10      RP 1.599928 0.599928      0.28449635    1.284496
#> 11     SCK 1.920000 0.920000      0.67525327    1.675253
#> 12     SCM 1.449996 0.449996      0.25492756    1.254928
```
