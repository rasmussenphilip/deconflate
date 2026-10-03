# Global dairy inputs from Rasmussen et al. (2024), at their means

Global (herd-size weighted) lactational incidence (prevalence for PTB)
from Table 2, pooled inter-disease odds ratios from Table 3, and raw
yield (% decrease) and fertility (% increase in calving interval)
impacts from Table 4. Pairs absent from Table 3 are treated as
independent, as in the paper. Culling impacts (hazard ratios) are
omitted; see
[`hr_to_risk()`](https://rasmussenphilip.github.io/deconflate/reference/hr_to_risk.md).

## Usage

``` r
example_global_dairy()
```

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
with outcomes `"yield"` and `"fertility"`.

## Details

The paper's Table 5 reports means of adjusted impacts over Monte Carlo
draws, which differ from adjusting the input means (e.g. displaced
abomasum yield: 1.18 vs 0.79).
