# Convert culling (or mortality) hazard ratios to period risks

Converts a hazard ratio for disease `i` into the period risk of the
event (e.g. culling within a year) among animals with and without the
disease, given the disease's prevalence and the overall period risk.

## Usage

``` r
hr_to_risk(
  hr,
  prevalence,
  overall_risk,
  method = c("proportional_hazards", "or_approx", "overall_odds")
)
```

## Arguments

- hr:

  Hazard ratio(s).

- prevalence:

  Disease probability (same length as `hr`, or recycled).

- overall_risk:

  Overall period risk of the event in the population.

- method:

  - `"proportional_hazards"` (default): assumes proportional hazards
    over the period, so `risk_exposed = 1 - (1 - risk_unexposed)^hr`,
    with `risk_unexposed` solved so that the population risk equals
    `overall_risk`.

  - `"or_approx"`: treats the hazard ratio as an odds ratio in the 2x2
    table of disease by event (as in Rasmussen et al. 2022, section
    2.4.4 and Table 6). Provided to reproduce the published values.

  - `"overall_odds"`: treats the hazard ratio as an odds ratio relative
    to the overall risk, `risk_exposed = hr * r / (hr * r + 1 - r)`,
    with the overall risk `r` as the reference (`risk_unexposed = r`),
    as in the loss calculations of Rasmussen et al. (2024). `prevalence`
    is not used.

## Value

A data frame with `risk_exposed`, `risk_unexposed` (the reference risk)
and `excess` (their difference).

## Examples

``` r
# Displaced abomasum, Rasmussen et al. (2022) Table 6: HR 3.83, prevalence 0.03,
# culling rate 0.27; published excess probability 0.31.
hr_to_risk(3.83, 0.03, 0.27, method = "or_approx")
#>     hr prevalence overall_risk risk_exposed risk_unexposed    excess
#> 1 3.83       0.03         0.27    0.5744262      0.2605848 0.3138415
hr_to_risk(3.83, 0.03, 0.27)
#>     hr prevalence overall_risk risk_exposed risk_unexposed    excess
#> 1 3.83       0.03         0.27    0.6799849      0.2573201 0.4226648
# Rasmussen et al. (2024): adjusted HR relative to the overall culling risk
hr_to_risk(2.75, NA, 0.27, method = "overall_odds")
#>     hr prevalence overall_risk risk_exposed risk_unexposed    excess
#> 1 2.75         NA         0.27    0.5042445           0.27 0.2342445
```
