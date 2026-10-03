# Convert an adjusted excess risk back to a hazard ratio

Inverse of the proportional-hazards conversion in
[`hr_to_risk()`](https://rasmussenphilip.github.io/deconflate/reference/hr_to_risk.md),
for reporting adjusted culling impacts as hazard ratios.

## Usage

``` r
excess_to_hr(excess, risk_unexposed)
```

## Arguments

- excess:

  Excess period risk among diseased animals.

- risk_unexposed:

  Period risk among animals without the disease.

## Value

Hazard ratio(s):
`log(1 - (risk_unexposed + excess)) / log(1 - risk_unexposed)`.
