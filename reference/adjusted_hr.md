# Convert adjusted culling impacts back to hazard ratios

Convert adjusted culling impacts back to hazard ratios

## Usage

``` r
adjusted_hr(
  result,
  conversion = NULL,
  outcome = "culling",
  method = c("proportional_hazards", "published", "excess_hr")
)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- conversion:

  The
  [`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md)
  used to build the impacts. Not needed for `method = "excess_hr"`.

- outcome:

  Outcome label of the culling impacts.

- method:

  - `"proportional_hazards"` (default): inverts the proportional-hazards
    conversion with each disease's unexposed risk
    ([`excess_to_hr()`](https://rasmussenphilip.github.io/deconflate/reference/excess_to_hr.md)).

  - `"published"`: rescales the raw hazard ratio by the ratio of
    adjusted to raw excess risk (Rasmussen et al. 2022, eq. 23).

  - `"excess_hr"`: the impacts are hazard ratios minus 1, adjusted
    directly (Rasmussen et al. 2024; see
    [`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)),
    so the adjusted hazard ratio is the adjusted impact plus 1.

## Value

A data frame with raw and adjusted impacts (`excess`, `excess_adjusted`:
excess risks, or HR - 1 for `"excess_hr"`) and hazard ratios (`hr`,
`hr_adjusted`).
