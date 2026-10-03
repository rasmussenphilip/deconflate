# Convert adjusted excess risks back to hazard ratios

Convert adjusted excess risks back to hazard ratios

## Usage

``` r
adjusted_hr(
  result,
  conversion,
  outcome = "culling",
  method = c("proportional_hazards", "published")
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
  used to build the impacts.

- outcome:

  Outcome label of the culling impacts.

- method:

  `"proportional_hazards"` (default): inverts the proportional-hazards
  conversion with each disease's unexposed risk
  ([`excess_to_hr()`](https://rasmussenphilip.github.io/deconflate/reference/excess_to_hr.md)).
  `"published"`: rescales the raw hazard ratio by the ratio of adjusted
  to raw excess risk (Rasmussen et al. 2022, eq. 23).

## Value

A data frame with raw and adjusted excess risks and hazard ratios.
