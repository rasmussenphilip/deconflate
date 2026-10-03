# Turn a hazard-ratio conversion into impacts

Turn a hazard-ratio conversion into impacts

## Usage

``` r
as_impacts(x, outcome = "culling", scale = c("absolute", "proportion"), ...)
```

## Arguments

- x:

  A
  [`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md)
  result.

- outcome:

  Outcome label.

- scale:

  How the excess risk enters the productivity gap:

  - `"absolute"` (default): excess risk in probability units, so the
    disease-free risk is `observed - sum(excess * P)`. Use with an
    observed risk given as a proportion (e.g. `0.27`).

  - `"proportion"`: excess risk treated as a proportional increase in
    the observed rate, so the disease-free rate is `observed / (1 + L)`,
    as in Rasmussen et al. (2022). Use this to reproduce the published
    results.

- ...:

  Passed to
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  (e.g. `source`).

## Value

A
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
object with direction `"increase"`.
