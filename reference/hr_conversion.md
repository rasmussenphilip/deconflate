# Convert culling hazard ratios into excess-risk impacts

Converts disease-specific culling (or mortality) hazard ratios into the
excess period risk among diseased animals, using
[`hr_to_risk()`](https://rasmussenphilip.github.io/deconflate/reference/hr_to_risk.md)
with each disease's probability from the model. Diseases without a
hazard ratio get `hr = 1` (no excess). Use
[`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md)
to turn the result into impacts for
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
and
[`adjusted_hr()`](https://rasmussenphilip.github.io/deconflate/reference/adjusted_hr.md)
to convert adjusted impacts back to hazard ratios.

## Usage

``` r
hr_conversion(
  diseases,
  hr,
  overall_risk,
  method = c("proportional_hazards", "or_approx")
)
```

## Arguments

- diseases:

  A
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  object.

- hr:

  Named numeric vector of hazard ratios (names = disease ids).

- overall_risk:

  Overall period risk of the event in the population (e.g. the annual
  culling rate as a proportion).

- method:

  `"proportional_hazards"` (default) or `"or_approx"` (the published
  approach of Rasmussen et al. 2022); see
  [`hr_to_risk()`](https://rasmussenphilip.github.io/deconflate/reference/hr_to_risk.md).

## Value

A `cm_hr` data frame with one row per disease.

## Examples

``` r
conv <- hr_conversion(example_supplement(), c(d1 = 2, d2 = 1.5), overall_risk = 0.27)
conv
#>   disease  hr prevalence overall_risk risk_exposed risk_unexposed    excess
#> 1      d1 2.0       0.10         0.27    0.4392845      0.2511906 0.1880939
#> 2      d2 1.5       0.15         0.27    0.3565998      0.2547177 0.1018821
#> 3      d3 1.0       0.20         0.27    0.2700000      0.2700000 0.0000000
#>                 method
#> 1 proportional_hazards
#> 2 proportional_hazards
#> 3 proportional_hazards
as_impacts(conv)
#>   disease outcome     value    scale       units direction adjusted_for source
#> 1      d1 culling 0.1880939 absolute probability  increase         <NA>   <NA>
#> 2      d2 culling 0.1018821 absolute probability  increase         <NA>   <NA>
#> 3      d3 culling 0.0000000 absolute probability  increase         <NA>   <NA>
#>   input_scale
#> 1    absolute
#> 2    absolute
#> 3    absolute
```
