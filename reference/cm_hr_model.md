# Combine a population with hazard ratios

Combine a population with hazard ratios

## Usage

``` r
cm_hr_model(population, hazard_ratios)
```

## Arguments

- population:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  (or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  whose impacts are ignored).

- hazard_ratios:

  A
  [`cm_hazard_ratios()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hazard_ratios.md)
  object with one value per disease.

## Value

A `cm_hr_model` object.
