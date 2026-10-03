# Combine a population with an impact vector

Combine a population with an impact vector

## Usage

``` r
cm_model(
  population,
  impacts = NULL,
  interactions = NULL,
  associations = NULL,
  three_way = NULL,
  missing_associations = c("independent", "unknown"),
  adjusted_associations = c("error", "use_as_marginal")
)
```

## Arguments

- population:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).
  For convenience, a
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object can be given instead, together with `associations`,
  `missing_associations` and `three_way`.

- impacts:

  A
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  object with one value per disease.

- interactions:

  Optional
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  object.

- associations, three_way, missing_associations, adjusted_associations:

  Used only when `population` is a
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object; see
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).

## Value

A `cm_model` object (which is also a `cm_population`).

## Examples

``` r
m <- cm_model(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
  associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
)
m
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
#>   Impacts: (unlabelled) [%]
#>   Estimands: crude: 3
```
