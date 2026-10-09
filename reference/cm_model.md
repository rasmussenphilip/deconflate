# Combine a population with an impact table

Combine a population with an impact table

## Usage

``` r
cm_model(
  population,
  impacts = NULL,
  interactions = NULL,
  associations = NULL,
  three_way = NULL,
  adjusted_associations = c("error", "use_as_marginal"),
  distributions = NULL
)
```

## Arguments

- population:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).
  For convenience, a
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object can be given instead, together with `associations` and
  `three_way`.

- impacts:

  A
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  object with one value per disease.

- interactions:

  Optional
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  object (additive impacts only).

- associations, three_way, adjusted_associations:

  Used only when `population` is a
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object; see
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).

- distributions:

  Optional named list of distributions
  ([distributions](https://rasmussenphilip.github.io/deconflate/reference/distributions.md))
  for uncertain inputs, used when
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  runs draws. Names are keys: `"prob:<disease>"` (the disease's `value`,
  on the scale it was entered), `"assoc:<d1>:<d2>"` (the association, on
  its own measure), `"three:<d1>:<d2>:<d3>"`, `"impact:<disease>"` and
  `"inter:<d1>:<d2>"`.
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  fills this in from the `dist` columns of the tables;
  [`cm_dist_table()`](https://rasmussenphilip.github.io/deconflate/reference/cm_dist_table.md)
  builds it from a table of keys. When `population` is a model and
  `distributions` is not given, its distributions are kept, except those
  of impacts or interactions that are replaced.

## Value

A `cm_model` object (which is also a `cm_population`).

## Examples

``` r
m <- cm_model(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
  associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)),
  distributions = list("impact:d1" = dist_normal(2.5, 0.5))
)
m
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (2 with an association, 1 unknown)
#>   Impacts: (unlabelled) [%] (additive)
#>   Estimands: crude: 3
#>   Uncertain inputs (with a distribution): 1
```
