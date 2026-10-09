# Describe the population: diseases and their associations

A population holds the disease probabilities, the pairwise associations
and, optionally, three-way association scenarios. Combine it with one
impact table in
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Usage

``` r
cm_population(
  diseases,
  associations = NULL,
  three_way = NULL,
  adjusted_associations = c("error", "use_as_marginal")
)
```

## Arguments

- diseases:

  A
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object.

- associations:

  Optional
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
  object.
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  needs at least one association; the sensitivity tools (e.g.
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md))
  also run without any.

- three_way:

  Optional
  [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md)
  object.

- adjusted_associations:

  Covariate-adjusted association measures are not marginal 2x2
  associations. `"error"` (default) rejects them; `"use_as_marginal"`
  uses them as if they were marginal, which is an approximation, and
  records this.

## Value

A `cm_population` object.

## Details

Pairs without an association are unknown: the global model fills in
their association from the others (see
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)).

## Examples

``` r
pop <- cm_population(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
)
pop
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 (2 with an association, 1 unknown)
```
