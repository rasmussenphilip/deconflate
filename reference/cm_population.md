# Describe the population: diseases and their associations

A population holds everything that is shared by all impact analyses: the
disease probabilities, the pairwise associations and, optionally,
three-way association scenarios. Combine it with an impact vector in
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
or with several in
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md).

## Usage

``` r
cm_population(
  diseases,
  associations = NULL,
  three_way = NULL,
  missing_associations = c("independent", "unknown"),
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

- three_way:

  Optional
  [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md)
  object (global model only).

- missing_associations:

  How to treat pairs without an association: `"independent"` (default;
  odds ratio 1, as in Rasmussen et al. 2022) or `"unknown"`
  (unconstrained). These are different assumptions.

- adjusted_associations:

  Covariate-adjusted association measures are not marginal 2x2
  associations. `"error"` (default) rejects them; `"use_as_marginal"`
  uses them as if they were marginal, which is an approximation, and
  records this.

## Value

A `cm_population` object.

## Examples

``` r
pop <- cm_population(
  cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
)
pop
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
```
