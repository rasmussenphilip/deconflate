# Batch sampler: several analyses on shared population draws

Builds one
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
per analysis of a
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
object. In
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md),
each draw of the disease probabilities and associations is shared by all
analyses, and each analysis draws its own impacts and interactions. Draw
identifiers are the same across analyses.

## Usage

``` r
cm_batch_sampler(
  analyses,
  diseases = list(),
  associations = list(),
  impacts = list(),
  interactions = list()
)
```

## Arguments

- analyses:

  A
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
  object.

- diseases, associations:

  Named lists of `cm_dist` for the shared population inputs (as in
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)).

- impacts, interactions:

  Named lists, one element per analysis, each a named list of `cm_dist`
  as in
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

## Value

A `cm_batch_sampler` object.
