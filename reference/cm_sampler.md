# Build a Monte Carlo sampler from input distributions

Creates a function that, for each draw, copies `model` and replaces
selected inputs with random draws from
[distributions](https://rasmussenphilip.github.io/deconflate/reference/distributions.md).
Inputs without a distribution keep their values. The result is passed to
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md).

## Usage

``` r
cm_sampler(
  model,
  diseases = list(),
  associations = list(),
  impacts = list(),
  interactions = list(),
  three_way = list()
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- diseases:

  Named list of `cm_dist`, names = disease ids.

- associations:

  Named list of `cm_dist`, names `"d1:d2"` (either order) matching rows
  of the model's associations.

- impacts:

  Named list of `cm_dist`, names = disease ids.

- interactions:

  Named list of `cm_dist`, names `"d1:d2"` (either order) matching rows
  of the model's interactions.

- three_way:

  Named list of `cm_dist`, names `"d1:d2:d3"` (any order) matching rows
  of the model's three-way terms. Three-way terms affect the global
  method only.

## Value

A function of the draw index returning a
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
with attribute `specs` (the distributions, keyed as in `params`). It
also accepts `u` (named uniforms, for stratified sampling) and `values`
(named input values that replace draws, for importance sampling and for
batch runs);
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
uses these.

## Details

Draws are on the scale the inputs were entered on:

- diseases: the `value` given to
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  (e.g. an incidence rate, converted to a probability as specified by
  its `type`);

- associations: the association measure (e.g. odds ratio);

- three-way terms: the ratio of conditional odds ratios (see
  [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md));

- impacts and interactions: their own units.

A draw that produces an impossible input (a probability outside (0, 1),
a non-positive odds ratio, a conditional probability below 0 or above 1)
is rejected by
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
and counted.

Keys (as in `params` of
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)):
`prob:<disease>`, `assoc:<d1>:<d2>`, `three:<d1>:<d2>:<d3>`,
`impact:<disease>` and `inter:<d1>:<d2>`.

## Examples

``` r
s <- cm_sampler(example_supplement(),
                associations = list("d1:d2" = dist_lognormal_ci(2, 1.4, 2.9)),
                impacts = list(d1 = dist_normal(2.5, 0.5)))
s(1)$associations
#>   disease1 disease2 measure    value adjusted adjusted_for source corrected n11
#> 1       d1       d2      OR 2.153526    FALSE         <NA>   <NA>     FALSE  NA
#> 2       d1       d3      OR 1.000000    FALSE         <NA>   <NA>     FALSE  NA
#> 3       d2       d3      OR 3.000000    FALSE         <NA>   <NA>     FALSE  NA
#>   n10 n01 n00
#> 1  NA  NA  NA
#> 2  NA  NA  NA
#> 3  NA  NA  NA
```
