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
  outcome_correlation = NULL
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

  Named list of `cm_dist`, names `"outcome:disease"`.

- interactions:

  Named list of `cm_dist`, names `"outcome:d1:d2"`.

- outcome_correlation:

  Optional correlation matrix with row and column names equal to outcome
  labels. Within each disease, the impacts on these outcomes are drawn
  with this correlation (Gaussian copula); marginal distributions are
  unchanged. Impacts of different diseases stay independent.

## Value

A function of the draw index returning a
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
with attributes `specs` (the distributions, keyed as in `params` of
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md))
and `correlated` (keys drawn through the copula). The function also
accepts `u` (named uniforms, for stratified sampling) and `values`
(named input values that replace draws, for importance sampling);
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
uses these.

## Details

Draws are on the scale the inputs were entered on:

- diseases: the `value` given to
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  (e.g. an incidence rate, converted to a probability as specified by
  its `type`);

- associations: the association measure (e.g. odds ratio);

- impacts and interactions: the input scale (e.g. percent, if entered as
  percent).

A draw that produces an impossible input (a probability outside (0, 1),
a non-positive odds ratio) is rejected by
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
and counted.

## Examples

``` r
m <- example_supplement()
s <- cm_sampler(m, associations = list("d1:d2" = dist_lognormal_ci(2, 1.4, 2.9)))
s(1)$associations
#>   disease1 disease2 measure    value adjusted adjusted_for source n11 n10 n01
#> 1       d1       d2      OR 2.153526    FALSE         <NA>   <NA>  NA  NA  NA
#> 2       d1       d3      OR 1.000000    FALSE         <NA>   <NA>  NA  NA  NA
#> 3       d2       d3      OR 3.000000    FALSE         <NA>   <NA>  NA  NA  NA
#>   n00
#> 1  NA
#> 2  NA
#> 3  NA
```
