# Several impact analyses on one population

Bundles several impact vectors (e.g. milk yield, calving interval and
welfare) that share the same diseases and associations. Each is adjusted
separately; Monte Carlo runs on a batch sampler
([`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md))
reuse the same disease and association draws for every analysis.

## Usage

``` r
cm_analyses(population, ..., interactions = list())
```

## Arguments

- population:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  (or a
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  whose impacts are dropped).

- ...:

  Named
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  objects, one per analysis.

- interactions:

  Optional named list of
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  objects, with names matching the analyses.

## Value

A `cm_analyses` object with `population` and `models` (a named list of
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
objects).

## Examples

``` r
pop <- example_supplement()
a <- cm_analyses(pop,
  yield = cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
  fertility = cm_impacts(c("d1", "d2", "d3"), c(1, 2, 0), units = "%"))
deconflate(a)
#> <cm_results> 2 analyses: yield, fertility
#> 
#> == yield ==
#> <cm_result> method: simultaneous; yield [%]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    2.143 -0.14270
#>       d2 5.0    3.387 -0.32258
#>       d3 7.5    6.934 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> 
#> == fertility ==
#> <cm_result> method: simultaneous; fertility [%]
#> 
#>  disease raw adjusted    change
#>       d1   1   0.7881 -0.211902
#>       d2   2   2.0119  0.005927
#>       d3   0  -0.3361        NA
#> 
#> Raw sum: 0.4; adjusted total: 0.3134
#> Diagnostics: residual 0.00e+00, condition number 1.53, sign changes 0
#> 
```
