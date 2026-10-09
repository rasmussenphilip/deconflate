# Plots

Base-graphics plots for the main result types:

- `plot(<cm_result>)`: raw and adjusted impacts by disease, with 95%
  intervals of the adjusted impacts when the result has draws.

- `plot(<cm_event_result>)`: adjusted hazard ratios by disease (with
  intervals when the result has draws), and the raw hazard and rate
  ratios.

- `plot_burden()`: each disease's share of the aggregate (including its
  share of interaction effects), or of the risk attributable to disease.

- `plot(<cm_screen>)`: the most influential scenarios from
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
  [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
  or
  [`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md).

- `plot(<cm_oat>)`: tornado plot from
  [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md).

## Usage

``` r
# S3 method for class 'cm_result'
plot(x, ...)

# S3 method for class 'cm_event_result'
plot(x, ...)

plot_burden(result, ...)

# S3 method for class 'cm_screen'
plot(x, top = 15, ...)

# S3 method for class 'cm_oat'
plot(x, top = 15, ...)

# S3 method for class 'cm_threshold'
plot(x, items = NULL, ...)
```

## Arguments

- x, result:

  The object to plot.

- ...:

  Passed to the underlying graphics function.

- top:

  Number of rows to show.

- items:

  For `cm_threshold` plots: the items to show (default all).

## Value

The input, invisibly.

## Examples

``` r
res <- deconflate(example_supplement())
plot(res)

plot_burden(res)
```
