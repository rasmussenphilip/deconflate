# Plots

Base-graphics plots for the main result types:

- `plot(<cm_result>)`: raw and adjusted impacts by disease;
  `plot(<cm_results>)` one panel per analysis.

- `plot_burden()`: each disease's share of the aggregate (including its
  share of interaction effects); for several analyses, one bar per
  analysis.

- `plot(<cm_mc>)`: Monte Carlo means and intervals of adjusted impacts.

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

# S3 method for class 'cm_results'
plot(x, ...)

plot_burden(result, ...)

# S3 method for class 'cm_mc'
plot(x, probs = c(0.025, 0.975), method = NULL, ...)

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

- probs:

  Interval bounds for Monte Carlo plots.

- method:

  For Monte Carlo runs with several methods: the method to show (default
  the first).

- top:

  Number of rows to show.

- items:

  For `cm_threshold` plots: the items to show (default all).

## Value

The input, invisibly.

## Examples

``` r
res <- deconflate(example_uk_dairy_2022())
#> Warning: Adjusted impacts change sign for MF. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
#> Warning: Adjusted impacts change sign for SCK. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.
plot(res$yield)

plot_burden(res)
```
