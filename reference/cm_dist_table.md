# Build distributions from an uncertainty table

Converts a table with columns `key`, `dist` and `p1`-`p4` (see the
Columns section of
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md))
into a named list of `cm_dist` objects. Keys are kept as given.

## Usage

``` r
cm_dist_table(x)
```

## Arguments

- x:

  A path to a CSV file or a data frame.

## Value

A named list of `cm_dist` objects (names = keys).

## Examples

``` r
cm_dist_table(data.frame(key = c("impact:d1", "assoc:d1:d2"),
                         dist = c("normal", "lognormal_ci"),
                         p1 = c(2.5, 2), p2 = c(0.5, 1.4), p3 = c(NA, 2.9)))
#> $`impact:d1`
#> <cm_dist> normal(mean = 2.5, sd = 0.5, lower = -Inf, upper = Inf), mean 2.5
#> 
#> $`assoc:d1:d2`
#> <cm_dist> lognormal(meanlog = 0.6931, sdlog = 0.1858), mean 2.035
#> 
```
