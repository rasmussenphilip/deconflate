# Contribution table

One row per disease with raw and adjusted impacts and the disease's
contribution to the aggregate (its Shapley share, including any
interaction share). With `valuation`, the attributed gap and its value
are added (see
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
and
[`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)).

## Usage

``` r
contribution_table(result, valuation = NULL)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- valuation:

  Optional list with `observed`, `direction`, `effect` and (optionally)
  `unit_value`, passed to
  [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
  and
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md).

## Value

A data frame.

## Examples

``` r
contribution_table(deconflate(example_supplement()),
                   list(observed = 10000, direction = "decrease", effect = "percent",
                        unit_value = 0.3))
#>   disease raw adjusted      change contribution      main interaction     share
#> 1      d1 2.5 2.143250 -0.14269981    0.2143250 0.2143250           0 0.1016130
#> 2      d2 5.0 3.387080 -0.32258408    0.5080619 0.5080619           0 0.2408757
#> 3      d3 7.5 6.934209 -0.07543874    1.3868419 1.3868419           0 0.6575113
#>         gap     value
#> 1  21.89431  6.568292
#> 2  51.90090 15.570271
#> 3 141.67238 42.501715
```
