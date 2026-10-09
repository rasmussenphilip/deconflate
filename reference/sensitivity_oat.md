# One-at-a-time sensitivity analysis

Varies each input by `variation` (e.g. +/- 20%) with everything else
fixed, as in Rasmussen et al. (2024), Fig. 7, and reports the resulting
range of the adjusted aggregate.

## Usage

``` r
sensitivity_oat(
  model,
  method = "simultaneous",
  variation = 0.2,
  inputs = c("prob", "assoc", "impact")
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- method:

  Adjustment method.

- variation:

  Relative change applied downwards and upwards.

- inputs:

  Which inputs to vary: `"prob"`, `"assoc"`, `"impact"`.

## Value

A `cm_oat` data frame with one row per input: the input key, totals at
the low and high values and the swing, sorted by swing.

## Examples

``` r
sensitivity_oat(example_supplement())
#>         input value total_low total_high      swing  rel_swing
#> 1   impact:d3  7.50  1.844070   2.374387 0.53031715 0.25142703
#> 2     prob:d3  0.20  1.851628   2.374637 0.52300893 0.24796215
#> 3   impact:d2  5.00  1.998423   2.220035 0.22161264 0.10506809
#> 4     prob:d2  0.15  2.026335   2.197276 0.17094138 0.08104449
#> 5 assoc:d2:d3  3.00  2.170224   2.062015 0.10820854 0.05130242
#> 6   impact:d1  2.50  2.063348   2.155110 0.09176175 0.04350488
#> 7     prob:d1  0.10  2.069218   2.149697 0.08047878 0.03815555
#> 8 assoc:d1:d3  1.00  2.138876   2.084124 0.05475254 0.02595856
#> 9 assoc:d1:d2  2.00  2.127728   2.093846 0.03388228 0.01606382
```
