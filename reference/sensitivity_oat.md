# One-at-a-time sensitivity analysis

Varies each input by `variation` (e.g. +/- 20%) with everything else
fixed, as in Rasmussen et al. (2024), Fig. 7, and reports the resulting
range of the total burden.

## Usage

``` r
sensitivity_oat(
  model,
  method = "simultaneous",
  variation = 0.2,
  inputs = c("prob", "assoc", "impact"),
  outcome = NULL,
  economics = NULL
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

- outcome, economics:

  Metric, as in
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).

## Value

A `cm_screen`-like data frame with one row per input: the input key,
totals at the low and high values and the swing, sorted by swing.

## Examples

``` r
sensitivity_oat(example_supplement())
#>             input value  total_low total_high        swing  rel_swing
#> 1 impact:yield:d3  7.50 0.01844070 0.02374387 0.0053031715 0.25142703
#> 2         prob:d3  0.20 0.01851628 0.02374637 0.0052300893 0.24796215
#> 3 impact:yield:d2  5.00 0.01998423 0.02220035 0.0022161264 0.10506809
#> 4         prob:d2  0.15 0.02026335 0.02197276 0.0017094138 0.08104449
#> 5     assoc:d2:d3  3.00 0.02170224 0.02062015 0.0010820854 0.05130242
#> 6 impact:yield:d1  2.50 0.02063348 0.02155110 0.0009176175 0.04350488
#> 7         prob:d1  0.10 0.02069218 0.02149697 0.0008047878 0.03815555
#> 8     assoc:d1:d3  1.00 0.02138876 0.02084124 0.0005475254 0.02595856
#> 9     assoc:d1:d2  2.00 0.02127728 0.02093846 0.0003388228 0.01606382
```
