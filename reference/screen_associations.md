# Screen disease pairs for influential associations

Re-runs the adjustment with each disease pair's association changed, one
pair at a time, and reports how much the aggregate and the ranking of
diseases change. Pairs without a specified association (independent by
default, or unknown) are set to each of `or_values`; specified
associations of any measure are multiplied by each of `multipliers`
(contingency tables via their odds ratio). Scenarios whose inputs are
invalid or infeasible are kept, with `NA` results and the reason, so
that no measure is skipped silently. This identifies associations worth
estimating, as in Rasmussen et al. (2022), Fig. 3 and Rasmussen et al.
(2024), Fig. 7C-D.

## Usage

``` r
screen_associations(
  model,
  method = "simultaneous",
  or_values = c(0.5, 2),
  multipliers = c(0.5, 2),
  pairs = NULL
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- method:

  Adjustment method.

- or_values:

  Odds ratios tried for unspecified pairs.

- multipliers:

  Factors applied to specified associations.

- pairs:

  Optional character vector of pairs (`"d1:d2"`); default all.

## Value

A `cm_screen` data frame sorted by the absolute relative change in the
total, with the scenario, total, change, relative change, the largest
shift in any disease's rank, the Spearman correlation of disease
contributions with the baseline, and `failed` (the reason a scenario
could not be run).

## Examples

``` r
screen_associations(example_supplement(), or_values = c(0.5, 3))
#> <cm_screen> 6 scenarios; baseline total 2.109
#>   pair    status scenario total  change rel_change max_rank_shift rank_corr
#>  d2:d3 specified OR x 0.5  2.31  0.1974     0.0936              0         1
#>  d2:d3 specified   OR x 2  1.94 -0.1643    -0.0779              0         1
#>  d1:d3 specified   OR x 2  2.01 -0.0963    -0.0457              0         1
#>  d1:d3 specified OR x 0.5  2.20  0.0859     0.0407              0         1
#>  d1:d2 specified   OR x 2  2.05 -0.0582    -0.0276              0         1
#>  d1:d2 specified OR x 0.5  2.16  0.0549     0.0260              0         1
#>  failed
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
```
