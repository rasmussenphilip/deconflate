# Screen disease pairs for influential associations

Re-runs the adjustment with each disease pair's association changed, one
pair at a time, and reports how much the total burden and the ranking of
diseases change. Pairs without a specified association (independent by
default, or unknown) are set to each of `or_values`; specified odds
ratios and risk ratios are multiplied by each of `multipliers`. This
identifies associations worth estimating, as in Rasmussen et al. (2022),
Fig. 3 and Rasmussen et al. (2024), Fig. 7C-D.

## Usage

``` r
screen_associations(
  model,
  method = "simultaneous",
  or_values = c(0.5, 2),
  multipliers = c(0.5, 2),
  pairs = NULL,
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

- or_values:

  Odds ratios tried for unspecified pairs.

- multipliers:

  Factors applied to specified odds ratios or risk ratios.

- pairs:

  Optional character vector of pairs (`"d1:d2"`) to screen; default all
  pairs.

- outcome:

  Outcome whose burden is the metric (ignored if `economics` is given);
  default the first outcome.

- economics:

  Optional list (`observed`, `unit_value`, `additional`) to use total
  monetary losses as the metric (see
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)).

## Value

A data frame sorted by the absolute relative change in the total, with
the scenario, total, change, relative change, the largest shift in any
disease's rank and the Spearman correlation of disease burdens with the
baseline. Infeasible scenarios have `NA` totals.

## Examples

``` r
screen_associations(example_supplement(), or_values = c(0.5, 3))
#> <cm_screen> 6 scenarios; baseline total 0.02109
#>   pair    status scenario  total    change rel_change max_rank_shift rank_corr
#>  d2:d3 specified OR x 0.5 0.0231  0.001974     0.0936              0         1
#>  d2:d3 specified   OR x 2 0.0194 -0.001643    -0.0779              0         1
#>  d1:d3 specified   OR x 2 0.0201 -0.000963    -0.0457              0         1
#>  d1:d3 specified OR x 0.5 0.0220  0.000859     0.0407              0         1
#>  d1:d2 specified   OR x 2 0.0205 -0.000582    -0.0276              0         1
#>  d1:d2 specified OR x 0.5 0.0216  0.000549     0.0260              0         1
```
