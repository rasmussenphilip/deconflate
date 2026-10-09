# Screen disease pairs for influential associations

Re-runs the adjustment with each disease pair's association changed, one
pair at a time, and reports how much the total and the ranking of
diseases change. Pairs without an association (unknown) are set to each
of `or_values`; given associations of any measure are multiplied by each
of `multipliers`. Scenarios whose inputs are invalid or infeasible are
kept, with `NA` results and the reason, so that no scenario is skipped
silently. This identifies associations worth estimating, as in Rasmussen
et al. (2022), Fig. 3 and Rasmussen et al. (2024), Fig. 7C-D.

## Usage

``` r
screen_associations(
  model,
  method = "auto",
  event_model = FALSE,
  overall_risk = NULL,
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

  As in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

- event_model, overall_risk:

  As in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  (a distribution of the overall risk is used at its mean).

- or_values:

  Odds ratios tried for pairs without an association.

- multipliers:

  Factors applied to given associations.

- pairs:

  Optional character vector of pairs (`"d1:d2"`); default all.

## Value

A `cm_screen` data frame sorted by the absolute relative change in the
total, with the scenario, total, change, relative change, the largest
shift in any disease's rank, the Spearman correlation of disease
contributions with the baseline, and `failed` (the reason a scenario
could not be run).

## Details

The model needs no associations: with none, the baseline has every
disease independent (adjusted impacts equal the raw ones), and the
screen shows which pairs would matter if they were associated. Each
scenario is adjusted as
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
would (point estimates; the method is chosen for the scenario's own
inputs, see `method`).

The total is the adjusted aggregate for additive impacts, and the risk
attributable to disease for event impacts (`event_model = TRUE`);
disease rankings use the contributions to it.

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

# Without any association estimates: which pairs would matter?
m <- example_supplement()
m$associations <- NULL
screen_associations(m, or_values = c(0.5, 2, 4))
#> <cm_screen> 9 scenarios; baseline total 2.5
#>   pair  status scenario total  change rel_change max_rank_shift rank_corr
#>  d2:d3 unknown   OR = 4  2.09 -0.4143    -0.1657              0         1
#>  d1:d3 unknown   OR = 4  2.28 -0.2223    -0.0889              0         1
#>  d2:d3 unknown   OR = 2  2.28 -0.2158    -0.0863              0         1
#>  d2:d3 unknown OR = 0.5  2.70  0.1984     0.0794              0         1
#>  d1:d2 unknown   OR = 4  2.35 -0.1488    -0.0595              0         1
#>  d1:d3 unknown   OR = 2  2.38 -0.1171    -0.0468              0         1
#>  d1:d3 unknown OR = 0.5  2.60  0.1042     0.0417              0         1
#>  d1:d2 unknown   OR = 2  2.43 -0.0733    -0.0293              0         1
#>  d1:d2 unknown OR = 0.5  2.56  0.0588     0.0235              0         1
#>  failed
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
```
