# Compare scenarios

Adjusts several versions of a model (e.g. with different association or
interaction assumptions) and tabulates total burden and disease rankings
side by side.

## Usage

``` r
compare_scenarios(
  ...,
  method = "simultaneous",
  outcome = NULL,
  economics = NULL
)
```

## Arguments

- ...:

  Named
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  objects.

- method:

  Adjustment method (`"global"` is used automatically for models with
  interactions).

- outcome, economics:

  Metric, as in
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).

## Value

A list with `totals` (one row per scenario) and `by_disease` (burden and
rank per disease and scenario).

## Examples

``` r
base <- example_supplement()
strong <- set_association(base, "d1", "d3", 3)
compare_scenarios(base = base, strong_d1_d3 = strong)
#> $totals
#>       scenario      total rel_to_first
#> 1         base 0.02109229   0.00000000
#> 2 strong_d1_d3 0.01959735  -0.07087606
#> 
#> $by_disease
#>       scenario disease       burden rank
#> 1         base      d1 0.0021432505    3
#> 2         base      d2 0.0050806194    2
#> 3         base      d3 0.0138684189    1
#> 4 strong_d1_d3      d1 0.0006530914    3
#> 5 strong_d1_d3      d2 0.0052786982    2
#> 6 strong_d1_d3      d3 0.0136655608    1
#> 
```
