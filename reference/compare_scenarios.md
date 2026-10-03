# Compare scenarios

Adjusts several versions of a model (e.g. with different association,
three-way or interaction assumptions) and tabulates the aggregate and
disease rankings side by side.

## Usage

``` r
compare_scenarios(..., method = "simultaneous", valuation = NULL)
```

## Arguments

- ...:

  Named
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  objects.

- method:

  Adjustment method (`"global"` is used automatically for models with
  interactions or three-way terms).

- valuation:

  Optional valuation list (see
  [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md)).

## Value

A list with `totals` (one row per scenario, with the reason in `failed`
when a scenario could not be adjusted or gave an undefined result) and
`by_disease` (contribution and rank per disease and scenario).

## Examples

``` r
base <- example_supplement()
strong <- set_association(base, "d1", "d3", 3)
compare_scenarios(base = base, strong_d1_d3 = strong)
#> $totals
#>       scenario    total failed rel_to_first
#> 1         base 2.109229   <NA>   0.00000000
#> 2 strong_d1_d3 1.959735   <NA>  -0.07087606
#> 
#> $by_disease
#>       scenario disease contribution rank
#> 1         base      d1   0.21432505    3
#> 2         base      d2   0.50806194    2
#> 3         base      d3   1.38684189    1
#> 4 strong_d1_d3      d1   0.06530914    3
#> 5 strong_d1_d3      d2   0.52786982    2
#> 6 strong_d1_d3      d3   1.36655608    1
#> 
```
