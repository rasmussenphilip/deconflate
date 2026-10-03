# Describe raw impact estimates (one impact vector)

One analysis adjusts one set of compatible, additive impact estimates:
one value per disease, all in the same units (e.g. kg of milk per cow,
percent of yield, days, euros, welfare scores). The engine does not
convert units; results come back in the units supplied. For several
types of impact, use one impact vector each
([`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)).

## Usage

``` r
cm_impacts(
  disease,
  value,
  estimand = "crude",
  adjusted_for = NA_character_,
  source = NA_character_,
  label = NULL,
  units = NULL
)
```

## Arguments

- disease:

  Character vector of disease ids (one row per disease).

- value:

  Numeric raw impacts. Every disease in the model needs a value (use 0
  for no impact).

- estimand:

  `"crude"` (default) or `"adjusted_linear"`, one per row or recycled.

- adjusted_for:

  For `"adjusted_linear"`: the diseases the estimate was adjusted for,
  separated by `";"` (e.g. `"LAM; CM"`), or `"all"`.

- source:

  Optional citation.

- label, units:

  Optional analysis-level metadata (e.g. `label = "milk yield loss"`,
  `units = "% of yield"`), carried into the results.

## Value

A `cm_impacts` data frame (attributes `label` and `units`).

## Estimands

Each value must be one of the supported estimands:

- `"crude"`: the difference in the outcome between animals with and
  without the disease (unadjusted for other diseases).

- `"adjusted_linear"`: the coefficient of the disease in an additive
  (linear) regression of the outcome on the disease and the diseases in
  `adjusted_for`, in the same source population. The adjustment then
  uses the population projection of the omitted diseases (see
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
  `adjusted_for = "all"` means every other disease in the model.

`adjusted_for` is only used with `estimand = "adjusted_linear"`; it is
never used to infer the estimand. Other adjusted estimands (e.g. matched
or propensity-score estimates) are not supported. The probabilities and
associations must describe the population the estimates come from.

## Examples

``` r
cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%")
#>   disease value estimand adjusted_for source
#> 1      d1   2.5    crude         <NA>   <NA>
#> 2      d2   5.0    crude         <NA>   <NA>
#> 3      d3   7.5    crude         <NA>   <NA>
```
