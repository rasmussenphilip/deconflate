# Describe raw culling (or mortality) hazard ratios

Hazard ratios are not additive impacts, so they are adjusted by a
separate adapter,
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md),
outside the additive engine.

## Usage

``` r
cm_hazard_ratios(
  disease,
  value,
  estimand = "crude",
  adjusted_for = NA_character_,
  source = NA_character_
)
```

## Arguments

- disease:

  Character vector of disease ids (one row per disease; use 1 for a
  disease with no effect).

- value:

  Positive hazard ratios.

- estimand:

  `"crude"` (unadjusted for other diseases) or `"adjusted"` (from a
  model that included the diseases in `adjusted_for`), one per row or
  recycled.

- adjusted_for:

  For `"adjusted"`: the diseases the estimate was adjusted for,
  separated by `";"`, or `"all"` (every other disease).

- source:

  Optional citation.

## Value

A `cm_hazard_ratios` data frame.

## Examples

``` r
cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3))
#>   disease value estimand adjusted_for source
#> 1      d1   1.5    crude         <NA>   <NA>
#> 2      d2   2.0    crude         <NA>   <NA>
#> 3      d3   1.3    crude         <NA>   <NA>
```
