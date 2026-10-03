# Describe raw (unadjusted) disease impact estimates

Describe raw (unadjusted) disease impact estimates

## Usage

``` r
cm_impacts(
  disease,
  value,
  outcome = "impact",
  scale = "proportion",
  units = NA_character_,
  direction = "decrease",
  adjusted_for = NA_character_,
  source = NA_character_
)
```

## Arguments

- disease:

  Character vector of disease ids.

- value:

  Numeric raw impact estimates.

- outcome:

  Character outcome label(s), e.g. `"yield"`, `"fertility"`. Each
  disease needs exactly one impact per outcome (use `0` for no impact).

- scale:

  The scale of `value`, constant within an outcome:

  - `"proportion"`: proportional change relative to the disease-free
    value (e.g. `0.025`);

  - `"percent"`: converted to a proportion;

  - `"absolute"`: in `units` (e.g. excess culling risk);

  - `"hazard_ratio"`: a hazard ratio (e.g. for culling or mortality),
    adjusted on the log scale; see
    [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
    and
    [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).

- units:

  Optional units label (required for `"absolute"`).

- direction:

  `"decrease"` if disease lowers the outcome (e.g. yield) or
  `"increase"` if it raises it (e.g. calving interval). Used by
  [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md).
  Set to `"increase"` for hazard ratios.

- adjusted_for:

  Diseases the raw estimate was already adjusted for, separated by `";"`
  (e.g. `"LAM; CM"`). Their conflation terms are removed from the
  adjustment for this impact.

- source:

  Optional citation.

## Value

A `cm_impacts` data frame.

## Examples

``` r
cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), outcome = "yield",
           scale = "percent")
#>   disease outcome value      scale units direction adjusted_for source
#> 1      d1   yield 0.025 proportion  <NA>  decrease         <NA>   <NA>
#> 2      d2   yield 0.050 proportion  <NA>  decrease         <NA>   <NA>
#> 3      d3   yield 0.075 proportion  <NA>  decrease         <NA>   <NA>
#>   input_scale
#> 1     percent
#> 2     percent
#> 3     percent
cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.3, 1.1), outcome = "culling",
           scale = "hazard_ratio")
#>   disease outcome value        scale        units direction adjusted_for source
#> 1      d1 culling   1.5 hazard_ratio hazard ratio  increase         <NA>   <NA>
#> 2      d2 culling   2.3 hazard_ratio hazard ratio  increase         <NA>   <NA>
#> 3      d3 culling   1.1 hazard_ratio hazard ratio  increase         <NA>   <NA>
#>    input_scale
#> 1 hazard_ratio
#> 2 hazard_ratio
#> 3 hazard_ratio
```
