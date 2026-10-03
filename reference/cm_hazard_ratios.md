# Describe raw culling (or mortality) hazard ratios

Hazard ratios are not additive impacts, so they are adjusted by a
separate model,
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
(the snapshot hazard-multiplier model), outside the additive engine. Its
estimands are named explicitly, because a published hazard ratio is not
automatically one of them.

## Usage

``` r
cm_hazard_ratios(
  disease,
  value,
  estimand,
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

  Required: `"snapshot_crude"` or `"snapshot_stratified"` (see
  Estimands), one per row or recycled. There is no default, so that the
  assumption is stated.

- adjusted_for:

  For `"snapshot_stratified"`: the diseases the estimate was stratified
  (adjusted) for, separated by `";"`, or `"all"` (every other disease).

- source:

  Optional citation.

## Value

A `cm_hazard_ratios` data frame.

## Estimands

- `"snapshot_crude"`: the ratio of the average hazard among animals with
  the disease to that among animals without it, at the start of
  follow-up, in the population described by the probabilities and
  associations (unadjusted for other diseases).

- `"snapshot_stratified"`: the same ratio within strata of the diseases
  in `adjusted_for`, combined across strata with Mantel-Haenszel-type
  weights; with `adjusted_for = "all"` it is the disease's own hazard
  multiplier.

A Cox hazard ratio estimated over follow-up is a different quantity (see
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md),
"What the snapshot model is not"). Entering one as a snapshot estimand
is an approximation, which is better when follow-up is short relative to
the hazards; it is the user's assumption, not a conversion the package
makes. A Cox model adjusted for other diseases estimates a conditional
coefficient, which matches `"snapshot_stratified"` with the same
adjustment set only under the snapshot model's assumptions.

## Examples

``` r
cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), estimand = "snapshot_crude")
#>   disease value       estimand adjusted_for source
#> 1      d1   1.5 snapshot_crude         <NA>   <NA>
#> 2      d2   2.0 snapshot_crude         <NA>   <NA>
#> 3      d3   1.3 snapshot_crude         <NA>   <NA>
```
