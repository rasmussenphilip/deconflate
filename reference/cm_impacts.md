# Describe raw impact estimates

One impact table holds one raw estimate per disease for one outcome.
There are two kinds:

## Usage

``` r
cm_impacts(
  disease,
  value,
  estimand = NULL,
  adjusted_for = NA_character_,
  source = NA_character_,
  label = NULL,
  units = NULL,
  measure = NULL
)
```

## Arguments

- disease:

  Character vector of disease ids (one row per disease).

- value:

  Numeric raw impacts. Every disease in the model needs a value (for no
  effect: 0 for additive impacts and risk differences, 1 for ratios).

- estimand:

  Additive impacts: `"crude"` (default) or `"adjusted_linear"`. Event
  impacts: `"snapshot_crude"` or `"snapshot_stratified"` (required). One
  per row or recycled.

- adjusted_for:

  For `"adjusted_linear"` and `"snapshot_stratified"`: the diseases the
  estimate was adjusted for, separated by `";"` (e.g. `"LAM; CM"`), or
  `"all"`.

- source:

  Optional citation.

- label, units:

  Optional table-level metadata (e.g. `label = "milk yield loss"`,
  `units = "% of yield"`), carried into the results.

- measure:

  `NULL` (default) for additive impacts; for event impacts, one measure
  per row or recycled (see Event impacts).

## Value

A `cm_impacts` data frame (attributes `label`, `units` and `kind`,
`"additive"` or `"event"`).

## Details

- **Additive impacts** (no `measure`): amounts in the outcome's own
  units, the same for every row (e.g. kg of milk, percent of yield,
  days, a welfare score). The package does not convert units; results
  come back in the units supplied. Adjust them with
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

- **Event impacts** (`measure` given): comparisons of the risk of an
  event (e.g. death or culling) between animals with and without the
  disease, as hazard ratios, rate ratios, risk ratios, odds ratios or
  risk differences. Adjust them with
  `deconflate(..., event_model = TRUE)`, which uses the snapshot hazard
  model.

For several outcomes, make one impact table each and adjust each in
turn.

## Estimands of additive impacts

- `"crude"` (default): the difference in the outcome between animals
  with and without the disease (unadjusted for other diseases).

- `"adjusted_linear"`: the coefficient of the disease in an additive
  (linear) regression of the outcome on the disease and the diseases in
  `adjusted_for`, in the same source population. The adjustment then
  uses the population projection of the omitted diseases (see
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
  `adjusted_for = "all"` means every other disease in the model.

Other adjusted estimands (e.g. matched or propensity-score estimates)
are not supported. The probabilities and associations must describe the
population the estimates come from.

## Event impacts

`measure` is one of `"HR"` (hazard ratio), `"rate_ratio"`, `"RR"` (risk
ratio over the period), `"OR"` (odds ratio of the event over the period)
or `"RD"` (risk difference over the period); measures can be mixed. The
estimand must be stated (there is no default, because it is an
assumption):

- `"snapshot_crude"`: the comparison of animals with and without the
  disease at the start of the period, in the population described by the
  probabilities and associations (unadjusted for other diseases);

- `"snapshot_stratified"`: the same comparison within strata of the
  diseases in `adjusted_for`, combined across strata with
  Mantel-Haenszel-type weights; with `adjusted_for = "all"` a hazard
  ratio is the disease's own hazard multiplier.

A Cox hazard ratio estimated over follow-up is not exactly a snapshot
estimand (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
"Event impacts"); entering one as such is an approximation and the
user's assumption. Risk ratios, odds ratios and risk differences refer
to the period of the overall risk given to
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Examples

``` r
cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%")
#>   disease value estimand adjusted_for measure source
#> 1      d1   2.5    crude         <NA>    <NA>   <NA>
#> 2      d2   5.0    crude         <NA>    <NA>   <NA>
#> 3      d3   7.5    crude         <NA>    <NA>   <NA>
# Event impacts: culling, with mixed measures
cm_impacts(c("d1", "d2", "d3"), c(1.5, 1.3, 0.04), measure = c("HR", "RR", "RD"),
           estimand = "snapshot_crude", label = "culling")
#>   disease value       estimand adjusted_for measure source
#> 1      d1  1.50 snapshot_crude         <NA>      HR   <NA>
#> 2      d2  1.30 snapshot_crude         <NA>      RR   <NA>
#> 3      d3  0.04 snapshot_crude         <NA>      RD   <NA>
```
