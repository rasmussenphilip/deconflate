# Culling (or mortality) attributable to disease

Converts adjusted hazard ratios into the part of an event's overall risk
(e.g. annual culling) that is attributable to disease, without counting
an animal with several diseases more than once, and allocates it to
diseases.

## Usage

``` r
attributable_risk(
  result,
  overall_risk,
  outcome = "culling",
  unit_value = NULL,
  joint = NULL,
  allocate = TRUE,
  max_present = 10L
)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result with a hazard-ratio outcome (see
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)).

- overall_risk:

  Overall period risk of the event, as a proportion (e.g. `0.25` for an
  annual culling rate of 25%).

- outcome:

  Outcome label.

- unit_value:

  Optional value per animal removed (e.g. replacement cost less salvage
  value); adds `value` columns.

- joint:

  Optional
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result; by default the result's own joint distribution (global method)
  or a new fit.

- allocate:

  Logical: allocate the attributable risk to diseases? This is the slow
  step for many co-occurring diseases.

- max_present:

  Passed to
  [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md).

## Value

A `cm_attributable` list with `summary` (overall, disease-free and
attributable risk, attributable fraction and value), `by_disease`
(adjusted hazard ratio, attributable risk, share and value),
`skipped_mass` and `baseline_hazard`.

## Details

Within the period, an animal with disease combination `x` has a constant
hazard `h0 * exp(sum_i beta[i] * x[i])`, with `beta = log(adjusted HR)`.
The baseline hazard `h0` is chosen so that the population risk, averaged
over the distribution of disease combinations
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)),
equals `overall_risk`:

`sum_x P(x) * (1 - exp(-h0 * exp(beta . x))) = overall_risk`.

The disease-free risk is `1 - exp(-h0)`, and the attributable risk is
`overall_risk - (1 - exp(-h0))`. An animal's risk cannot exceed 1, so
the attributable risk is smaller than the sum of per-disease excess
risks when diseases co-occur.

The attributable risk is allocated to diseases by Shapley values over
disease combinations
([`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)),
with the loss `1 - exp(-h0 * exp(beta . x)) - (1 - exp(-h0))`.

The model is consistent with `method = "global"` in
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md):
there, the adjusted hazard ratios reproduce the raw ones exactly over
the same joint distribution. Hazard ratios from the other methods are
used as they are.

## Examples

``` r
m <- example_supplement()
m$impacts <- combine_impacts(m$impacts,
  cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), outcome = "culling",
             scale = "hazard_ratio"))
res <- deconflate(m, method = "global")
attributable_risk(res, overall_risk = 0.25, unit_value = 1300)
#> <cm_attributable> outcome: culling
#>   Overall risk 0.25; disease-free risk 0.2134; attributable 0.03665 (14.7% of the overall risk)
#>   Value: 47.64
#> 
#>  disease hr_adjusted attributable  share  value
#>       d1       1.383     0.007371 0.2011  9.582
#>       d2       1.891     0.023488 0.6409 30.534
#>       d3       1.144     0.005791 0.1580  7.528
```
