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
  unit_value = NULL,
  joint = NULL,
  allocate = TRUE,
  max_present = Inf
)
```

## Arguments

- result:

  A
  [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
  result.

- overall_risk:

  Overall period risk of the event, as a proportion (e.g. `0.25` for an
  annual culling rate of 25%).

- unit_value:

  Optional value per animal removed (e.g. replacement cost less salvage
  value); adds `value` columns.

- joint:

  Optional
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result; by default the result's own joint distribution or a new fit.
  It is checked against the population.

- allocate:

  Logical: allocate the attributable risk to diseases?

- max_present:

  Passed to
  [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
  (default: no skipping).

## Value

A `cm_attributable` list with `summary` (overall, disease-free and
attributable risk, attributable fraction, value, and any unallocated
part), `by_disease` and `baseline_hazard`.

## Details

Within the period, an animal with disease combination `x` has a constant
hazard `h0 * exp(sum_i beta[i] * x[i])`, with `beta = log(adjusted HR)`.
The baseline hazard `h0` is chosen so that the population risk, averaged
over the distribution of disease combinations
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)),
equals `overall_risk`. The disease-free risk is `1 - exp(-h0)`, and the
attributable risk is `overall_risk - (1 - exp(-h0))`. An animal's risk
cannot exceed 1, so the attributable risk is smaller than the sum of
per-disease excess risks when diseases co-occur.

The attributable risk is allocated to diseases by Shapley values over
all disease combinations
([`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)).
Any combinations skipped with `max_present` are reported as unallocated.

This uses the same snapshot model as
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
(constant hazards within the period, multiplicative hazard ratios, no
change in the mixture of diseases over the period).

## Examples

``` r
hr <- cm_hr_model(example_supplement(), cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3)))
attributable_risk(deconflate_hr(hr), overall_risk = 0.25, unit_value = 1300)
#> <cm_attributable> snapshot hazard-multiplier model
#>   Overall risk 0.25; disease-free risk 0.2134; attributable 0.03665 (14.7% of the overall risk)
#>   Value: 47.64
#> 
#>  disease hr_adjusted attributable  share  value
#>       d1       1.383     0.007371 0.2011  9.582
#>       d2       1.891     0.023488 0.6409 30.534
#>       d3       1.144     0.005791 0.1580  7.528
```
