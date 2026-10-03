# Adjust hazard ratios for comorbidity (hazard-ratio adapter)

Adjusts raw hazard ratios of culling (or mortality) for the hazard
ratios of associated diseases. This is a separate adapter from the
additive engine
([`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)),
because hazard ratios combine multiplicatively.

## Usage

``` r
deconflate_hr(
  model,
  method = c("snapshot", "first_order", "published"),
  joint = NULL,
  warn = TRUE,
  ...
)
```

## Arguments

- model:

  A
  [`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md).

- method:

  `"snapshot"`, `"first_order"` or `"published"`.

- joint:

  Optional
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result, checked against the population. The snapshot method uses it;
  the other methods keep it in the result for
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).

- warn:

  Logical: warn when adjusted hazard ratios cross 1 or are not positive
  and finite?

- ...:

  Passed to
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md).

## Value

A `cm_hr_result` with `adjusted` (raw and adjusted hazard ratios),
`diagnostics`, `joint`, `model` and `method`.

## Methods

- `"snapshot"` (default): a snapshot hazard-multiplier model. An
  animal's hazard is `h0 * exp(sum_i beta_i * D_i)`, and the raw hazard
  ratio of disease `i` is taken to be the ratio of the average hazard
  multiplier among animals with and without `i`, over the fitted
  distribution of disease combinations
  ([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md))
  at the start of follow-up. The `beta`s are solved so that these ratios
  equal the raw hazard ratios. For an adjusted estimate, the ratio is
  computed within strata of its adjustment set and combined across
  strata with Mantel-Haenszel-type weights; with `adjusted_for = "all"`
  the hazard ratio is used as it is.

- `"first_order"`: the log-linear approximation, `log(HR_raw) = A beta`,
  with `A` as in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  (pairwise tables only).

- `"published"`: Rasmussen et al. (2024): `HR - 1` adjusted with eq. 16
  (`"snapshot_crude"` hazard ratios only). For reproduction and
  comparison.

## What the snapshot model is not

A Cox hazard ratio estimated over follow-up is not, in general, the
snapshot ratio: animals with high hazards leave first, so the mixture of
disease combinations among survivors changes over time, and the marginal
hazard ratio changes with it. The snapshot model therefore does not give
an exact de-conflation of published Cox coefficients. It is exact for
its own estimand (instantaneous marginal ratios at baseline), and a
reasonable approximation when follow-up is short relative to the hazards
or the diseases are rare. Results are hazard ratios; turn them into
culling attributable to disease with
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).

## Examples

``` r
hr <- cm_hr_model(example_supplement(),
                  cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3),
                                   estimand = "snapshot_crude"))
deconflate_hr(hr)$adjusted
#>   disease raw adjusted      change       estimand adjusted_for
#> 1      d1 1.5 1.383153 -0.07789802 snapshot_crude         <NA>
#> 2      d2 2.0 1.890932 -0.05453423 snapshot_crude         <NA>
#> 3      d3 1.3 1.143987 -0.12000996 snapshot_crude         <NA>
deconflate_hr(hr, method = "first_order")$adjusted
#>   disease raw adjusted      change       estimand adjusted_for
#> 1      d1 1.5 1.402925 -0.06471659 snapshot_crude         <NA>
#> 2      d2 2.0 1.887440 -0.05628013 snapshot_crude         <NA>
#> 3      d3 1.3 1.169124 -0.10067407 snapshot_crude         <NA>
```
