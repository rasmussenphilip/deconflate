# Adjust raw impact estimates for comorbidity

Raw impact estimates (comparisons of animals with and without a disease)
are treated as conflations of the disease's own impact and the impacts
of associated diseases. Under additive impacts, the raw impact of
disease `i` is

## Usage

``` r
deconflate(
  model,
  method = c("simultaneous", "published", "global"),
  joint = NULL,
  warn = TRUE,
  ...
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  with impacts.

- method:

  - `"simultaneous"` (default): solves the additive equations exactly
    (`(I + t(E)) m = m_raw`) using the pairwise 2x2 tables. No joint
    distribution is needed.

  - `"published"`: the proportional approximation of Rasmussen et al.
    (2022), eq. 16:
    `m[i] = m_raw[i]^2 / (m_raw[i] + sum_k E[k, i] * m_raw[k])`.
    Provided for reproduction and comparison; it can mask incompatible
    inputs because it cannot return negative impacts.

  - `"global"`: fits the maximum-entropy joint distribution
    ([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md);
    the "iterative" model) and solves the full equations, including any
    interactions in the model. Unknown pairs are left unconstrained.
    Without interactions and unknown pairs, it equals `"simultaneous"`
    for additive outcomes.

- joint:

  Optional pre-computed
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result for `method = "global"`.

- warn:

  Logical: warn when adjusted impacts change sign or the method produces
  non-finite values?

- ...:

  Passed to
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md).

## Value

A `cm_result` with elements `adjusted` (long data frame of raw and
adjusted impacts), `diagnostics`, `conflation` (the matrix `A` per
outcome, with `m_raw = A m + offset`), `joint_pairs` (matrix of
`P(j and k)`), the `model`, the `joint` fit (global only) and `method`.

## Details

`m_raw[i] = m[i] + sum_k E[k, i] * m[k]`

and, with pairwise interactions `delta`, additionally

`+ sum_k delta[i, k] * P(k | i) + sum_{j < k; j, k != i} delta[j, k] * (P(j, k | i) - P(j, k | not i))`.

Impacts flagged with `adjusted_for` in
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
have the corresponding conflation terms set to zero. This is not yet
supported for hazard-ratio outcomes with `method = "global"`.

The diagnostics report, per outcome: the maximum absolute difference
between the supplied raw impacts and those reconstructed from the
adjusted impacts under the model (on the log scale for hazard ratios,
and on the `HR - 1` scale for the published method); the number of
adjusted impacts whose sign differs from the raw impact (for hazard
ratios: on which side of 1 they lie); and the condition number of the
conflation matrix (or of the Jacobian, for the exact hazard model). For
the published method a non-zero reconstruction residual is expected: it
measures the approximation error.

## Hazard ratios

Impacts on the `"hazard_ratio"` scale (see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md))
are adjusted under a multiplicative (Cox-type) model, in which an
animal's hazard is `h0 * exp(sum_i beta[i] * D[i])` and the adjusted
hazard ratio of disease `i` is `exp(beta[i])`:

- `"global"` solves the model exactly: the raw hazard ratio of disease
  `i` is the ratio of the average hazard among animals with and without
  `i`, `E[exp(beta . D) | i] / E[exp(beta . D) | not i]`, over the
  fitted distribution of disease combinations. The equations are solved
  by Newton's method.

- `"simultaneous"` uses the first-order (log-linear) version,
  `log HR_raw = (I + t(E)) beta`, which needs no joint distribution.

- `"published"` reproduces Rasmussen et al. (2024): `HR - 1` is adjusted
  with eq. 16 and 1 is added back.

The adjusted values in the result are hazard ratios for all three
methods. Use
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
to turn them into culling (or mortality) attributable to disease.

## Examples

``` r
res <- deconflate(example_supplement())
res$adjusted
#>   outcome disease   raw   adjusted      change      scale             units
#> 1   yield      d1 0.025 0.02143250 -0.14269981 proportion units/animal/year
#> 2   yield      d2 0.050 0.03387080 -0.32258408 proportion units/animal/year
#> 3   yield      d3 0.075 0.06934209 -0.07543874 proportion units/animal/year
#>   direction
#> 1  decrease
#> 2  decrease
#> 3  decrease

# Culling hazard ratios: exact multiplicative model
m <- example_supplement()
m$impacts <- combine_impacts(m$impacts,
  cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), outcome = "culling",
             scale = "hazard_ratio"))
deconflate(m, method = "global")$adjusted
#>   outcome disease   raw   adjusted      change        scale             units
#> 1   yield      d1 0.025 0.02143250 -0.14269981   proportion units/animal/year
#> 2   yield      d2 0.050 0.03387080 -0.32258408   proportion units/animal/year
#> 3   yield      d3 0.075 0.06934209 -0.07543874   proportion units/animal/year
#> 4 culling      d1 1.500 1.38315297 -0.07789802 hazard_ratio      hazard ratio
#> 5 culling      d2 2.000 1.89093155 -0.05453423 hazard_ratio      hazard ratio
#> 6 culling      d3 1.300 1.14398705 -0.12000996 hazard_ratio      hazard ratio
#>   direction
#> 1  decrease
#> 2  decrease
#> 3  decrease
#> 4  increase
#> 5  increase
#> 6  increase
```
