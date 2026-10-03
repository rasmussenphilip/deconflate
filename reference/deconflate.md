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
    ([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md))
    and solves the full equations, including any interactions in the
    model. Unknown pairs are left unconstrained. Without interactions
    and unknown pairs, it equals `"simultaneous"`.

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
have the corresponding conflation terms set to zero.

The diagnostics report, per outcome: the maximum absolute difference
between the supplied raw impacts and those reconstructed from the
adjusted impacts under the additive (or interaction) model; the number
of adjusted impacts whose sign differs from the raw impact; and the
condition number of the conflation matrix. For the published method a
non-zero reconstruction residual is expected: it measures the
approximation error.

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
```
