# Adjust raw impact estimates for comorbidity

Raw impact estimates (comparisons of animals with and without a disease)
are treated as conflations of the disease's own impact and the impacts
of associated diseases. For one vector of additive impacts `b` (in any
units), the raw estimates satisfy

## Usage

``` r
deconflate(model, ...)

# S3 method for class 'cm_model'
deconflate(
  model,
  method = c("simultaneous", "published", "global"),
  joint = NULL,
  warn = TRUE,
  feasibility = c("screen", "lp", "none"),
  ...
)

# S3 method for class 'cm_analyses'
deconflate(
  model,
  method = c("simultaneous", "published", "global"),
  joint = NULL,
  ...
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  or a
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
  object (each analysis is adjusted in turn).

- ...:

  Passed to
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  (method `"global"`), or to the method for each analysis.

- method:

  - `"simultaneous"` (default): solves `raw = A b` exactly using the
    pairwise 2x2 tables. No joint distribution is needed.

  - `"published"`: the proportional approximation of Rasmussen et al.
    (2022), eq. 16:
    `b[i] = raw[i]^2 / (raw[i] + sum_{k != i} A[i, k] raw[k])`, for
    crude estimates only. Provided for reproduction and comparison; it
    can mask incompatible inputs, and its value is undefined when the
    denominator is zero.

  - `"global"`: fits the maximum-entropy joint distribution
    ([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md);
    the "iterative" model) and solves the equations including any
    interactions. Unknown pairs are left unconstrained. Without
    interactions and unknown pairs, it equals `"simultaneous"`. For more
    than about 20 diseases, pass `backend = "sampled"` (through `...`,
    or fit the joint with
    [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
    and pass it as `joint`).

- joint:

  Optional
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result for `method = "global"`. It is checked against the model
  (diseases, probabilities, associations and three-way terms);
  impact-only changes do not require a refit.

- warn:

  Logical: warn when adjusted impacts change sign or are not finite?

- feasibility:

  For the pairwise methods: `"screen"` (default; triple screen), `"lp"`
  (exact check, needs `lpSolve`) or `"none"`.

## Value

A `cm_result` with elements:

- `adjusted`: raw and adjusted impacts, relative change and estimand;

- `totals`: the sum of `p_i * raw_i` (the naive aggregate), the adjusted
  aggregate and its interaction part;

- `contributions`: per disease, main and interaction contributions,
  their total and share;

- `diagnostics`: reconstruction residual, rank and condition number of
  `A`, sign changes and the feasibility check;

- `conflation` (`A` and `offset`), `interactions` (matrix of `delta`),
  `joint_pairs` (matrix of `P(j and k)`), `model`, `joint`, `method`,
  `label` and `units`. For a
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
  object, a named list of results (class `cm_results`).

## Details

`raw = A %*% b + offset`,

where `A` depends on the disease probabilities and associations (and on
each estimate's estimand), and `offset` holds the contribution of any
pairwise interactions. Results are in the units of the impacts supplied.

## Estimands and the conflation matrix

For a crude estimate of disease `i`, `A[i, k] = P(k | i) - P(k | not i)`
(the excess probability; Rasmussen et al. 2022, eq. 14). For a
coefficient from an additive regression adjusted for the diseases `S`
(`estimand = "adjusted_linear"` in
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)),
`A[i, k]` is the coefficient of `D_i` in the population linear
projection of `D_k` on `D_i` and `D_S`, and `A[i, S] = 0`. These
coefficients depend only on the pairwise tables. They assume the
probabilities and associations describe the population the estimates
come from. An adjustment set that is collinear with the disease is not
identifiable and is rejected.

With interactions `delta` (global method), the offset of disease `i` is
the corresponding regression coefficient of the interaction burden
`sum_{j<k} delta[j, k] D_j D_k`, which needs the joint distribution of
disease combinations.

## Aggregate and contributions

The expected aggregate impact per animal is
`sum_i p_i b_i + sum_{j<k} P(j and k) delta[j, k]`. Each disease's
contribution is its own term plus half of each interaction term it is
involved in (the closed-form Shapley value); contributions add up to the
aggregate. Shares are `NA` when the aggregate is zero.

## Number of diseases

The `"simultaneous"` and `"published"` methods use the pairwise tables
only, so they work for any number of diseases (the triple screen checks
n(n-1)(n-2)/6 triples). The `"global"` method needs the joint
distribution: the exact backend of
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
enumerates 2^n combinations (about 20 diseases at most); the sampled
backend fits the same model by Monte Carlo for more. The exact LP
feasibility check is limited to 14 diseases.

## Feasibility

Pairwise tables can each be valid while no population has all of them
(an invertible `A` does not mean the inputs are feasible). By default
the pairwise methods screen every triple of diseases (a necessary
condition; see
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md));
`feasibility = "lp"` runs the exact check. The global method fits the
joint distribution, and stops if it does not converge.

## Examples

``` r
res <- deconflate(example_supplement())
res$adjusted
#>   disease raw adjusted      change estimand adjusted_for
#> 1      d1 2.5 2.143250 -0.14269981    crude         <NA>
#> 2      d2 5.0 3.387080 -0.32258408    crude         <NA>
#> 3      d3 7.5 6.934209 -0.07543874    crude         <NA>
res$totals
#>   raw_sum adjusted_total interaction_total
#> 1     2.5       2.109229                 0

# A coefficient adjusted for one other disease
m <- example_supplement()
m$impacts <- cm_impacts(c("d1", "d2", "d3"), c(2.2, 5, 7.5),
                        estimand = c("adjusted_linear", "crude", "crude"),
                        adjusted_for = c("d2", NA, NA), units = "%")
deconflate(m)$adjusted
#>   disease raw adjusted      change        estimand adjusted_for
#> 1      d1 2.2 2.354369  0.07016762 adjusted_linear           d2
#> 2      d2 5.0 3.370814 -0.32583725           crude         <NA>
#> 3      d3 7.5 6.936927 -0.07507646           crude         <NA>
```
