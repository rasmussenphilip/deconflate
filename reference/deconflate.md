# Adjust raw impact estimates for comorbidity

Raw impact estimates (comparisons of animals with and without a disease)
are treated as conflations of the disease's own impact and the impacts
of associated diseases. `deconflate()` removes the part of each estimate
that belongs to other diseases, for one impact table at a time: additive
impacts (the default), or event impacts (`event_model = TRUE`). When
inputs have distributions, it also runs draws (`n_draws`) and reports
intervals.

## Usage

``` r
deconflate(model, ...)

# S3 method for class 'cm_model'
deconflate(
  model,
  method = c("auto", "simultaneous", "global"),
  event_model = FALSE,
  overall_risk = NULL,
  n_draws = 1000,
  seed = NULL,
  sampling = c("random", "lhs"),
  lhs_replicates = 10L,
  joint = NULL,
  warn = TRUE,
  feasibility = c("screen", "lp", "none"),
  ...
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  (e.g. from
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md))
  with at least one association.

- ...:

  Passed to
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  (e.g. `backend = "sampled"`).

- method:

  `"auto"` (default), `"simultaneous"` or `"global"` (see Methods). Not
  used for event impacts, which always use the snapshot model.

- event_model:

  `FALSE` (default) for additive impacts; `TRUE` for event impacts (an
  impact table with a `measure` column).

- overall_risk:

  With `event_model = TRUE` (required): the overall risk of the event in
  the population over the period, as a proportion (e.g. `0.25`), or a
  distribution (e.g. `dist_beta(250, 750)`).

- n_draws:

  Number of draws for the uncertainty (default 1000; 0 for point
  estimates only; otherwise at least 2). Used only when some input has a
  distribution.

- seed:

  Optional random seed for the draws and for a sampled joint
  distribution (one is chosen and stored in the result otherwise).

- sampling:

  `"random"` (default) or `"lhs"` (Latin hypercube, in `lhs_replicates`
  independent blocks; the Monte Carlo error is then estimated from the
  block means).

- lhs_replicates:

  Number of Latin hypercube blocks (at least 2; at most `n_draws / 2`
  are used).

- joint:

  Optional
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result for the global method or the event model. It is checked against
  the model (diseases, probabilities, associations and three-way terms).

- warn:

  Logical: warn when adjusted impacts change sign or are not finite?

- feasibility:

  For the simultaneous method: `"screen"` (default; triple screen),
  `"lp"` (exact check, needs `lpSolve`) or `"none"`.

## Value

For additive impacts, a `cm_result` with elements:

- `adjusted`: raw and adjusted impacts, relative change and estimand
  (with draws, also `lower` and `upper`);

- `totals`: the sum of `p_i * raw_i` (the naive aggregate), the adjusted
  aggregate and its interaction part;

- `contributions`: per disease, main and interaction contributions,
  their total and share;

- `diagnostics`: reconstruction residual, rank and condition number of
  `A`, sign changes and the feasibility check;

- `unknown_pairs`: pairs without an association and the odds ratios the
  global fit gave them;

- `draws`: the uncertainty (summary, accepted draws, rejections, seed),
  or `NULL`;

- `notes`, `conflation` (`A` and `offset`), `interactions` (matrix of
  `delta`), `joint_pairs` (matrix of `P(j and k)`), `model`, `joint`,
  `method`, `label` and `units`.

For event impacts, a `cm_event_result` with `adjusted` (raw estimates
and adjusted hazard ratios), `attributable` (`summary`: overall,
disease-free and attributable risk and the attributable fraction;
`by_disease`: the Shapley allocation; `baseline_hazard`), `diagnostics`,
`unknown_pairs`, `draws`, `notes`, `joint`, `model`, `method` and
`overall_risk`. Risks are proportions of animals with the event during
the period.

## Additive impacts

For one vector of additive impacts `b` (in any units), the raw estimates
satisfy

`raw = A %*% b + offset`,

where `A` depends on the disease probabilities and associations (and on
each estimate's estimand), and `offset` holds the contribution of any
pairwise interactions. Results are in the units of the impacts supplied.

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

With interactions `delta`, the offset of disease `i` is the
corresponding regression coefficient of the interaction burden
`sum_{j<k} delta[j, k] D_j D_k`, which needs the joint distribution of
disease combinations.

The expected aggregate impact per animal is
`sum_i p_i b_i + sum_{j<k} P(j and k) delta[j, k]`. Each disease's
contribution is its own term plus half of each interaction term it is
involved in (the closed-form Shapley value); contributions add up to the
aggregate. Shares are `NA` when the aggregate is zero.

## Methods

- The simultaneous method solves `raw = A b` exactly from the pairwise
  2x2 tables. No joint distribution is needed, so it works for any
  number of diseases.

- The global method fits the maximum-entropy distribution of disease
  combinations
  ([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md);
  the "iterative" model) and solves the same equations, including any
  interactions. Pairs without an association are unknown: the fit fills
  them in from the other associations.

The two give the same results when every pair has an association, there
are no interactions and no three-way terms. `method = "auto"` (the
default) uses the simultaneous method then, and the global method
otherwise, with a note giving the reason; `"simultaneous"` behaves the
same way (the global method is used when it is needed), and `"global"`
always uses the global method. The global method's exact backend
enumerates 2^n combinations; with more than 20 diseases the sampled
backend is used (see
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)).

The proportional approximation of Rasmussen et al. (2022, eq. 16) is not
a method here: it is kept for comparison and reproduction in
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
and the `reproduce_*()` functions.

## Event impacts

With `event_model = TRUE`, the impact table holds event impacts (hazard
ratios, rate ratios, risk ratios, odds ratios or risk differences; see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)),
which are adjusted with the snapshot hazard model. Within the period of
the overall risk, an animal with disease combination `d` has a constant
hazard `h0 * exp(sum_i beta_i d_i)`, so its risk over the period is
`R(d) = 1 - exp(-h0 * exp(sum_i beta_i d_i))`; `h0` is set so that the
population risk, averaged over the joint distribution of disease
combinations, equals `overall_risk`. The `beta`s are solved so that each
raw estimate is reproduced:

- a hazard ratio or rate ratio by the ratio of the average hazard
  multipliers among animals with and without the disease;

- a risk ratio, odds ratio or risk difference by the ratio, odds ratio
  or difference of the average risks `R` among animals with and without
  the disease;

at the start of the period (crude), or within strata of the estimate's
adjustment set, combined with Mantel-Haenszel-type weights (stratified).
The adjusted hazard ratios `exp(beta)` are each disease's own hazard
multiplier. The risk attributable to disease is `overall_risk` minus the
disease-free risk `1 - exp(-h0)`; an animal's risk cannot exceed 1, so
it is smaller than the sum of per-disease excess risks when diseases
co-occur. It is allocated to diseases by Shapley values over disease
combinations
([`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)).

What the snapshot model is not: a Cox hazard ratio estimated over
follow-up is not, in general, the snapshot ratio, because animals with
high hazards leave first and the mixture of disease combinations among
survivors changes. The model is exact for its own estimands and a
reasonable approximation when follow-up is short relative to the hazards
or the diseases are rare. Risk-based estimates (risk ratios, odds
ratios, risk differences) must refer to the same period as
`overall_risk`.

## Uncertainty

Inputs with a distribution (the `dist` columns of the input tables, or
`distributions` in
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md);
and `overall_risk` given as a distribution) are drawn `n_draws` times;
each draw is adjusted in the same way as the central estimate. The
central estimate uses the point values (`value`, and the mean of an
`overall_risk` distribution); the draws give 95% intervals (2.5% and
97.5% quantiles), means, Monte Carlo standard errors and stability
checks (`$draws$summary`). Draws whose values cannot hold together (e.g.
associations that no population can have at once) are rejected and
counted (`$draws$rejections`), not replaced; notes in `$notes` report a
high rejection share and limited Monte Carlo precision. Without
distributions, no draws are run.

## Feasibility

Pairwise tables can each be valid while no population has all of them
(an invertible `A` does not mean the inputs are feasible). The
simultaneous method screens every triple of diseases (a necessary
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

# Uncertain inputs
m2 <- cm_model(example_supplement(), example_supplement()$impacts,
               distributions = list("impact:d1" = dist_normal(2.5, 0.5),
                                    "assoc:d2:d3" = dist_lognormal_ci(3, 2, 4.5)))
deconflate(m2, n_draws = 200, seed = 1)
#> <cm_result> method: simultaneous; yield [%]
#> 
#>  disease raw adjusted lower upper   change
#>       d1 2.5    2.143 1.321 3.053 -0.14270
#>       d2 5.0    3.387 2.749 3.843 -0.32258
#>       d3 7.5    6.934 6.832 7.062 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109 (95% interval 2.002 to 2.236)
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).

# Event impacts (culling): hazard ratios and a risk ratio
cull <- cm_model(example_supplement(),
                 cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.2), measure = c("HR", "HR", "RR"),
                            estimand = "snapshot_crude", label = "culling"))
deconflate(cull, event_model = TRUE, overall_risk = 0.25)
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure raw adjusted_hr
#>       d1      HR 1.5       1.382
#>       d2      HR 2.0       1.908
#>       d3      RR 1.2       1.098
#> 
#> Overall risk 0.25; disease-free risk 0.2148; attributable to disease 0.03524 (14.1% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable  share
#>       d1     0.007366 0.2091
#>       d2     0.023888 0.6780
#>       d3     0.003981 0.1130
#> 
#> Diagnostics: residual 5.27e-16, condition number 1.72, sign changes 0
#> 
#> Notes:
#> * No input has a distribution, so no draws were run: the results are point
#>   estimates.
```
