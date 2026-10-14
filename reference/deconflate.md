# Adjust raw impact estimates for comorbidity

Raw impact estimates (comparisons of animals with and without a disease)
are conflations of the disease's own impact and the impacts of
associated diseases. `deconflate()` removes the part of each estimate
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
  event_model = FALSE,
  overall_risk = NULL,
  n_draws = 1000,
  seed = NULL,
  sampling = c("random", "lhs"),
  lhs_replicates = 10L,
  pool = c("random", "fixed"),
  joint = NULL,
  warn = TRUE,
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
  (e.g. `backend = "sampled"`). With more than 20 diseases the sampled
  backend is used.

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

- pool:

  `"random"` (default) or `"fixed"`: how several association rows about
  the same quantity are pooled.

- joint:

  Optional
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result. It is checked against the model (diseases, probabilities and
  association rows).

- warn:

  Logical: warn when adjusted impacts change sign or are not finite?

## Value

For additive impacts, a `cm_result` with elements:

- `adjusted`: one row per disease and per interaction: the raw value (of
  its row, when it has exactly one), the adjusted value, the relative
  change and the number of rows (with draws, also `lower` and `upper`);

- `totals`: the sum of `p_d * raw_d` (the naive aggregate; `NA` when a
  disease has several rows), the adjusted aggregate and its interaction
  part;

- `contributions`: per disease, main and interaction contributions,
  their total and share;

- `estimates`: every input row (impacts and associations) with its
  value, fitted value and standardised residual (and, for impact rows,
  the SD used as weight; 0 for an exact row);

- `heterogeneity`: pooled association rows and the fit of the impact
  rows (Q, df, p, I^2, tau^2), or `NULL`;

- `diagnostics`: the computation used (`"full model"` or
  `"pairwise shortcut"`), reconstruction residual, rank and condition
  number of the row equations, sign changes and the feasibility check;

- `unknown_pairs`: pairs that no row mentions and the odds ratios the
  fit gave them;

- `draws`: the uncertainty (summary, accepted draws, rejections, seed),
  or `NULL`;

- `notes`, `conflation` (the row equations), `joint_pairs` (matrix of
  `P(j and k)`), `model`, `joint`, `method`, `label` and `units`.

For event impacts, a `cm_event_result` with `adjusted` (adjusted hazard
ratios of diseases and interactions), `attributable` (`summary`:
overall, disease-free and attributable risk and the attributable
fraction; `by_disease`: the Shapley allocation; `baseline_hazard`),
`estimates`, `heterogeneity`, `diagnostics`, `unknown_pairs`, `draws`,
`notes`, `joint`, `model`, `method` and `overall_risk`. Risks are
proportions of animals with the event during the period.

## The model

The package has one model. The joint distribution of disease
combinations is the maximum-entropy distribution that matches the
disease probabilities and every association row
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md));
pairs that no row mentions are unknown and are filled in by the fit
(listed in `$unknown_pairs`). Every impact row is then an equation in
the unknown impacts: an impact for each disease and an interaction for
each combination that appears in the impact table (see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
for what each row means).

When the impacts are additive and the inputs are only crude pair
associations (one for every pair) and impact rows about single diseases,
the pairwise 2x2 tables give exactly the same answer (the simultaneous
solution of Rasmussen et al. 2022, eq. 14, generalised to adjusted
estimates), so the package uses them ("pairwise shortcut", recorded in
the result), and fits the joint distribution once to confirm that the
associations can hold together (up to 20 diseases). The draws then use
the shortcut too.

The proportional approximation of Rasmussen et al. (2022, eq. 16) is not
part of the model: it is kept for comparison and reproduction in
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
and the `reproduce_*()` functions.

## Additive impacts

With impacts `b` and interactions `delta`, the outcome of an animal with
disease combination `c` is
`y(c) = sum_d b_d c_d + sum_T delta_T all_T(c)`. Each impact row is a
linear function of `b` and `delta`: for a crude row, the difference in
mean outcome between animals with all the row's diseases and animals
with none of them (within `among`); for an adjusted row, a regression
coefficient. Results are in the units of the impacts.

The expected aggregate impact per animal is
`sum_d p_d b_d + sum_T P(T) delta_T`. Each disease's contribution is its
own term plus an equal share of each interaction it is part of (the
closed-form Shapley value); contributions add up to the aggregate.
Shares are `NA` when the aggregate is zero.

## Several rows per item

Impact rows are fitted together: with one row per unknown the equations
are solved exactly; with more rows, by weighted least squares with
weights 1/SD^2 (SDs from each row's distribution; a row without a
distribution is exact and is reproduced exactly). The result reports
each row's fitted value and standardised residual (`$estimates`) and the
heterogeneity statistic Q (`$heterogeneity`). Two exact rows that
contradict each other are an error.

Association rows that describe the same quantity are pooled on the log
odds ratio scale of their 2x2 table (the logit of a prevalence, the log
of an `OR_ratio`, or, for rows of different kinds about the same
quantity, the logit of its probability), inverse-variance weighted, with
random effects (DerSimonian-Laird) by default and `pool = "fixed"` as an
option; Q, I^2 and tau^2 are in `$heterogeneity`. A row without a
distribution is exact and wins. In the draws, each row is drawn from its
own distribution and pooled with the central weights, and a
between-study term keeps the pooled variance at its random-effects
value.

## Event impacts

With `event_model = TRUE`, the impact table holds event impacts (hazard
ratios, rate ratios, risk ratios, odds ratios or risk differences; see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)),
which are adjusted with the snapshot hazard model. Within the period of
the overall risk, an animal with disease combination `c` has a constant
hazard `h0 * exp(sum_d beta_d c_d + sum_T gamma_T all_T(c))`, so its
risk over the period is `R(c) = 1 - exp(-h0 * exp(...))`; `h0` is set so
that the population risk equals `overall_risk`. The `beta`s and `gamma`s
(interactions on the log-hazard scale) are solved so that each raw
estimate is reproduced: a hazard or rate ratio by the ratio of the
average hazard multipliers of the two groups the row compares; a risk
ratio, odds ratio or risk difference by those of their average risks;
crude, or within strata of the row's adjustment set (combined with
Mantel-Haenszel-type weights). The adjusted hazard ratios `exp(beta)`
are each disease's own hazard multiplier. The risk attributable to
disease is `overall_risk` minus the disease-free risk `1 - exp(-h0)`,
allocated to diseases by Shapley values over disease combinations
([`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)).

What the snapshot model is not: a Cox hazard ratio estimated over
follow-up is not, in general, the snapshot ratio, because animals with
high hazards leave first and the mixture of disease combinations among
survivors changes. The model is exact for its own comparisons and a
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
checks (`$draws$summary`). Draws whose values cannot hold together are
rejected and counted (`$draws$rejections`), not replaced; notes in
`$notes` report a high rejection share and limited Monte Carlo
precision. Without distributions, no draws are run.

## Errors and warnings

Errors and warnings are condition objects with classes (see
[deconflate-package](https://rasmussenphilip.github.io/deconflate/reference/deconflate-package.md)).
Those about particular inputs also hold them as fields, so that code can
use them without reading the message:

- `table` and `row` (or `rows`): the input rows concerned, e.g. an
  impact row that is not identifiable, or exact rows that contradict
  each other (with their `residuals` for impact rows);

- `terms`: diseases and interactions, e.g. those the impact rows do not
  determine (`deconflate_singular`), or those whose adjusted impact
  changes sign (the warning `deconflate_sign_change`);

- `triples`: the conflicting triples when the crude pair rows cannot all
  hold (`deconflate_infeasible`), as in
  [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md);

- `pairs`: the pairs without a crude pair row, for the calculations that
  need every pair (`deconflate_unknown_pairs`);

- `residual`: how far a fit was from converging
  (`deconflate_nonconvergence`).

## Examples

``` r
res <- deconflate(example_supplement())
res$adjusted
#>   disease    kind raw adjusted      change rows
#> 1      d1 disease 2.5 2.143250 -0.14269981    1
#> 2      d2 disease 5.0 3.387080 -0.32258408    1
#> 3      d3 disease 7.5 6.934209 -0.07543874    1
res$totals
#>   raw_sum adjusted_total interaction_total
#> 1     2.5       2.109229                 0

# An estimate adjusted for one other disease, and a second estimate for d1
m <- example_supplement()
m$impacts <- cm_impacts(c("d1", "d1", "d2", "d3"), c(2.2, 2.6, 5, 7.5),
                        adjusted_for = c("d2", NA, NA, NA), units = "%")
m$distributions <- list("impacts:1" = dist_normal(2.2, 0.4),
                        "impacts:2" = dist_normal(2.6, 0.5))
deconflate(m, n_draws = 0)$estimates
#>          table row    item measure among adjusted_for study value  sd   fitted
#> 1      impacts   1      d1    <NA>  <NA>           d2  <NA>   2.2 0.4 2.157380
#> 2      impacts   2      d1    <NA>  <NA>         <NA>  <NA>   2.6 0.5 2.667119
#> 3      impacts   3      d2    <NA>  <NA>         <NA>  <NA>   5.0 0.0 5.000000
#> 4      impacts   4      d3    <NA>  <NA>         <NA>  <NA>   7.5 0.0 7.500000
#> 5 associations   1 d1 ~ d2      OR  <NA>         <NA>  <NA>   2.0  NA 2.000000
#> 6 associations   2 d1 ~ d3      OR  <NA>         <NA>  <NA>   1.0  NA 1.000000
#> 7 associations   3 d2 ~ d3      OR  <NA>         <NA>  <NA>   3.0  NA 3.000000
#>   weight std_residual
#> 1     NA    0.1065495
#> 2     NA   -0.1342381
#> 3     NA           NA
#> 4     NA           NA
#> 5      1           NA
#> 6      1           NA
#> 7      1           NA

# Uncertain inputs
m2 <- cm_model(example_supplement(), example_supplement()$impacts,
               distributions = list("impact:d1" = dist_normal(2.5, 0.5),
                                    "assoc:d2:d3" = dist_lognormal_ci(3, 2, 4.5)))
deconflate(m2, n_draws = 200, seed = 1)
#> <cm_result> pairwise shortcut; yield [%]
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
                            label = "culling"))
deconflate(cull, event_model = TRUE, overall_risk = 0.25)
#> <cm_event_result> culling; snapshot
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
#> Diagnostics: residual 6.76e-12, condition number 1.72, sign changes 0
#> 
#> Notes:
#> * No input has a distribution, so no draws were run: the results are point
#>   estimates.

# Errors hold the inputs they are about
bad <- cm_model(cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                              cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05))),
                cm_impacts(c("a", "b", "c"), c(1, 1, 1)))
e <- tryCatch(deconflate(bad), deconflate_infeasible = function(e) e)
e$triples
#>   disease1 disease2 disease3      gap
#> 1        a        b        c 0.225884
```
