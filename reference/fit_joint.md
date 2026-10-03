# Fit the maximum-entropy distribution of disease combinations

The global model needs the probability of every combination of diseases.
Pairwise tables do not determine it (for three or more diseases there
are 2^n - 1 free probabilities and n(n+1)/2 constraints), so the package
uses the maximum-entropy distribution that matches the disease
probabilities and every specified pairwise table: a log-linear model
with main effects and the two-way terms of the constrained pairs, plus
any three-way terms from
[`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md).
This is an assumption, not something the pairwise evidence identifies.
Pairs whose measure is `"unknown"` are not constrained; pairs set to
`"independent"` are constrained to an odds ratio of 1.

## Usage

``` r
fit_joint(
  model,
  tol = 1e-10,
  max_iter = 10000L,
  max_diseases = 20L,
  backend = c("exact", "sampled"),
  n_samples = 50000L,
  n_chains = 1000L,
  burn_in = 50L,
  fit_iter = 300L,
  calibrate = TRUE,
  seed = NULL
)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- tol:

  Convergence tolerance of the IPF: maximum absolute deviation of any
  constrained cell probability, recomputed on the returned distribution.

- max_iter:

  Maximum number of full IPF sweeps.

- max_diseases:

  Safety limit on the number of diseases for the exact backend (the
  table has 2^n cells).

- backend:

  `"exact"` (default) or `"sampled"`.

- n_samples:

  Sampled backend: number of combinations drawn after calibration.

- n_chains:

  Sampled backend: number of parallel Gibbs chains.

- burn_in:

  Sampled backend: Gibbs sweeps discarded before calibration and before
  the final draws.

- fit_iter:

  Sampled backend: maximum number of moment-matching iterations.

- calibrate:

  Sampled backend: adjust the sample weights to match the targets
  exactly?

- seed:

  Optional random seed (sampled backend).

## Value

A `cm_joint` object: `cells` (0/1 matrix of combinations), `prob`
(probability of each combination), `converged`, `iterations`,
`max_residual`, `backend`, the `targets` it was fitted to (used by
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
to check that a supplied joint matches the model) and, for the sampled
backend, `diagnostics`.

## Exact backend

`backend = "exact"` enumerates all 2^n combinations. Starting from
independence (the product of marginal probabilities; slide 12 of the
2023 seminar), with the three-way terms put into the starting
distribution, iterative proportional fitting (IPF) rescales the
combination probabilities until every pairwise table and marginal is
matched. IPF keeps every log-linear term of its start that the
constraints do not fix, so the result is the maximum-entropy
distribution with the requested three-way terms. Enumeration limits this
backend to about 20 diseases.

## Sampled backend

`backend = "sampled"` fits the same log-linear model without enumerating
combinations, for larger numbers of diseases:

1.  Calibration: the model's parameters (a main effect per disease and
    an interaction per constrained pair; three-way terms fixed at the
    log of their ratios) are fitted by Monte Carlo moment matching.
    Gibbs sampling from many parallel chains estimates the probabilities
    and pairwise tables implied by the current parameters, and the
    parameters are updated (by the differences in logits and log odds
    ratios between the targets and the estimates) until the estimates
    match the targets within Monte Carlo error. The interaction
    parameters are conditional log odds ratios; they are not set equal
    to the marginal log odds ratios of the inputs (those are only the
    starting values).

2.  Sampling: `n_samples` combinations are drawn from the calibrated
    model.

3.  Optionally (`calibrate = TRUE`, the default), the weights of the
    sampled combinations are adjusted by IPF so that the disease
    probabilities and pairwise tables match the targets exactly on the
    sample (raking). The remaining Monte Carlo error then affects only
    higher-order structure.

The result has the same form as the exact one (`cells` are the distinct
sampled combinations, `prob` their weights), so the global method,
interaction offsets, Shapley allocation, the snapshot hazard model and
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
work with it. Its `diagnostics` give the constraint residuals of the raw
sample and after raking, R-hat and effective sample sizes of the fitted
probabilities across chains, and the Monte Carlo errors. Results are
approximations with Monte Carlo error: compare runs with different
seeds, and with the exact backend where n allows.

## Convergence

If the pairwise tables are jointly infeasible, the exact IPF does not
converge. Non-convergence after `max_iter` sweeps suggests, but does not
prove, infeasibility (see
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)).
For the sampled backend, a calibration that does not converge is
reported as unresolved, not as evidence that the inputs are infeasible.

## Examples

``` r
j <- fit_joint(example_supplement())
j
#> <cm_joint> 3 diseases, 8 combinations (exact backend)
#>   Converged: TRUE after 8 sweeps (max residual 2.48e-11)
#>   Constrained pairs: 3
# \donttest{
js <- fit_joint(example_supplement(), backend = "sampled", n_samples = 20000, seed = 1)
js$diagnostics$summary
#>   n_samples n_chains n_unique fit_iterations fit_converged residual_sample
#> 1     20000     1000        8             17          TRUE     0.001220607
#>   residual_calibrated max_rhat  min_ess    max_mcse
#> 1        5.075046e-11 1.004167 17279.93 0.003048884
# }
```
