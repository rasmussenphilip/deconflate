# Fit the maximum-entropy distribution of disease combinations

Starting from independence (the product of marginal probabilities; slide
12 of the 2023 seminar), iterative proportional fitting (IPF) repeatedly
rescales the 2^n combination probabilities so that every specified
pairwise 2x2 table is matched. From an independence start, IPF converges
to the maximum-entropy distribution subject to these constraints: a
log-linear model with all two-way terms and no explicit three-way or
higher log-linear interaction terms. Three-disease co-occurrence still
departs from independence through the pairwise terms.

## Usage

``` r
fit_joint(model, tol = 1e-10, max_iter = 10000L, max_diseases = 20L)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- tol:

  Convergence tolerance: maximum absolute deviation of any constrained
  cell probability, recomputed on the returned distribution.

- max_iter:

  Maximum number of full IPF sweeps.

- max_diseases:

  Safety limit on the number of diseases (the table has 2^n cells).

## Value

A `cm_joint` object: `cells` (2^n x n 0/1 matrix), `prob` (probability
of each combination), `converged`, `iterations`, `max_residual` and the
`targets` it was fitted to (used by
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
to check that a supplied joint matches the model).

## Details

Three-way terms from
[`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md)
are put into the starting distribution. IPF keeps every log-linear term
of its starting distribution that the constraints do not fix, so the
fitted distribution matches all pairwise tables and has the requested
three-way terms.

Pairs whose measure is `"unknown"` are not constrained. Pairs set to
`"independent"` are constrained to an odds ratio of 1.

The constraint system is underdetermined for three or more diseases
(2^n - 1 free probabilities, n(n+1)/2 constraints), which is why the
maximum-entropy criterion (an assumption, not something the pairwise
evidence identifies) is needed to pick one solution. If the pairwise
tables are jointly infeasible, IPF does not converge; non-convergence
after `max_iter` sweeps suggests, but does not prove, infeasibility (see
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)).
