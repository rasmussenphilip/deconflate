# Reweight Monte Carlo draws (importance sampling)

Reweights accepted draws from
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
so that summaries reflect a different input distribution, without
re-running the adjustment. Weights are self-normalised, so densities
need only be known up to a constant.
[`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
builds the log-ratio from distributions; use `cm_reweight()` directly
for custom targets.

## Usage

``` r
cm_reweight(mc, log_ratio)
```

## Arguments

- mc:

  A `cm_mc` object.

- log_ratio:

  A function taking `mc$params` (one row per accepted draw) and
  returning `log(target density / proposal density)` for each row.

## Value

`mc` with updated `weights` and `ess` (effective sample size,
`1 / sum(w^2)`). Importance weights from the original run (see the
`proposal` argument of
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md))
are kept.

## Details

The sampling (proposal) distribution must cover the scenario: give the
sampler wider spread than the base case (e.g.
[`dist_mixture()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)).
Reweighting cannot recover rejected draws: if the scenario puts weight
where draws were infeasible, results are conditional on feasibility.
Always check the effective sample size.
