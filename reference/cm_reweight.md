# Reweight Monte Carlo draws (importance sampling)

Reweights accepted draws from
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
so that summaries reflect a different input distribution, without
re-running the adjustment. Weights are self-normalised, so densities
need only be known up to a constant.
[`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
builds the log-ratio from distributions and checks supports; use
`cm_reweight()` directly for custom targets.

## Usage

``` r
cm_reweight(mc, log_ratio)
```

## Arguments

- mc:

  A `cm_mc` object.

- log_ratio:

  A function taking `mc$params` (one row per accepted draw) and
  returning
  `log(target density / density of the inputs' own distributions)` for
  each row. Importance weights from the original run are kept.

## Value

`mc` with updated `weights` and `ess` (Kish effective sample size,
`1 / sum(w^2)`).

## Details

Reweighting cannot recover rejected draws: if the target puts weight
where draws were rejected, results are conditional on acceptance. Always
check the effective sample size.
