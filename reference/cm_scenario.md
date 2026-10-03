# Scenario analysis by reweighting Monte Carlo draws

Replaces the distributions of some inputs with scenario distributions
and reweights the existing draws accordingly (see
[`cm_reweight()`](https://rasmussenphilip.github.io/deconflate/reference/cm_reweight.md)).
Keys are those of `mc$params`, e.g. `"assoc:LAM:SCK"` or
`"impact:yield:LAM"`.

## Usage

``` r
cm_scenario(mc, changes)
```

## Arguments

- mc:

  A `cm_mc` object from a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)-based
  run.

- changes:

  Named list of `cm_dist` objects (the scenario distributions).

## Value

The reweighted `mc` (with `ess`).

## Details

Inputs drawn with an outcome correlation (copula) cannot be reweighted
individually, because changing a marginal also changes the copula
coordinates; re-run
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
for such scenarios.
