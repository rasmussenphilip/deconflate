# Scenario analysis by reweighting Monte Carlo draws

Replaces the distributions of some inputs with scenario distributions
and reweights the existing draws accordingly (see
[`cm_reweight()`](https://rasmussenphilip.github.io/deconflate/reference/cm_reweight.md)).
Keys are those of `mc$params`, e.g. `"assoc:LAM:SCK"` or `"impact:LAM"`.
Each scenario distribution must lie within the support of the
distribution the input was sampled from.

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
