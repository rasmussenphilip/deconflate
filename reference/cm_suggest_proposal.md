# Suggest an importance-sampling proposal for an unstable estimate

Finds the sampled input most associated with the most extreme draws of
one adjusted impact, and builds a defensive mixture proposal for it:
half the input's own distribution and half a uniform distribution over
the range of that input in the extreme draws. Re-running
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
with this proposal samples the region that produces extreme values more
often and down-weights it, which reduces the Monte Carlo error when the
mean exists. Weights are bounded by `1 / (1 - weight)`, so the run
cannot be dominated by a single draw.

## Usage

``` r
cm_suggest_proposal(
  mc,
  outcome,
  disease,
  method = NULL,
  top = 0.02,
  weight = 0.5
)
```

## Arguments

- mc:

  A `cm_mc` object from a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)-based
  run.

- outcome, disease:

  The estimate to stabilise.

- method:

  Adjustment method (default: the first in `mc`).

- top:

  Fraction of draws treated as extreme.

- weight:

  Weight of the uniform component in the mixture.

## Value

A named list with one `cm_dist`, for the `proposal` argument of
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md).
The attribute `"explanation"` describes it.

## Details

It cannot help when the mean does not exist (stability `"no_mean"`; see
[`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)).
