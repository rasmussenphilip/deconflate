# Suggest an importance-sampling proposal for an unstable estimate

Finds the sampled input most associated with the most extreme draws of
one adjusted impact, and builds a defensive mixture proposal for it:
part the input's own distribution and part a uniform distribution over
the range of that input in the extreme draws (within the input's
support). The mixture always covers the input's support and bounds the
weights by `1 / (1 - weight)`. It may reduce the Monte Carlo error when
the extreme values come from that region and the mean exists; it does
not guarantee a finite variance or better precision, so compare standard
errors.

## Usage

``` r
cm_suggest_proposal(mc, disease, method = NULL, top = 0.02, weight = 0.5)
```

## Arguments

- mc:

  A `cm_mc` object from a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)-based
  run.

- disease:

  The disease whose adjusted impact is unstable.

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
