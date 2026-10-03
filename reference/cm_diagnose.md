# Diagnose unstable Monte Carlo estimates

Lists the estimates that
[`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md)
flags as unstable, with the reason and a suggested remedy:

- `"no_mean"`: with the published approximation, `m^2 / (m + c)` has a
  pole where `m + c = 0`. If `m + c` takes both signs across draws, the
  adjusted impact's distribution has tails so heavy that its mean does
  not exist, and no sampling scheme makes it converge. Report the median
  or a trimmed mean, or use the exact (`"simultaneous"`) method.
  (Strictly, with unbounded input distributions such as the normal, the
  pole is always inside the support; the flag is raised when draws
  actually reach it.)

- `"heavy_tail"`: a few extreme draws dominate the variance. If they
  come from an identifiable region of one input, importance sampling
  ([`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md)
  and the `proposal` argument of
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md))
  samples that region more often and down-weights it, which reduces the
  Monte Carlo error.

- `"imprecise"`: the Monte Carlo standard error is large relative to the
  mean. Use more draws, or Latin hypercube sampling
  (`sampling = "lhs"`).

## Usage

``` r
cm_diagnose(mc, what = c("adjusted", "loss", "total"))
```

## Arguments

- mc:

  A `cm_mc` object.

- what:

  Passed to
  [`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md).

## Value

A data frame (class `cm_diagnosis`) with one row per flagged estimate:
outcome, disease, method, stability, detail and suggestion.
