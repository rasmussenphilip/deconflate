# Diagnose unstable Monte Carlo estimates

Lists the estimates that
[`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md)
flags, with the reason and a suggested remedy:

- `"possible_pole"`: with the published approximation, `m^2 / (m + c)`
  has a pole where `m + c = 0` and `m != 0`. When the sampled inputs
  reach both sides of it, the estimate's distribution can have tails so
  heavy that its mean does not exist; no sampling scheme then makes the
  mean converge. Report quantiles, or use the exact (`"simultaneous"`)
  method. (With `c = 0`, e.g. independent diseases, the formula reduces
  to `m` and is not flagged.)

- `"heavy_tail"`: a few extreme draws dominate the variance. Importance
  sampling may help if they come from one region of one input. The exact
  method can also be heavy-tailed when sampled inputs make the
  conflation matrix nearly singular.

- `"imprecise"`: the Monte Carlo standard error is large relative to the
  mean. Use more draws, or Latin hypercube sampling.

- `"non_finite"`: draws gave non-finite results and were rejected.

## Usage

``` r
cm_diagnose(mc, what = c("adjusted", "contribution", "total"))
```

## Arguments

- mc:

  A `cm_mc` object.

- what:

  Passed to
  [`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md).

## Value

A data frame (class `cm_diagnosis`) with one row per flagged estimate:
item, method, stability, detail and suggestion.
