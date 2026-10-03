# Summarise Monte Carlo results

Summarises adjusted impacts, contributions or aggregates over the
accepted draws, with Monte Carlo standard errors and stability checks.
When some estimates look unstable, a message explains why and suggests
what to do (see
[`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)).

## Usage

``` r
# S3 method for class 'cm_mc'
summary(
  object,
  what = c("adjusted", "contribution", "total", "rejections"),
  probs = c(0.025, 0.5, 0.975),
  trim = 0.05,
  diagnose = TRUE,
  ...
)
```

## Arguments

- object:

  A `cm_mc` object.

- what:

  `"adjusted"` (adjusted impacts by disease), `"contribution"`
  (contributions to the aggregate by disease), `"total"` (naive and
  adjusted aggregate) or `"rejections"` (counts by type).

- probs:

  Quantiles to report.

- trim:

  Fraction trimmed from each tail for `trimmed_mean` (a different
  estimand from the mean).

- diagnose:

  Logical: print suggestions when estimates look unstable?

- ...:

  Unused.

## Value

A data frame with (weighted) means, SDs, the Monte Carlo standard error
of the mean (`mcse`), quantiles, a trimmed mean and stability
diagnostics:

- `mcse`: for independent draws, the self-normalised importance-sampling
  estimate `sqrt(sum(w^2 (x - mean)^2))` (with equal weights,
  `sd / sqrt(n)`); for Latin hypercube runs, the standard error of the
  pooled mean as a ratio estimator over the R replicate blocks, with
  `S_b` the weighted sum and `W_b` the weight of block b:
  `sqrt(sum((S_b - mean * W_b)^2) / (R (R - 1))) / mean(W_b)`. Blocks
  whose draws were all rejected, or have zero weight, count in R with
  `S_b = W_b = 0`. With equal weights and every block present, this is
  the standard deviation of the block means over `sqrt(R)`;

- `rel_mcse`: `mcse` relative to the absolute mean;

- `tail_share`: the share of the variance contributed by the most
  extreme 1% of draws;

- `stability`: `"ok"`, `"imprecise"` (`rel_mcse` above 5%),
  `"insufficient_info"` (precision cannot be assessed: fewer than two
  Latin hypercube blocks with positive weight, or a mean of zero with a
  positive standard error; `mcse` gives the absolute precision where
  available), `"heavy_tail"` (`tail_share` above 60%) or
  `"possible_pole"` (the published approximation's denominator changes
  sign within the sampled inputs in a way the raw impact does not
  explain, so the estimate has a pole inside the input distribution and
  its mean may not exist). A missing precision is never reported as
  `"ok"`.
