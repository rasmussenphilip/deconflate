# Summarise Monte Carlo results

Summarises adjusted impacts (or losses) over the accepted draws, and
checks each estimate's stability. When some estimates look unstable, a
message explains why and suggests a remedy (see
[`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)).

## Usage

``` r
# S3 method for class 'cm_mc'
summary(
  object,
  what = c("adjusted", "loss", "total", "rejections"),
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

  `"adjusted"` (adjusted impacts by outcome and disease), `"loss"`
  (monetary losses by outcome and disease), `"total"` (total losses per
  outcome and overall) or `"rejections"`.

- probs:

  Quantiles to report.

- trim:

  Fraction trimmed from each tail for `trimmed_mean`.

- diagnose:

  Logical: print suggestions when estimates look unstable?

- ...:

  Unused.

## Value

A data frame of (weighted) means, SDs, the Monte Carlo standard error of
the mean (`mcse`, based on the effective sample size), quantiles, a
trimmed mean and stability diagnostics:

- `rel_mcse`: `mcse` relative to the absolute mean;

- `tail_share`: the share of the variance contributed by the most
  extreme 1% of draws;

- `stability`: `"ok"`, `"imprecise"` (`rel_mcse` above 5%),
  `"heavy_tail"` (`tail_share` above 60%) or `"no_mean"` (the published
  approximation divides by a quantity that changes sign across draws, so
  the mean does not exist).
