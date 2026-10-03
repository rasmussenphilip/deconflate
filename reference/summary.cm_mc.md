# Summarise Monte Carlo results

Summarise Monte Carlo results

## Usage

``` r
# S3 method for class 'cm_mc'
summary(
  object,
  what = c("adjusted", "loss", "total", "rejections"),
  probs = c(0.025, 0.5, 0.975),
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

- ...:

  Unused.

## Value

A data frame of (weighted) means, SDs, quantiles and the Monte Carlo
standard error of the mean (based on the effective sample size).
