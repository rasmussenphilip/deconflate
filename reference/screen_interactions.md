# Screen disease pairs for influential impact interactions

Adds a pairwise interaction of each size in `values` to each disease
pair, one at a time, re-runs the global adjustment and reports the
change in total burden and rankings. The joint distribution is fitted
once and reused. Use it to see which assumed interactions would matter,
before looking for evidence on them.

## Usage

``` r
screen_interactions(
  model,
  outcome,
  values = c(-0.01, 0.01),
  pairs = NULL,
  economics = NULL
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- outcome:

  Outcome to add interactions to.

- values:

  Interaction sizes (proportion scale; positive synergistic, negative
  antagonistic).

- pairs:

  Optional character vector of pairs (`"d1:d2"`).

- economics:

  Optional economics list; otherwise the burden of `outcome` is the
  metric.

## Value

A `cm_screen` data frame as in
[`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).
