# Screen disease pairs for influential impact interactions

Adds a pairwise interaction of each size in `values` to each disease
pair, one at a time, re-runs the global adjustment and reports the
change in the aggregate and rankings. The joint distribution is fitted
once and reused. Interactions apply to additive impacts only.

## Usage

``` r
screen_interactions(model, values, pairs = NULL, ...)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- values:

  Interaction sizes (same units as the impacts; positive synergistic,
  negative antagonistic).

- pairs:

  Optional character vector of pairs (`"d1:d2"`).

- ...:

  Passed to
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  (e.g. `backend = "sampled"` for many diseases).

## Value

A `cm_screen` data frame as in
[`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).
