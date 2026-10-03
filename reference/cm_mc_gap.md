# Productivity gaps over Monte Carlo draws (optional helper)

Applies
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
(and a unit value) to every accepted draw of a Monte Carlo run, and
summarises the gap and its value.

## Usage

``` r
cm_mc_gap(
  mc,
  observed,
  direction = c("decrease", "increase"),
  effect = c("proportion", "percent", "absolute"),
  unit_value = NULL,
  seed = NULL
)
```

## Arguments

- mc:

  A `cm_mc` object.

- observed:

  Observed mean of the outcome (a number or a `cm_dist`, drawn per
  draw).

- direction, effect:

  As in
  [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md).

- unit_value:

  Optional value per unit of gap (a number or a `cm_dist`, drawn per
  draw).

- seed:

  Optional seed for drawing `observed` and `unit_value`.

## Value

A list with `draws` (gap and value by draw, method and disease, and the
totals) and `summary` (weighted mean and quantiles of the total gap and
value per method).
