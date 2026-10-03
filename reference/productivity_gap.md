# Productivity gap of one adjusted impact vector (optional helper)

Converts the adjusted aggregate impact of one analysis into the gap
between the observed mean of an outcome and its disease-free value, and
attributes the gap to diseases. This generalises Rasmussen et al.
(2022), eqs. 17-22. It is a separate step from the adjustment, because
it needs to know what the impacts mean.

## Usage

``` r
productivity_gap(
  result,
  observed,
  direction = c("decrease", "increase"),
  effect = c("proportion", "percent", "absolute")
)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- observed:

  Observed mean of the outcome.

- direction:

  `"decrease"` or `"increase"`: what disease does to the outcome.

- effect:

  `"proportion"`, `"percent"` or `"absolute"`: how the impacts relate to
  the outcome.

## Value

A list with `summary` (observed, disease-free value, gap and aggregate)
and `attribution` (gap attributed to each disease, split into main and
interaction parts).

## Details

With aggregate `L` and observed mean `x`:

- `effect = "proportion"` (or `"percent"`, divided by 100 first):
  impacts are proportional changes relative to the disease-free value.
  For an outcome that disease decreases (e.g. yield) the disease-free
  value is `x / (1 - L)`; for one it increases (e.g. calving interval),
  `x / (1 + L)`.

- `effect = "absolute"`: impacts are in the outcome's units; the
  disease-free value is `x + L` (decrease) or `x - L` (increase).

Each disease's part of the gap is its contribution (see
[`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md))
times the same factor, computed directly rather than by dividing by `L`,
so it stays defined when contributions cancel.

## Examples

``` r
res <- deconflate(example_supplement(), method = "published")
productivity_gap(res, observed = 10000, direction = "decrease", effect = "percent")
#> $summary
#>   observed disease_free      gap  aggregate direction  effect
#> 1    10000     10215.66 215.6617 0.02111089  decrease percent
#> 
#> $attribution
#>   disease       gap  gap_main gap_interaction
#> 1      d1  21.09535  21.09535               0
#> 2      d2  56.68610  56.68610               0
#> 3      d3 137.88024 137.88024               0
#> 
```
