# Attribute the aggregate to diseases (Shapley allocation)

The expected aggregate impact is
`sum_i p_i b_i + sum_{j < k} delta[j, k] P(j and k)`. Removing any
disease involved in a term removes that term, so each disease's Shapley
value is its own term plus an equal share of every interaction term it
is involved in. The shares add up to the aggregate. This is the
`contributions` table of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Usage

``` r
attribute_burden(result)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

## Value

A data frame with, per disease, the main contribution, the share of
interaction terms, the total and the share of the aggregate (`NA` when
the aggregate is zero).
