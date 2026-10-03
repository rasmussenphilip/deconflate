# Attribute the aggregate burden to diseases (Shapley allocation)

For outcomes on the proportion scale, the expected proportional loss is

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

A data frame with, per outcome and disease, the main-effect burden, the
disease's share of interaction burden, the total and the fraction of
`L`.

## Details

`L = sum_i m[i] P(i) + sum_{j < k} delta[j, k] P(j and k)`.

Removing any disease involved in a term removes that term. Under this
accounting convention, the Shapley value of each disease is its own term
plus an equal share of every interaction term it is involved in:
`s[i] = m[i] P(i) + 1/2 sum_k delta[i, k] P(i and k)`. The shares sum to
`L`.
