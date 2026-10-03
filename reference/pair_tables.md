# Pairwise 2x2 tables for every disease pair

Pairwise 2x2 tables for every disease pair

## Usage

``` r
pair_tables(model)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Value

A data frame with one row per unordered pair: the measure used, whether
it was specified or defaulted, the joint probability `p11` (`NA` for
unknown pairs), the implied odds ratio, and the excess probabilities
`ep_2_given_1 = P(d2 | d1) - P(d2 | not d1)` and `ep_1_given_2`.
