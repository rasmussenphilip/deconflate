# Excess-probability matrix from pairwise tables

Excess-probability matrix from pairwise tables

## Usage

``` r
excess_matrix(model)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Value

An n x n matrix `E` with `E[k, i] = P(k | i) - P(k | not i)`, the excess
probability of disease `k` among animals with disease `i` (`ep_ki` in
Rasmussen et al. 2022, eq. 14). The diagonal is 0. Errors if any pair is
`"unknown"`: pairwise methods need every pair specified.
