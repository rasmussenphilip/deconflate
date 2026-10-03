# Loss functions for Shapley attribution

Builders for the `loss` argument of
[`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md).

## Usage

``` r
loss_additive(main, terms = NULL)

loss_multiplicative(main)
```

## Arguments

- main:

  Named numeric vector of per-disease impacts.

- terms:

  Optional list of terms, each a list with `diseases` (character vector
  of two or more ids) and `value`.

## Value

A function of a named 0/1 vector.

## Details

- `loss_additive()`: `sum_i m[i] x[i]` plus higher-order terms, each
  applying only when all its diseases are present.

- `loss_multiplicative()`: `1 - prod_i (1 - m[i])^x[i]`, i.e.
  proportional impacts that compound.

## Examples

``` r
f <- loss_additive(c(a = 0.02, b = 0.03),
                   terms = list(list(diseases = c("a", "b"), value = 0.01)))
f(c(a = 1, b = 1))
#> [1] 0.06
```
