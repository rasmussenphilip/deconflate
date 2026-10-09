# Summarise an adjustment result

Summarise an adjustment result

## Usage

``` r
# S3 method for class 'cm_result'
summary(object, ...)
```

## Arguments

- object:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- ...:

  Unused.

## Value

A `summary.cm_result` list with `method`, `label`, `units`, `totals`
(naive and adjusted aggregate and the reduction), `diagnostics` and
`contributions`.
