# Summarise an adjustment result

Summarise an adjustment result

## Usage

``` r
# S3 method for class 'cm_result'
summary(object, valuation = NULL, ...)
```

## Arguments

- object:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- valuation:

  Optional valuation list (see
  [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)).

- ...:

  Unused.

## Value

A `summary.cm_result` list with `method`, `label`, `units`, `totals`
(naive and adjusted aggregate, their difference, and the gap and value
with `valuation`), `diagnostics` and `contributions`.
