# Summarise an adjustment result

Summarise an adjustment result

## Usage

``` r
# S3 method for class 'cm_result'
summary(object, economics = NULL, ...)
```

## Arguments

- object:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- economics:

  Optional economics list (see
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)).

- ...:

  Unused.

## Value

A `summary.cm_result` list with `method`, `diagnostics`, `totals`
(expected loss per outcome, raw and adjusted, plus gaps and values if
`economics` is given) and `contributions`
([`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)).
