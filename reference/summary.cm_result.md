# Summarise an adjustment result

Summarise an adjustment result

## Usage

``` r
# S3 method for class 'cm_result'
summary(object, ...)

# S3 method for class 'cm_event_result'
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
(additive impacts: naive and adjusted aggregate and the reduction; event
impacts: overall, disease-free and attributable risk), `diagnostics`,
`contributions` (see
[`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md))
and `notes`.
