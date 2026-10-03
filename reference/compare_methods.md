# Compare adjustment methods

Runs
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
with several methods and returns adjusted impacts side by side.

## Usage

``` r
compare_methods(model, methods = c("published", "simultaneous", "global"), ...)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- methods:

  Methods to compare.

- ...:

  Passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Value

A data frame with one column of adjusted impacts per method.
