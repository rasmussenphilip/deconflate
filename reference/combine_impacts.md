# Combine impact tables

Stacks several
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
objects (e.g. yield and fertility impacts built separately, or culling
impacts from
[`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md))
into one, and re-validates the result.

## Usage

``` r
combine_impacts(...)
```

## Arguments

- ...:

  `cm_impacts` objects.

## Value

A `cm_impacts` data frame.
