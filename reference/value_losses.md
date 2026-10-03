# Value a productivity gap in monetary terms (optional helper)

Multiplies a gap (from
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md))
by a unit value, as in Rasmussen et al. (2022), Tables 9-10. Lump-sum
costs that are not part of any adjustment (e.g. veterinary expenditure)
can be added with `additional`; they are reported separately and
included in `total`.

## Usage

``` r
value_losses(gap, unit_value, additional = 0)
```

## Arguments

- gap:

  A
  [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
  result.

- unit_value:

  Value per unit of gap (e.g. milk price per kg).

- additional:

  Optional named numeric vector of lump-sum costs.

## Value

A list with `by_disease` (gap and value), `value` (the gap's value),
`additional` and `total`.

## Examples

``` r
res <- deconflate(example_supplement(), method = "published")
value_losses(productivity_gap(res, 10000, "decrease", "percent"), unit_value = 0.30)
#> $by_disease
#>   disease       gap     value
#> 1      d1  21.09535  6.328605
#> 2      d2  56.68610 17.005829
#> 3      d3 137.88024 41.364071
#> 
#> $value
#> [1] 64.69851
#> 
#> $additional
#> [1] 0
#> 
#> $total
#> [1] 64.69851
#> 
```
