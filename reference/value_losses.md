# Value productivity gaps in monetary terms

Multiplies each outcome's gap (from
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md))
by a unit value and returns disease-specific and total losses, as in
Rasmussen et al. (2022), Tables 9-10.

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

  Named numeric vector: value per unit of gap for each outcome. For
  example, with yield in kg the milk price per kg; with calving interval
  in days the value of a day's lost milk; with culling in percentage
  points the replacement cost per percentage point.

- additional:

  Optional named numeric vector of costs added to the total as lump sums
  (e.g. veterinary expenditure per animal).

## Value

A list with `by_disease` (outcome, disease, gap, value), `by_outcome`
(gap and value per outcome), `by_disease_total` (value per disease
summed over outcomes) and `total` (including `additional`).

## Examples

``` r
res <- deconflate(example_supplement(), method = "published")
value_losses(productivity_gap(res, c(yield = 10000)), c(yield = 0.30))
#> $by_disease
#>   outcome disease       gap     value
#> 1   yield      d1  21.09535  6.328605
#> 2   yield      d2  56.68610 17.005829
#> 3   yield      d3 137.88024 41.364071
#> 
#> $by_outcome
#>   outcome      gap    value
#> 1   yield 215.6617 64.69851
#> 
#> $by_disease_total
#>   disease     value
#> 3      d3 41.364071
#> 2      d2 17.005829
#> 1      d1  6.328605
#> 
#> $total
#> [1] 64.69851
#> 
```
