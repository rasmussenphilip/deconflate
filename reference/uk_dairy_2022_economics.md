# Economic inputs for the UK dairy example (Rasmussen et al. 2022, Table 1)

Valuation inputs for each analysis of
[`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md),
for
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
and
[`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)
(or the `valuation` argument of
[`summary.cm_result()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_result.md)
and
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)):
milk yield 8737 kg/cow/year valued at GBP 0.3022/kg; calving interval
401 days, each day valued at lifetime daily yield (13 kg) times the milk
price; culling rate 27% per year (impacts as proportional increases, as
in the paper) valued at the replacement price per percentage point (GBP
13.3536). Veterinary expenditure of GBP 71.09 per cow per year is a
separate lump sum.

## Usage

``` r
uk_dairy_2022_economics()
```

## Value

A list with `valuation` (one valuation list per analysis) and
`additional` (lump-sum costs).
