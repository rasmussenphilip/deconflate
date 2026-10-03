# Economic inputs for the UK dairy example (Rasmussen et al. 2022, Table 1)

Observed means and unit values for
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
and
[`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md):
milk yield 8737 kg/cow/year valued at GBP 0.3022/kg; calving interval
401 days, each day valued at lifetime daily yield (13 kg) times the milk
price; culling rate 27% per year valued at the replacement price (GBP
1335.36); and veterinary expenditure of GBP 71.09 per cow per year added
as a lump sum.

## Usage

``` r
uk_dairy_2022_economics(culling_scale = c("proportion", "absolute"))
```

## Arguments

- culling_scale:

  Must match the `culling_scale` used in
  [`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md):
  with `"proportion"` the culling rate is in percent (27) and valued per
  percentage point; with `"absolute"` it is a proportion (0.27) valued
  per unit.

## Value

A list with `observed`, `unit_value` and `additional`.
