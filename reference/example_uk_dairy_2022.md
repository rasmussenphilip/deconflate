# UK dairy example from Rasmussen et al. (2022)

Thirteen endemic diseases and conditions of UK dairy cattle, with
cow-level prevalence (Table 2), inter-disease odds ratios (Table 3, all
other pairs independent), and impacts on milk yield (Table 4, %
decrease), calving interval (Table 5, % increase) and culling (Table 6,
hazard ratios converted to excess annual culling risk). Use
[`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md)
for the observed means and unit values (Table 1).

## Usage

``` r
example_uk_dairy_2022(
  culling_method = c("or_approx", "proportional_hazards"),
  culling_scale = c("proportion", "absolute"),
  yield_sck = 3.05
)
```

## Arguments

- culling_method:

  Hazard-ratio conversion: `"or_approx"` (published) or
  `"proportional_hazards"`.

- culling_scale:

  `"proportion"` (published) or `"absolute"`; see
  [`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md).

- yield_sck:

  Raw yield impact of subclinical ketosis, in percent.

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
with outcomes `"yield"`, `"fertility"` and `"culling"`. The hazard-ratio
conversion is attached as attribute `"culling_conversion"`, for
[`adjusted_hr()`](https://rasmussenphilip.github.io/deconflate/reference/adjusted_hr.md).

## Details

Disease ids: CO cystic ovary, DA displaced abomasum, DYS dystocia, FAS
fasciolosis, GIN gastrointestinal nematodes, LAM lameness, MAS mastitis,
MET metritis, MF milk fever, NEO neosporosis, PTB paratuberculosis, RP
retained placenta, SCK subclinical ketosis.

## Reproducibility of the published tables

With the defaults (`culling_method = "or_approx"`,
`culling_scale = "proportion"`) and `method = "published"`:

- fertility (calving interval) reproduces Tables 8-10 (disease-free
  value 375.09 days; gap value GBP 101.79 per cow per year);

- culling reproduces the adjusted hazard ratios of Table 8 to within
  0.01 and the gap to within 0.5% (22.55% vs 22.56%; GBP 59.48 vs
  59.27);

- yield does not reproduce Tables 8-10 exactly from the printed Table 4.
  The tables imply different yield inputs, e.g. subclinical ketosis 340
  kg / 8737 kg = 3.89% (note i of Table 4) rather than the printed
  3.05%, and Table 10 implies further differences. Use
  `yield_sck = 100 * 340 / 8737` to apply the 3.89% value.

## Examples

``` r
m <- example_uk_dairy_2022()
eco <- uk_dairy_2022_economics()
res <- deconflate(m, method = "published")
gaps <- productivity_gap(res, eco$observed)
gaps$summary
#>     outcome observed disease_free        gap total_loss_fraction
#> 1     yield     8737   9299.20724 562.207239          0.06045754
#> 2 fertility      401    375.09323  25.906774          0.06906756
#> 3   culling       27     22.54568   4.454317          0.19756853
value_losses(gaps, eco$unit_value, eco$additional)$total
#> [1] 402.2476
adjusted_hr(res, attr(m, "culling_conversion"), method = "published")
#>    disease   hr     excess excess_adjusted hr_adjusted
#> 1       CO 1.00 0.00000000      0.00000000    1.000000
#> 2       DA 3.83 0.31384146      0.21890546    2.671438
#> 3      DYS 1.90 0.14205130      0.11895883    1.591128
#> 4      FAS 1.00 0.00000000      0.00000000    1.000000
#> 5      GIN 1.00 0.00000000      0.00000000    1.000000
#> 6      LAM 3.40 0.25564849      0.22551778    2.999276
#> 7      MAS 2.78 0.21306880      0.16920019    2.207627
#> 8      MET 2.20 0.17385764      0.12260834    1.551490
#> 9       MF 2.50 0.20567043      0.15552270    1.890436
#> 10     NEO 1.60 0.09889327      0.09889327    1.600000
#> 11     PTB 2.40 0.19637306      0.13437261    1.642253
#> 12      RP 1.00 0.00000000      0.00000000    1.000000
#> 13     SCK 2.10 0.15726427      0.09665490    1.290664
```
