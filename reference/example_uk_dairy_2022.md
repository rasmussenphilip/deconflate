# UK dairy example from Rasmussen et al. (2022)

Thirteen endemic diseases and conditions of UK dairy cattle, with
cow-level prevalence (Table 2), inter-disease odds ratios (Table 3), and
the raw impacts of one outcome: milk yield (Table 4, % decrease),
calving interval (Table 5, % increase) or culling (Table 6, hazard
ratios, as event impacts (crude snapshot comparisons, an assumption);
use `deconflate(..., event_model = TRUE, overall_risk = 0.27)`, the
annual culling rate of the paper). Use
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
for the published tables, including the economic values.

## Usage

``` r
example_uk_dairy_2022(
  outcome = c("yield", "fertility", "culling"),
  yield_sck = 3.05
)
```

## Arguments

- outcome:

  `"yield"`, `"fertility"` or `"culling"`.

- yield_sck:

  Raw yield impact of subclinical ketosis, in percent (the printed 3.05,
  or 100 \* 340 / 8737 = 3.89 from note i of Table 4).

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Details

Pairs not in Table 3 are unknown, as for any input in this package: the
model fills them in. The paper treated them as independent;
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
does the same. The paper converted the hazard ratios to excess annual
culling risks by treating them as odds ratios; that historical
calculation is only available inside
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).

Disease ids: CO cystic ovary, DA displaced abomasum, DYS dystocia, FAS
fasciolosis, GIN gastrointestinal nematodes, LAM lameness, MAS mastitis,
MET metritis, MF milk fever, NEO neosporosis, PTB paratuberculosis, RP
retained placenta, SCK subclinical ketosis.

## Examples

``` r
uk <- example_uk_dairy_2022()
deconflate(uk)
#> Warning: Adjusted impacts change sign for MF. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimates, and whether they come from populations with different comorbidity patterns.
#> <cm_result> full model; milk yield loss [% decrease]
#> 
#>  disease  raw adjusted     change
#>       CO 0.00  -0.7247         NA
#>       DA 4.04   1.9547 -5.162e-01
#>      DYS 4.05   2.8034 -3.078e-01
#>      FAS 7.33   7.3300 -4.441e-16
#>      GIN 3.28   3.2800  1.110e-15
#>      LAM 5.54   4.9227 -1.114e-01
#>      MAS 4.57   3.8983 -1.470e-01
#>      MET 3.95   1.8763 -5.250e-01
#>       MF 0.41  -0.6895 -2.682e+00
#>      NEO 4.20   4.2000  0.000e+00
#>      PTB 5.90   4.0784 -3.087e-01
#>       RP 7.38   6.3352 -1.416e-01
#>      SCK 3.05   1.6417 -4.617e-01
#> 
#> Raw sum: 7.168; adjusted total: 5.843
#> Diagnostics: residual 1.78e-15, condition number 2.07, sign changes 1 (MF)
#> Unknown pairs: 59 without an association; the fit gave them odds ratios from 1 to 1.38 (see $unknown_pairs).
#> 
#> Notes:
#> * No input has a distribution, so no draws were run: the results are point
#>   estimates.
```
