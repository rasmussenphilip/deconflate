# Reproducing the published analyses

``` r

library(deconflate)
```

This vignette reproduces the calculations in Rasmussen et al. (2022,
*Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy
Sci.* 107:6945-6970), and explains where the package’s results differ
from the published tables.

## Supplementary File example (2022)

``` r

res <- deconflate(example_supplement(), method = "published")
res$adjusted[, c("disease", "raw", "adjusted")]
#>   disease   raw   adjusted
#> 1      d1 0.025 0.02065001
#> 2      d2 0.050 0.03699294
#> 3      d3 0.075 0.06748473
productivity_gap(res, c(yield = 10000))
#> $summary
#>   outcome observed disease_free      gap total_loss_fraction
#> 1   yield    10000     10215.66 215.6617          0.02111089
#> 
#> $attribution
#>   outcome disease       gap  gap_main gap_interaction
#> 1   yield      d1  21.09535  21.09535               0
#> 2   yield      d2  56.68610  56.68610               0
#> 3   yield      d3 137.88024 137.88024               0
```

The published gap of 10,225 units (20, 61 and 143 units per disease)
comes from rounding the adjusted impacts to 0.02, 0.04 and 0.07 before
computing the gap. Without rounding, the published method gives 10,215.7
units.

## UK dairy example (2022, Tables 8-10)

[`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md)
contains the 13 diseases, the 19 odds ratios and the yield, fertility
and culling inputs of Tables 2-6. Culling hazard ratios are converted to
excess annual culling risk by treating them as odds ratios, as in the
paper (`culling_method = "or_approx"`).

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
```

Adjusted culling hazard ratios (Table 8, eq. 23):

``` r

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

**What reproduces.**

- Fertility matches Tables 9-10: a disease-free calving interval of
  375.09 days, valued at GBP 101.79 per cow per year.
- The adjusted culling hazard ratios match Table 8 to within 0.01.
- The culling gap is within 0.5% of the published value.

**What does not.** Yield cannot be reproduced exactly from the printed
Table 4.

- The tables imply different inputs. For subclinical ketosis, note i of
  Table 4 (340 kg per lactation) gives 340 / 8737 = 3.89%, not the
  printed 3.05%. With 3.89% the adjusted values match Table 8 (SCK 2.66
  vs 2.65, lameness 4.75 vs 4.76).
- Table 10’s yield column implies further small differences in the
  inputs, for example for fasciolosis.

``` r

m389 <- example_uk_dairy_2022(yield_sck = 100 * 340 / 8737)
y <- deconflate(m389, method = "published")$adjusted
y[y$outcome == "yield", c("disease", "raw", "adjusted")]
#>    disease        raw     adjusted
#> 1       CO 0.00000000 0.0000000000
#> 2       DA 0.04040000 0.0227781090
#> 3      DYS 0.04050000 0.0293974969
#> 4      FAS 0.07330000 0.0733000000
#> 5      GIN 0.03280000 0.0328000000
#> 6      LAM 0.05540000 0.0475073442
#> 7      MAS 0.04570000 0.0371265733
#> 8      MET 0.03950000 0.0238990351
#> 9       MF 0.00410000 0.0009302367
#> 10     NEO 0.04200000 0.0420000000
#> 11     PTB 0.05900000 0.0443065737
#> 12      RP 0.07380000 0.0608140196
#> 13     SCK 0.03891496 0.0266151273
```

### Corrections the package makes available

- **Exact solution.** `method = "simultaneous"` solves the additive
  equations instead of approximating them.
- **Hazard ratios.** `culling_method = "proportional_hazards"` converts
  hazard ratios to risks under proportional hazards, instead of treating
  them as odds ratios. This gives larger excess risks.
- **Additive culling gap.** `culling_scale = "absolute"` subtracts the
  excess risk from the observed culling rate. The published approach
  instead treats excess risk as a proportional increase in the rate.

``` r

m2 <- example_uk_dairy_2022(culling_method = "proportional_hazards",
                            culling_scale = "absolute")
eco2 <- uk_dairy_2022_economics("absolute")
res2 <- suppressWarnings(deconflate(m2))
summary(res2, economics = eco2)
#> Comorbidity adjustment (method: simultaneous)
#> 
#>    outcome raw_loss adjusted_loss reduction observed disease_free      gap
#>    culling  0.30580       0.23648    0.2267     0.27    3.352e-02   0.2365
#>  fertility  0.07573       0.06448    0.1485   401.00    3.767e+02  24.2919
#>      yield  0.07168       0.05982    0.1654  8737.00    9.293e+03 555.9154
#>   value
#>  315.79
#>   95.43
#>  168.00
#> 
#> Total value (including additional costs): 650.31
#> 
#> Sign changes (raw impacts smaller than associated diseases imply):
#>   yield: MF
#>   fertility: SCK
```

## Global dairy analysis (2024, Table 5)

Table 5 reports the means of comorbidity-adjusted impacts over Monte
Carlo draws of the inputs in Tables 2-4.
[`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
encodes those distributions. In Tables 2-4, the central value of a PERT
distribution is its mode. Run with many more draws for stable estimates;
200 are used here to keep the vignette fast.

``` r

mc <- cm_monte_carlo(sampler_global_dairy(), 200, method = "published", seed = 2024)
s <- summary(mc)
y <- s[s$outcome == "yield", c("disease", "mean", "q0.025", "q0.975")]
y$mean_pct <- round(100 * y$mean, 2)
y$table5 <- c(CK = 0.03, CM = 1.36, DA = 1.18, DYS = 3.48, LAM = 2.62, MET = 2.87,
              MF = 0.07, OC = 2.59, PTB = 3.37, RP = 2.30, SCK = 7.11, SCM = 5.58)[y$disease]
y[, c("disease", "mean_pct", "table5")]
#>    disease mean_pct table5
#> 2       CK     0.03   0.03
#> 4       CM     1.40   1.36
#> 6       DA     1.12   1.18
#> 8      DYS     3.63   3.48
#> 10     LAM     2.61   2.62
#> 12     MET     2.82   2.87
#> 14      MF     0.07   0.07
#> 16      OC     2.53   2.59
#> 18     PTB     3.25   3.37
#> 20      RP     2.25   2.30
#> 22     SCK     7.10   7.11
#> 24     SCM     5.49   5.58
```

The yield means agree closely with Table 5. For fertility impacts whose
raw distributions extend below zero (e.g. metritis, subclinical ketosis
and subclinical mastitis), the published approximation is unstable, so
their means depend on how such draws are handled.

The simultaneous solution at the input means shows the inconsistency
directly. Several adjusted impacts change sign, while the total burden
stays similar:

``` r

res24 <- suppressWarnings(deconflate(example_global_dairy()))
res24$diagnostics[, c("outcome", "n_sign_changes", "sign_changes")]
#>     outcome n_sign_changes          sign_changes
#> 1     yield              4        CK, CM, DA, MF
#> 2 fertility              5 CK, DA, LAM, SCK, SCM
```
