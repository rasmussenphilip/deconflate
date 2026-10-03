# Reproducing the published analyses

``` r

library(deconflate)
```

This vignette reproduces the calculations in Rasmussen et al. (2022,
*Prev. Vet. Med.* 203:105617) and Rasmussen et al. (2024, *J. Dairy
Sci.* 107:6945-6970), and explains where the package’s results differ
from the published tables.
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
and
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
use the published (eq. 16) method and the hazard-ratio conversions of
the papers. Those conversions are not offered for new analyses (see
[`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md)).

## Supplementary File example (2022)

The impacts are entered in percent, so the adjusted impacts are in
percent:

``` r

res <- deconflate(example_supplement(), method = "published")
res$adjusted[, c("disease", "raw", "adjusted")]
#>   disease raw adjusted
#> 1      d1 2.5 2.065001
#> 2      d2 5.0 3.699294
#> 3      d3 7.5 6.748473
productivity_gap(res, 10000, "decrease", "percent")
#> $summary
#>   observed disease_free      gap  aggregate direction  effect
#> 1    10000     10215.66 215.6617 0.02111089  decrease percent
#> 
#> $attribution
#>   disease       gap  gap_main gap_interaction
#> 1      d1  21.09535  21.09535               0
#> 2      d2  56.68610  56.68610               0
#> 3      d3 137.88024 137.88024               0
```

The published gap of 10,225 units (20, 61 and 143 units per disease)
comes from rounding the adjusted impacts to 2%, 4% and 7% before
computing the gap. Without rounding, the published method gives 10,215.7
units (21.1, 56.7 and 137.9 units per disease).

## UK dairy example (2022, Tables 8-10)

[`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md)
contains the 13 diseases and the 19 odds ratios of Tables 2-3, and the
yield and fertility inputs of Tables 4-5 as two analyses. The culling
hazard ratios of Table 6 are attached as the attribute
`"hazard_ratios"`, for the hazard-ratio adapter. The paper’s culling
analysis, which converted the hazard ratios to excess annual culling
risks by treating them as odds ratios and recovered adjusted hazard
ratios by rescaling (eq. 23), exists only inside
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md).
[`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md)
holds the observed means and unit values of Table 1 for yield and
fertility, and the veterinary expenditure (the paper’s culling valuation
belongs to its historical conversion and is used only inside the
reproduction).
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
runs the whole calculation:

``` r

r22 <- reproduce_rasmussen_2022()
r22
#> <cm_reproduction> Rasmussen et al. (2022)
#> 
#> Gaps and values:
#>   analysis observed disease_free     gap  value
#>      yield     8737      9299.00 562.200 169.90
#>  fertility      401       375.10  25.910 101.80
#>    culling       27        22.55   4.454  59.48
#> 
#> Total including veterinary: 402.25
#> 
#> Comparison with Table 9:
#>                quantity  table9 package
#>      yield_disease_free 9306.00 9299.00
#>             yield_value  172.00  169.90
#>  fertility_disease_free  375.10  375.10
#>         fertility_value  101.80  101.80
#>    culling_disease_free   22.56   22.55
#>           culling_value   59.27   59.48
#>   total_with_veterinary  404.20  402.20
#> 
#> Adjusted culling hazard ratios (eq. 23) vs Table 8:
#>  disease   hr table8 package
#>       CO 1.00   1.00   1.000
#>       DA 3.83   2.68   2.671
#>      DYS 1.90   1.60   1.591
#>      FAS 1.00   1.00   1.000
#>      GIN 1.00   1.00   1.000
#>      LAM 3.40   3.00   2.999
#>      MAS 2.78   2.22   2.208
#>      MET 2.20   1.55   1.551
#>       MF 2.50   1.89   1.890
#>      NEO 1.60   1.60   1.600
#>      PTB 2.40   1.64   1.642
#>       RP 1.00   1.00   1.000
#>      SCK 2.10   1.29   1.291
r22$total
#> [1] 402.2476
```

### Comparison with Table 9

Using the inputs as printed in Tables 2-6:

|  | Published (Table 9) | Package |
|----|----|----|
| Yield potential (kg/cow/year) | 9306.32 | 9299.21 |
| Yield gap value (GBP/cow/year) | 172.05 | 169.90 |
| Calving interval potential (days) | 375.09 | 375.09 |
| Fertility gap value (GBP/cow/year) | 101.79 | 101.78 |
| Culling rate potential (%) | 22.56 | 22.55 |
| Culling gap value (GBP/cow/year) | 59.27 | 59.48 |
| Total including veterinary costs (GBP/cow/year) | 404.21 | 402.25 |

**What reproduces.**

- Fertility matches Table 9 (375.09 days, GBP 101.79 vs 101.78 after
  rounding).
- The adjusted culling hazard ratios match Table 8 to within 0.01, and
  the culling gap is within 0.5% of the published value.

**What does not.**

- **Yield.** The adjusted yield impacts cannot be reproduced exactly
  from the printed Table 4. For subclinical ketosis, note i of Table 4
  (340 kg per lactation) gives 340 / 8737 = 3.89%, not the printed
  3.05%. With 3.89%, most adjusted values match Table 8 (SCK 2.66 vs
  2.65, lameness 4.75 vs 4.76), but displaced abomasum does not (2.28 vs
  2.49). Table 10’s yield column implies further small input
  differences, for example for fasciolosis (7.04 vs 7.33).
- **Allocation of the fertility gap.** Table 10 allocates GBP 15.09 to
  cystic ovary and GBP 5.80 to metritis. Allocating the gap in
  proportion to each disease’s adjusted impact times its prevalence, as
  described in the paper (eq. 22), gives GBP 14.76 and GBP 6.04. The
  total is the same.

The original spreadsheets for the 2022 analysis are no longer available,
so the source of these differences cannot be checked. The published
corrigendum does not explain them: it corrects Table 7 only (the
dystocia loss taken from the literature, GBP 20.97 to GBP 31.92, and the
directly aggregated total, now GBP 660.26), which is not an input to the
comorbidity adjustment.

Adjusted yield impacts (percent) with the printed value of 3.05% and
with 3.89% for subclinical ketosis:

``` r

r389 <- reproduce_rasmussen_2022(yield_sck = 100 * 340 / 8737)
data.frame(disease = r22$adjusted$disease,
           sck_3.05 = r22$adjusted$yield,
           sck_3.89 = r389$adjusted$yield)
#>    disease   sck_3.05   sck_3.89
#> 1       CO 0.00000000 0.00000000
#> 2       DA 2.36694153 2.27781090
#> 3      DYS 2.93974969 2.93974969
#> 4      FAS 7.33000000 7.33000000
#> 5      GIN 3.28000000 3.28000000
#> 6      LAM 4.83062939 4.75073442
#> 7      MAS 3.76265428 3.71265733
#> 8      MET 2.40913756 2.38990351
#> 9       MF 0.09988662 0.09302367
#> 10     NEO 4.20000000 4.20000000
#> 11     PTB 4.43065737 4.43065737
#> 12      RP 6.10014912 6.08140196
#> 13     SCK 1.91867269 2.66151273
```

### What the package offers instead

**Exact solution.** `method = "simultaneous"` solves the additive
equations instead of approximating them. The example has the yield and
fertility analyses, valued with
[`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md):

``` r

uk <- example_uk_dairy_2022()
uk
#> <cm_analyses> 2 analyses on one population
#>   Diseases: 13 (CO, DA, DYS, FAS, GIN, LAM, MAS, MET, MF, NEO, PTB, RP, SCK)
#>   Disease pairs: 78 [independent (default): 59; specified: 19]
#>   - yield [% decrease]
#>   - fertility [% increase]
eco <- uk_dairy_2022_economics()
cmp <- compare_methods(uk, methods = c("published", "simultaneous"),
                       valuation = eco$valuation)
cmp$totals
#>    analysis       method raw_sum adjusted_total       gap     value
#> 1     yield    published  7.1678       6.045754 562.20724 169.89903
#> 2     yield simultaneous  7.1678       5.982142 555.91542 167.99764
#> 3 fertility    published  7.5734       6.906756  25.90677 101.77735
#> 4 fertility simultaneous  7.5734       6.448457  24.29186  95.43301
tapply(cmp$totals$value, cmp$totals$method, sum)
#>    published simultaneous 
#>     271.6764     263.4306
```

The published rows reproduce the yield and fertility values of the
reproduction above (GBP 169.90 and 101.78 per cow per year). The
simultaneous method gives GBP 168.00 for yield and 95.43 for fertility.
The published total from `reproduce_rasmussen_2022()$total` (GBP 402.25)
also includes culling, valued with the paper’s conversion, and
veterinary costs. That culling conversion is not available for new
analyses, so no simultaneous total including culling is given here.

**Hazard ratios.** For new analyses, the culling hazard ratios are
adjusted with the hazard-ratio adapter, and the culling attributable to
disease is computed without counting a cow twice. The hazard ratios of
Table 6 are attached to the example:

``` r

hr_uk <- cm_hr_model(uk$population, attr(uk, "hazard_ratios"))
compare_methods(hr_uk, overall_risk = 0.27, joint = fit_joint(uk$population))
#> <cm_comparison> methods: published, first_order, snapshot
#> Units: hazard ratio
#> 
#> Adjusted values:
#>  disease  raw published first_order snapshot
#>       CO 1.00      1.00       0.878    0.878
#>       DA 3.83      3.09       2.530    2.490
#>      DYS 1.90      1.74       1.790    1.780
#>      FAS 1.00      1.00       1.000    1.000
#>      GIN 1.00      1.00       1.000    1.000
#>      LAM 3.40      3.18       3.130    3.110
#>      MAS 2.78      2.45       2.390    2.270
#>      MET 2.20      1.80       1.760    1.690
#>       MF 2.50      2.11       2.010    1.960
#>      NEO 1.60      1.60       1.600    1.600
#>      PTB 2.40      1.88       1.620    1.620
#>       RP 1.00      1.00       0.795    0.734
#>      SCK 2.10      1.62       1.440    1.340
#> 
#> Totals:
#>       method overall_risk disease_free_risk attributable
#>    published         0.27           0.08776       0.1822
#>  first_order         0.27           0.09699       0.1730
#>     snapshot         0.27           0.10230       0.1677
```

The attributable risk times the replacement price (GBP 1335.36) gives
its value per cow and year;
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
does this with `unit_value`.

## Global dairy analysis (2024, Table 5)

Table 5 reports the means of comorbidity-adjusted impacts over 50,000
Monte Carlo draws per analysis.
[`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
encodes the input distributions. By default (`inputs = "analysis"`) it
follows the published analysis code rather than the printed tables:

- **Fixed disease probabilities.** Prevalence and incidence were not
  drawn in the adjustment. Each probability is `1 - exp(-incidence)` at
  the unrounded global mean, except paratuberculosis (a prevalence) and
  subclinical mastitis, which was entered without conversion (0.409
  rather than 0.336).
- **PERT distributions** use the central value in Tables 2-4 as the
  mode.
- **Unrounded impact parameters**, e.g. a clinical ketosis yield impact
  of 0.4322% rather than 0.43%.
- **Culling** was adjusted on the hazard ratio scale: HR - 1 was
  adjusted, and 1 was added back. Normal standard deviations were scaled
  by (HR - 1) / HR, and metritis used a PERT distribution (mode 1.12)
  rather than the normal distribution in Table 4 (mean 1.05). This
  analysis (`culling_hr_minus_1`) exists only inside
  [`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md),
  for reproduction;
  [`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)
  and
  [`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
  have the yield and fertility analyses only.

`sampler_global_dairy(inputs = "tables")` uses Tables 2-4 as printed
instead.
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
runs the published method on the sampler, with the historical culling
analysis added, and sets the means beside Table 5 (yield and fertility
in percent, culling as hazard ratios). Run it with many more draws for
stable estimates; 200 are used here to keep the vignette fast.

``` r

r24 <- reproduce_rasmussen_2024(n_draws = 200)
cmp24 <- r24$comparison
cmp24[cmp24$analysis == "yield", c("disease", "table5", "mean", "median", "stability")]
#>    disease table5       mean     median stability
#> 1       CK   0.03 0.03442976 0.03156076        ok
#> 2       CM   1.36 1.33089121 1.37155476        ok
#> 3       DA   1.18 1.16913885 0.98908208 imprecise
#> 4      DYS   3.48 3.46250651 3.45262679        ok
#> 5      LAM   2.62 2.56465441 2.52503159        ok
#> 6      MET   2.87 2.91014635 2.87194129        ok
#> 7       MF   0.07 0.07037444 0.06873653        ok
#> 8       OC   2.59 2.51844622 2.56329222        ok
#> 9      PTB   3.37 3.32022926 3.29087156        ok
#> 10      RP   2.30 2.23600434 2.22201614        ok
#> 11     SCK   7.11 7.14187599 7.06581819        ok
#> 12     SCM   5.58 5.69548041 5.54342870        ok
cmp24[cmp24$analysis == "culling_hr", c("disease", "table5", "mean", "median", "stability")]
#>    disease table5     mean   median     stability
#> 25      CK   1.18 1.184261 1.182960            ok
#> 26      CM   1.90 1.894408 1.877761            ok
#> 27      DA   2.75 2.774447 2.644105            ok
#> 28     DYS   1.18 1.193927 1.125758 possible_pole
#> 29     LAM   1.40 1.385253 1.379647            ok
#> 30     MET   1.03 1.027202 1.016015     imprecise
#> 31      MF   2.64 2.673999 2.654076            ok
#> 32      OC   1.51 1.478334 1.471560            ok
#> 33     PTB   2.07 2.058094 2.050252            ok
#> 34      RP   1.29 1.291959 1.282293            ok
#> 35     SCK   1.67 1.667363 1.669365            ok
#> 36     SCM   1.25 1.252106 1.251805            ok
```

With 200,000 draws (`inst/validation/reference_2024_analysis.py`), the
means are:

| Disease | Yield, Table 5 | Yield, package | Culling HR, Table 5 | Culling HR, package |
|----|----|----|----|----|
| CK | 0.03 | 0.04 | 1.18 | 1.18 |
| CM | 1.36 | 1.36 | 1.90 | 1.90 |
| DA | 1.18 | 1.18 | 2.75 | 2.74 |
| DYS | 3.48 | 3.49 | 1.18 | 1.20 |
| LAM | 2.62 | 2.60 | 1.40 | 1.39 |
| MET | 2.87 | 2.86 | 1.03 | 1.03 |
| MF | 0.07 | 0.07 | 2.64 | 2.64 |
| OC | 2.59 | 2.51 | 1.51 | 1.48 |
| PTB | 3.37 | 3.32 | 2.07 | 2.06 |
| RP | 2.30 | 2.30 | 1.29 | 1.29 |
| SCK | 7.11 | 7.11 | 1.67 | 1.67 |
| SCM | 5.58 | 5.58 | 1.25 | 1.25 |

**Negative odds-ratio draws.** About 3% of draws of the lameness:ovarian
cyst odds ratio (normal, mean 2.63, SD 1.44) are negative, as are about
1% of lameness:paratuberculosis draws. The analysis code used these
draws as they were, taking the absolute value of a negative discriminant
in eq. 11. A negative odds ratio has no valid 2x2 table, so the package
truncates normal odds-ratio distributions at zero. This explains the
paratuberculosis yield difference (3.38 without truncation). It does not
explain ovarian cyst (2.51 with or without truncation): in the analysis
output, the spread of the adjusted ovarian cyst yield impact is several
times larger than for the other diseases, which suggests that a few
extreme draws drive its mean.

**Fertility.** Fertility means for clinical mastitis, dystocia,
lameness, milk fever, ovarian cyst and retained placenta agree with
Table 5 to within about 0.05 (for example clinical mastitis 6.10 vs
6.09). Clinical ketosis (0.30 vs 0.34) and paratuberculosis (4.08 vs
4.23) differ more. Where the raw impact distribution extends below zero
(displaced abomasum, metritis, subclinical ketosis and subclinical
mastitis), the published approximation `m^2 / (m + c)` is unstable:
draws in which `m + c` is close to zero produce extreme values. For
subclinical ketosis the Monte Carlo mean does not settle even with
200,000 draws, and for subclinical mastitis the archived inputs give
0.52% rather than the 0.04% in Table 5. Medians, or the exact
(simultaneous) solution, are more robust summaries. The `stability`
column flags estimates whose means are unreliable (see
[`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md)):

``` r

cmp24[cmp24$analysis == "fertility", c("disease", "table5", "mean", "median", "stability")]
#>    disease table5       mean     median     stability
#> 13      CK   0.34  0.3111188  0.2747364            ok
#> 14      CM   6.09  6.1548342  6.1297197            ok
#> 15      DA   0.78  1.0023297  0.3848882    heavy_tail
#> 16     DYS   1.11  1.0743470  0.9363284            ok
#> 17     LAM   1.86  1.9423520  1.7641199            ok
#> 18     MET  11.22 10.5990086 10.4920940 possible_pole
#> 19      MF   1.06  1.0573668  1.0456397            ok
#> 20      OC   9.03  9.0885437  8.4399206            ok
#> 21     PTB   4.23  3.7642570  3.7627776 possible_pole
#> 22      RP   3.74  3.6555091  3.5283150            ok
#> 23     SCK   0.39  0.2131768  0.5177367 possible_pole
#> 24     SCM   0.04  0.4843150  0.2616170     imprecise
```

The simultaneous solution at the central values shows the inconsistency
directly. Several adjusted impacts change sign, while the aggregates
move less (about 3% for yield and 18% for fertility):

``` r

res24 <- deconflate(example_global_dairy(), warn = FALSE)
vapply(res24, function(r) r$diagnostics$sign_changes, character(1))
#>                   yield               fertility 
#>        "CK, CM, DA, MF" "CK, DA, LAM, SCK, SCM"
pub24 <- deconflate(example_global_dairy(), method = "published")
rbind(published = vapply(pub24, function(r) r$totals$adjusted_total, numeric(1)),
      simultaneous = vapply(res24, function(r) r$totals$adjusted_total, numeric(1)))
#>                 yield fertility
#> published    7.280553  4.794033
#> simultaneous 7.494192  3.933738
```

**Culling losses.** The adjusted culling impacts of the 2024 analysis
are adjusted HR - 1; adding 1 gives the adjusted hazard ratios. The
central values are in `r24$central$culling_hr_minus_1`, and they equal
the `"published"` method of
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md):

``` r

cu <- r24$central$culling_hr_minus_1$adjusted
hr_pub <- deconflate_hr(example_global_dairy_hr(), method = "published")
data.frame(disease = cu$disease, hr = 1 + cu$raw, adjusted_hr = 1 + cu$adjusted,
           deconflate_hr = hr_pub$adjusted$adjusted)
#>    disease       hr adjusted_hr deconflate_hr
#> 1       CK 1.500100    1.177580      1.177580
#> 2       CM 2.300000    1.903941      1.903941
#> 3       DA 2.851179    2.197930      2.197930
#> 4      DYS 1.258143    1.098377      1.098377
#> 5      LAM 1.744976    1.380683      1.380683
#> 6      MET 1.116444    1.012411      1.012411
#> 7       MF 2.999886    2.647637      2.647637
#> 8       OC 1.620000    1.458644      1.458644
#> 9      PTB 2.310508    2.047235      2.047235
#> 10      RP 1.599928    1.284496      1.284496
#> 11     SCK 1.920000    1.675253      1.675253
#> 12     SCM 1.449996    1.254928      1.254928
```

For the losses, the paper converted each adjusted hazard ratio to an
excess culling risk relative to the overall culling rate `r`,
`HR * r / (HR * r + 1 - r) - r`, which treats the hazard ratio as an
odds ratio and counts cows with several diseases more than once. For new
analyses use
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
(see
[`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md)).
