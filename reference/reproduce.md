# Reproduce the published analyses

Recompute the tables of Rasmussen et al. (2022) and Rasmussen et al.
(2024) with the published (eq. 16) method and the conversions used in
the papers, and set them beside the printed values.

## Usage

``` r
reproduce_rasmussen_2022(yield_sck = 3.05)

reproduce_rasmussen_2024(
  n_draws = 5000,
  seed = 2024,
  inputs = c("analysis", "tables")
)
```

## Arguments

- yield_sck:

  Raw yield impact of subclinical ketosis in percent (see
  [`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md)).

- n_draws:

  Number of Monte Carlo draws.

- seed:

  Random seed.

- inputs:

  `"analysis"` or `"tables"` (see
  [`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)).

## Value

A `cm_reproduction` object: a list with `adjusted` (adjusted values per
disease), `comparison` (package beside the published values), and for
2022 `gaps`, `values` and `total`, or for 2024 the Monte Carlo run `mc`.

## Details

These functions exist to document and check the published numbers. Two
of the conversions they use are not offered for new analyses:

- 2022: culling hazard ratios were converted to excess annual culling
  risks by treating them as odds ratios, and adjusted hazard ratios were
  recovered by rescaling (eq. 23);

- 2024: hazard ratios minus 1 were adjusted as additive impacts.

For new analyses of hazard ratios use
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
and
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).

## 2022 (Tables 8-10)

`reproduce_rasmussen_2022()` adjusts yield and calving interval for the
UK example
([`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md))
and the historical culling analysis, computes the productivity gaps and
their values with the paper's economic inputs (Table 1, and its culling
valuation), adds veterinary expenditure, and back-converts the culling
impacts to hazard ratios. (The gap and valuation steps exist only here:
the package itself does not value impacts.) Fertility and the culling
hazard ratios reproduce the paper; yield does not reproduce exactly from
the printed Table 4 (see
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).

## 2024 (Table 5)

`reproduce_rasmussen_2024()` runs
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
with the samplers of
[`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
plus the historical culling analysis (HR - 1, named
`culling_hr_minus_1`) and the published method, and reports the means
beside Table 5. Culling is reported on the hazard-ratio scale (1 +
adjusted HR - 1). Table 5 used 50,000 draws; use at least several
thousand for stable means. Some fertility means are unstable for the
published method (see the `stability` column).

## Examples

``` r
r22 <- reproduce_rasmussen_2022()
r22$comparison
#>                 quantity  table9    package
#> 1     yield_disease_free 9306.32 9299.20724
#> 2            yield_value  172.05  169.89903
#> 3 fertility_disease_free  375.09  375.09323
#> 4        fertility_value  101.79  101.77735
#> 5   culling_disease_free   22.56   22.54568
#> 6          culling_value   59.27   59.48117
#> 7  total_with_veterinary  404.21  402.24755
# \donttest{
r24 <- reproduce_rasmussen_2024(n_draws = 500)
r24$comparison[r24$comparison$analysis == "yield", ]
#>    analysis disease table5       mean     median         mcse stability
#> 1     yield      CK   0.03 0.03462749 0.03072734 0.0008424078        ok
#> 2     yield      CM   1.36 1.36906335 1.36655094 0.0223938872        ok
#> 3     yield      DA   1.18 1.17239741 0.95209126 0.0473520348        ok
#> 4     yield     DYS   3.48 3.45118251 3.44400596 0.0446067755        ok
#> 5     yield     LAM   2.62 2.57534992 2.54697374 0.0324972483        ok
#> 6     yield     MET   2.87 2.91654349 2.88149871 0.0449770720        ok
#> 7     yield      MF   0.07 0.06997760 0.06876974 0.0003451446        ok
#> 8     yield      OC   2.59 2.51416504 2.50470235 0.0246006448        ok
#> 9     yield     PTB   3.37 3.29105837 3.27534508 0.0325060426        ok
#> 10    yield      RP   2.30 2.34976783 2.32481004 0.0410583147        ok
#> 11    yield     SCK   7.11 7.15404733 7.13358710 0.0510119767        ok
#> 12    yield     SCM   5.58 5.61239650 5.54026950 0.0534154829        ok
# }
```
