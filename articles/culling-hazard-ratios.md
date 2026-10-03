# Culling and hazard ratios

``` r

library(deconflate)
```

Culling (and mortality) impacts are usually reported as hazard ratios
(HRs) from survival models. A raw HR for one disease is conflated with
the HRs of associated diseases, as a yield impact is, but HRs combine
multiplicatively and are not additive impacts. They are therefore
adjusted by a separate adapter, outside the additive engine of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md):

- [`cm_hazard_ratios()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hazard_ratios.md)
  describes the raw HRs (one per disease; use 1 for no effect) and
  states their estimand: `"snapshot_crude"` or `"snapshot_stratified"`;
- [`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md)
  combines them with a population;
- [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
  adjusts them;
- [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  turns adjusted HRs into the part of the overall culling risk that
  disease causes, and allocates it to diseases.

This is a separate model with its own assumptions (the snapshot
hazard-multiplier model, below), not a general conversion of published
hazard ratios. Its estimands are defined by the model:

- `"snapshot_crude"`: the ratio of the average hazard among animals with
  the disease to that among animals without it, at the start of
  follow-up, in the population described by the probabilities and
  associations;
- `"snapshot_stratified"`: the same ratio within strata of the diseases
  in `adjusted_for`, combined across strata.

A published Cox HR is neither of these. Entering one as a snapshot ratio
is an approximation that you state, not a conversion the package makes:
it is reasonable when follow-up is short relative to the hazards, and it
should be reported and checked (see “What the snapshot model is not”
below and in
[`?deconflate_hr`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)).
For this reason the estimand has no default. The names `"crude"` and
`"adjusted"` used before version 0.3.0 are not accepted:

``` r

tryCatch(cm_hazard_ratios("d1", 1.5, estimand = "crude"),
         deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "Hazard-ratio estimands are now named \"snapshot_crude\" and \"snapshot_stratified\" (deconflate 0.3.0); see ?cm_hazard_ratios for what they mean."
```

## Methods

``` r

pop <- example_supplement()   # its yield impacts are not used here
hr <- cm_hr_model(pop, cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), estimand = "snapshot_crude"))
hr
#> <cm_hr_model>
#> <cm_population>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [specified: 3]
#>   Hazard ratios: 3 (0 snapshot_stratified)
compare_methods(hr)
#> <cm_comparison> methods: published, first_order, snapshot
#> Units: hazard ratio
#> 
#> Adjusted values:
#>  disease raw published first_order snapshot
#>       d1 1.5      1.41        1.40     1.38
#>       d2 2.0      1.91        1.89     1.89
#>       d3 1.3      1.19        1.17     1.14
```

- `"snapshot"` (the default) is a snapshot hazard-multiplier model. A
  cow’s hazard is `h0 * exp(sum_i beta_i * D_i)`, and the raw HR of
  disease `i` is taken to be the ratio of the average hazard multiplier
  among cows with and without `i`, over the fitted distribution of
  disease combinations
  ([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md))
  at the start of follow-up. The `beta`s are solved so that these ratios
  equal the raw HRs, and `exp(beta)` are the adjusted HRs.
- `"first_order"` is the log-linear approximation
  `log(HR_raw) = A beta`, with the conflation matrix `A` of
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).
  It needs only the pairwise tables.
- `"published"` is the approach of Rasmussen et al. (2024): HR - 1
  adjusted with eq. 16 and 1 added back. It is kept for reproduction and
  comparison.

``` r

deconflate_hr(hr)
#> <cm_hr_result> method: snapshot
#> 
#>  disease raw adjusted   change       estimand
#>       d1 1.5    1.383 -0.07790 snapshot_crude
#>       d2 2.0    1.891 -0.05453 snapshot_crude
#>       d3 1.3    1.144 -0.12001 snapshot_crude
#> 
#> Diagnostics:
#>  max_reconstruction_residual n_sign_changes condition_number
#>                     2.22e-16              0             1.68
#>                                       feasibility
#>  joint distribution fitted (max residual 2.5e-11)
```

## Source hazard ratios adjusted for other diseases

If a source HR comes from a model that already included other diseases,
the closest snapshot estimand is `"snapshot_stratified"`, with those
diseases in `adjusted_for` (ids separated by `";"`), or `"all"` for
every other disease. In the snapshot model, the ratio for such an
estimate is computed within strata of its adjustment set and combined
across strata with Mantel-Haenszel-type weights. With
`adjusted_for = "all"`, the HR is the disease’s own hazard multiplier
and is used as it is. A Cox coefficient adjusted for the same diseases
matches this only under the snapshot model’s assumptions:

``` r

hr_adj <- cm_hr_model(pop, cm_hazard_ratios(
  c("d1", "d2", "d3"), c(1.5, 2.0, 1.3),
  estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_stratified"),
  adjusted_for = c("d2", NA, "all")
))
deconflate_hr(hr_adj)$adjusted
#>   disease raw adjusted       change            estimand adjusted_for
#> 1      d1 1.5 1.509402  0.006268033 snapshot_stratified           d2
#> 2      d2 2.0 1.821704 -0.089148021      snapshot_crude         <NA>
#> 3      d3 1.3 1.300000  0.000000000 snapshot_stratified          all
```

The published approach is defined for `"snapshot_crude"` HRs only:

``` r

compare_methods(hr_adj)$failed
#>                                                                published 
#> "The published (2024) approach is defined for crude hazard ratios only."
```

## What the snapshot model is not

A Cox HR estimated over follow-up is not, in general, the snapshot
ratio: cows with high hazards leave first, so the mixture of disease
combinations among survivors changes over time, and the marginal HR
changes with it. The snapshot model is exact for its own estimand
(instantaneous marginal ratios at baseline). For published Cox
coefficients it is an approximation, which is reasonable when follow-up
is short relative to the hazards or the diseases are rare.

## Culling attributable to disease

[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
finds the baseline hazard for which the population culling risk,
averaged over the distribution of disease combinations, equals the
observed rate, and compares it with the risk of a disease-free cow. A
cow’s risk cannot exceed 1, so a cow with several diseases is counted
once. The attributable risk is allocated to diseases by Shapley values
over all disease combinations
([`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md));
the allocation adds up to the attributable risk.

``` r

attributable_risk(deconflate_hr(hr), overall_risk = 0.25, unit_value = 1300)
#> <cm_attributable> snapshot hazard-multiplier model
#>   Overall risk 0.25; disease-free risk 0.2134; attributable 0.03665 (14.7% of the overall risk)
#>   Value: 47.64
#> 
#>  disease hr_adjusted attributable  share  value
#>       d1       1.383     0.007371 0.2011  9.582
#>       d2       1.891     0.023488 0.6409 30.534
#>       d3       1.144     0.005791 0.1580  7.528
```

`unit_value` is the value of a cow removed (e.g. replacement price less
salvage value).

## The global dairy example

[`example_global_dairy_hr()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy_hr.md)
contains the culling HRs of Rasmussen et al. (2024), on the population
of
[`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md).
They are published Cox HRs, entered as `"snapshot_crude"`; that is an
assumption, and the results below hold under it.
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
with `overall_risk` adds the attributable risk per method (without the
allocation). The global average replacement rate is 23.66% (Rasmussen et
al. 2024, Table 1). The joint distribution of the 12 diseases is fitted
once and passed on with `joint`, so that it is not refitted for each
method:

``` r

hr_gd <- example_global_dairy_hr()
j_gd <- fit_joint(hr_gd$population)
compare_methods(hr_gd, overall_risk = 0.2366, joint = j_gd)
#> <cm_comparison> methods: published, first_order, snapshot
#> Units: hazard ratio
#> 
#> Adjusted values:
#>  disease  raw published first_order snapshot
#>       CK 1.50      1.18       0.962    0.982
#>       CM 2.30      1.90       1.840    1.820
#>       DA 2.85      2.20       1.830    1.850
#>      DYS 1.26      1.10       1.120    1.050
#>      LAM 1.74      1.38       1.290    1.220
#>      MET 1.12      1.01       0.737    0.720
#>       MF 3.00      2.65       2.620    2.660
#>       OC 1.62      1.46       1.550    1.560
#>      PTB 2.31      2.05       2.020    2.020
#>       RP 1.60      1.28       1.230    1.180
#>      SCK 1.92      1.68       1.740    1.710
#>      SCM 1.45      1.25       1.250    1.210
#> 
#> Totals:
#>       method overall_risk disease_free_risk attributable
#>    published       0.2366            0.1076       0.1290
#>  first_order       0.2366            0.1144       0.1222
#>     snapshot       0.2366            0.1200       0.1166
```

The published column reproduces the paper’s adjustment at the central
values (Table 5 reports Monte Carlo means). The methods agree within
about 0.1 for most diseases, and differ most where comorbidity is
heaviest (clinical ketosis, lameness, displaced abomasum, metritis). For
clinical ketosis and metritis, the snapshot and first-order models give
adjusted HRs below 1: under the multiplicative model, the raw HRs are
smaller than the associated diseases alone would produce.
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
warns about this:

``` r

glob <- deconflate_hr(hr_gd, joint = j_gd)
#> Warning: Adjusted hazard ratios cross 1 for CK, MET. The raw hazard ratios are
#> smaller than the associated diseases alone would produce; check the estimands
#> and source populations.
ar <- attributable_risk(glob, overall_risk = 0.2366, unit_value = 1299.33 - 785.86)
ar
#> <cm_attributable> snapshot hazard-multiplier model
#>   Overall risk 0.2366; disease-free risk 0.12; attributable 0.1166 (49.3% of the overall risk)
#>   Value: 59.86
#> 
#>  disease hr_adjusted attributable      share    value
#>       CK      0.9821   -9.666e-05 -0.0008292 -0.04963
#>       CM      1.8194    2.981e-02  0.2557145 15.30613
#>       DA      1.8477    2.843e-03  0.0243863  1.45967
#>      DYS      1.0505    4.715e-04  0.0040444  0.24209
#>      LAM      1.2230    7.989e-03  0.0685346  4.10223
#>      MET      0.7198   -4.730e-03 -0.0405736 -2.42859
#>       MF      2.6601    5.124e-03  0.0439548  2.63097
#>       OC      1.5580    8.424e-03  0.0722672  4.32565
#>      PTB      2.0165    1.455e-02  0.1247870  7.46929
#>       RP      1.1817    3.403e-03  0.0291948  1.74749
#>      SCK      1.7107    3.590e-02  0.3079742 18.43420
#>      SCM      1.2109    1.289e-02  0.1105452  6.61683
```

(`unit_value` is the replacement price less the culled-cow price,
Rasmussen et al. 2024, Table 1.)

The model attributes about 11.7 percentage points of the 23.7% culling
rate to disease. The 2024 paper instead converted each adjusted HR to an
excess risk relative to the overall rate,
`HR * r / (HR * r + 1 - r) - r`, and summed these excess risks weighted
by prevalence. With the published adjusted HRs, that sum is about 14.8
percentage points, because it counts cows with several diseases more
than once.

## Culling alongside the additive analyses

Hazard ratios are kept out of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
and the contribution tables. Values from the two routes can be put side
by side. Here the yield analysis follows the paper (published method;
global average yield of 5013 kg and milk price of 0.5981 per kg,
Rasmussen et al. 2024, Table 1):

``` r

yres <- deconflate(example_global_dairy()$models$yield, method = "published")
yv <- value_losses(productivity_gap(yres, 5013, "decrease", "percent"), unit_value = 0.5981)
c(yield = yv$value, culling = ar$summary$value)
#>     yield   culling 
#> 235.43177  59.85632
```

## Methods that are not available for new analyses

The published analyses converted hazard ratios in ways that the package
does not offer for new analyses:

- Rasmussen et al. (2022) treated each HR as an odds ratio of a 2x2
  table of disease by culling, to obtain an excess annual culling risk,
  and recovered adjusted HRs by rescaling (eq. 23);
- Rasmussen et al. (2024) adjusted HR - 1 as if it were an additive
  impact, and valued culling with the excess risk relative to the
  overall rate given above.

These conversions are used only inside
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
and
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md),
so that the published tables can be recomputed (see
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).

## Caveats

- The models assume a constant hazard within the period and that the HRs
  apply for the whole period. Most source HRs come from models in which
  disease status changes during the lactation.
- Hazard ratios are not collapsible: even without confounding, an
  average HR differs from the HR within subgroups. The snapshot model
  accounts for this given the distribution of disease combinations; the
  first-order and published methods do not.
- The snapshot model and
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  use the whole distribution of disease combinations, so they depend on
  the maximum-entropy assumption and on any three-way terms
  ([`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)).
- Culling and death compete. Make sure the source HRs are of the same
  kind (cause-specific or subdistribution).
