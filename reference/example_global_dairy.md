# Global dairy inputs from Rasmussen et al. (2024), at their central values

The 12 diseases of Rasmussen et al. (2024): global (herd-size weighted)
lactational incidence (prevalence for PTB), pooled inter-disease odds
ratios (Table 3, pairs not listed are independent), and raw impacts on
yield (% decrease), fertility (% increase in calving interval) and,
optionally, culling (hazard ratio minus 1). Central values are means of
normal distributions and modes of PERT distributions.

## Usage

``` r
example_global_dairy(
  inputs = c("analysis", "tables"),
  culling = TRUE,
  culling_scale = c("excess_hr", "hazard_ratio")
)
```

## Arguments

- inputs:

  `"analysis"` or `"tables"`; see the section above.

- culling:

  Include the culling outcome?

- culling_scale:

  `"excess_hr"` (HR - 1, as in the paper) or `"hazard_ratio"`; see the
  section on culling.

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
with outcomes `"yield"`, `"fertility"` and, if `culling = TRUE`,
`"culling"`.

## Two versions of the inputs

- `inputs = "analysis"` (default) uses the inputs of the published
  analysis code (1st revision), which reproduce Table 5:

  - disease probabilities are fixed at `1 - exp(-incidence)`, using the
    unrounded global mean incidence. Subclinical mastitis (SCM) was
    entered without conversion (0.4094 rather than 0.3360); this is kept
    here (`type = "probability"`);

  - impacts are the unrounded central values (e.g. clinical ketosis
    yield 0.4322% rather than 0.43%);

  - culling is entered as HR - 1, as in the analysis. Metritis uses the
    analysis value (PERT mode 1.116) rather than Table 4 (normal, mean
    1.05).

- `inputs = "tables"` uses Tables 2-4 as printed, with incidence
  converted for every disease except PTB.

## Culling

With `culling_scale = "excess_hr"` (default, as in the analysis),
culling impacts are hazard ratios minus 1 on the `"absolute"` scale, so
that
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
adjusts the excess hazard ratio. Convert results back with
`adjusted_hr(res, method = "excess_hr")`. These impacts are not excess
risks, so do not pass them to
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md);
the paper converts adjusted hazard ratios to excess culling risk with
`hr_to_risk(..., method = "overall_odds")`.

With `culling_scale = "hazard_ratio"`, culling impacts are hazard ratios
(see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)).
`method = "published"` then gives the same adjusted hazard ratios as the
paper, while `"simultaneous"` and `"global"` use the multiplicative
model;
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
gives the culling attributable to disease.

The paper's Table 5 reports means of adjusted impacts over Monte Carlo
draws
([`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)),
which differ from adjusting the central values (e.g. displaced abomasum
yield: 1.18 vs 0.79).

## Examples

``` r
res <- deconflate(example_global_dairy(), method = "published")
adjusted_hr(res, method = "excess_hr")
#>    disease       hr   excess excess_adjusted hr_adjusted
#> 1       CK 1.500100 0.500100      0.17758038    1.177580
#> 2       CM 2.300000 1.300000      0.90394063    1.903941
#> 3       DA 2.851179 1.851179      1.19792993    2.197930
#> 4      DYS 1.258143 0.258143      0.09837667    1.098377
#> 5      LAM 1.744976 0.744976      0.38068270    1.380683
#> 6      MET 1.116444 0.116444      0.01241058    1.012411
#> 7       MF 2.999886 1.999886      1.64763745    2.647637
#> 8       OC 1.620000 0.620000      0.45864412    1.458644
#> 9      PTB 2.310508 1.310508      1.04723492    2.047235
#> 10      RP 1.599928 0.599928      0.28449635    1.284496
#> 11     SCK 1.920000 0.920000      0.67525327    1.675253
#> 12     SCM 1.449996 0.449996      0.25492756    1.254928

# Culling hazard ratios under the multiplicative model
m <- example_global_dairy(culling_scale = "hazard_ratio")
compare_methods(m, methods = c("published", "simultaneous"))
#> <cm_comparison> methods: published, simultaneous
#> 
#> Adjusted impacts:
#>    outcome disease unit    raw published simultaneous
#>      yield      CK    %  0.432    0.0247      -5.3700
#>      yield      CM    %  3.250    1.3300      -0.3380
#>      yield      DA    %  2.840    0.7940      -2.6300
#>      yield     DYS    %  4.920    3.5700       4.3100
#>      yield     LAM    %  4.810    2.5300       1.9800
#>      yield     MET    %  5.610    2.8400       2.7900
#>      yield      MF    %  0.537    0.0690      -1.5400
#>      yield      OC    %  3.750    2.6400       3.2200
#>      yield     PTB    %  4.300    3.2300       3.9400
#>      yield      RP    %  4.200    2.2600       2.4700
#>      yield     SCK    %  8.400    7.1000       8.2800
#>      yield     SCM    %  6.290    5.5900       6.5800
#>  fertility      CK    %  1.450    0.3260      -1.5600
#>  fertility      CM    %  8.420    6.2000       7.5600
#>  fertility      DA    %  1.080    0.1530      -3.3300
#>  fertility     DYS    %  2.400    1.0600       0.9910
#>  fertility     LAM    %  3.300    1.1200      -2.0400
#>  fertility     MET    % 14.700   10.8000      13.0000
#>  fertility      MF    %  2.410    1.0600       2.0000
#>  fertility      OC    %  9.690    7.8500       9.0400
#>  fertility     PTB    %  5.350    3.9900       4.7400
#>  fertility      RP    %  6.760    3.7100       2.5500
#>  fertility     SCK    %  1.120    0.3490      -0.0399
#>  fertility     SCM    %  0.264    0.0322      -1.2500
#>    culling      CK   HR  1.500    1.1800       0.9620
#>    culling      CM   HR  2.300    1.9000       1.8400
#>    culling      DA   HR  2.850    2.2000       1.8300
#>    culling     DYS   HR  1.260    1.1000       1.1200
#>    culling     LAM   HR  1.740    1.3800       1.2900
#>    culling     MET   HR  1.120    1.0100       0.7370
#>    culling      MF   HR  3.000    2.6500       2.6200
#>    culling      OC   HR  1.620    1.4600       1.5500
#>    culling     PTB   HR  2.310    2.0500       2.0200
#>    culling      RP   HR  1.600    1.2800       1.2300
#>    culling     SCK   HR  1.920    1.6800       1.7400
#>    culling     SCM   HR  1.450    1.2500       1.2500
#> 
#> Sign changes (adjusted impact on the other side of zero, or of 1 for hazard ratios):
#>   simultaneous: yield CK, yield CM, yield DA, yield MF, fertility CK, fertility DA, fertility LAM, fertility SCK, fertility SCM, culling CK, culling MET
#> 
#> Totals:
#>    outcome       method raw_loss adjusted_loss
#>  fertility    published  0.07472       0.04794
#>      yield    published  0.09931       0.07281
#>  fertility simultaneous  0.07472       0.03934
#>      yield simultaneous  0.09931       0.07494
```
