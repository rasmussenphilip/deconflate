# Using your own data

``` r

library(deconflate)
```

Inputs can be kept in CSV files (one per table) or typed in R as data
frames.
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
checks every table, reports all problems at once, and returns a model
(and a Monte Carlo sampler if you give an uncertainty table).

## Start from the template

[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set of files to a folder:

``` r

dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 5 files to /tmp/RtmpcwW4DG/my-inputs
list.files(dir)
#> [1] "associations.csv" "diseases.csv"     "impacts.csv"      "interactions.csv"
#> [5] "uncertainty.csv"
```

The five files are:

| File | Contents | Required columns |
|----|----|----|
| `diseases.csv` | one row per disease | `id`, `value` |
| `associations.csv` | one row per associated pair | `disease1`, `disease2` |
| `impacts.csv` | one row per disease and outcome | `disease`, `outcome`, `value` |
| `interactions.csv` | optional impact interactions | `disease1`, `disease2`, `value`, `outcome` |
| `uncertainty.csv` | optional input distributions | `key`, `dist` |

Only `diseases.csv` is required. The optional columns (e.g. `type`,
`measure`, `scale`, `direction`, `units`, `adjusted_for`, `source`) are
described in
[`?cm_read_inputs`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md).
Here is the template’s impacts table:

``` r

read.csv(file.path(dir, "impacts.csv"))
#>   disease outcome value        scale direction units adjusted_for
#> 1     LAM   yield  4.81      percent  decrease    NA           NA
#> 2     SCK   yield  8.40      percent  decrease    NA           NA
#> 3     MET   yield  5.61      percent  decrease    NA           NA
#> 4     LAM culling  1.74 hazard_ratio              NA           NA
#> 5     SCK culling  1.92 hazard_ratio              NA           NA
#> 6     MET culling  1.50 hazard_ratio              NA           NA
#>                source
#> 1 Illustrative values
#> 2 Illustrative values
#> 3 Illustrative values
#> 4 Illustrative values
#> 5 Illustrative values
#> 6 Illustrative values
```

Yield impacts are percent decreases; culling impacts are hazard ratios
(`scale = hazard_ratio`), adjusted on the log scale (see
[`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md)).
Every disease needs a row for every outcome; use 0 for no impact.

Two more example sets are installed with the package:

- `global_dairy_2024`: the 12 diseases, 38 odds ratios, yield, fertility
  and culling impacts, and input distributions of Rasmussen et
  al. (2024), as used in the published analysis (culling as hazard
  ratios);
- `example_with_errors`: a small set with deliberate mistakes, to see
  how problems are reported.

``` r

system.file("extdata", package = "deconflate")
#> [1] "/home/runner/work/_temp/Library/deconflate/extdata"
```

## Read the files

``` r

inp <- cm_read_inputs(dir = dir)
inp
#> <cm_inputs>
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
#>   Outcomes: yield, culling
#>   Uncertain inputs: 5 (use $sampler with cm_monte_carlo())
```

`inp$model` is an ordinary model:

``` r

res <- deconflate(inp$model, method = "global")
res$adjusted[, c("outcome", "disease", "raw", "adjusted")]
#>   outcome disease    raw   adjusted
#> 1   yield     LAM 0.0481 0.02870596
#> 2   yield     SCK 0.0840 0.07819114
#> 3   yield     MET 0.0561 0.03169881
#> 4 culling     LAM 1.7400 1.53155018
#> 5 culling     SCK 1.9200 1.78606195
#> 6 culling     MET 1.5000 1.14080183
```

Compare the methods side by side:

``` r

compare_methods(inp$model)
#> <cm_comparison> methods: published, simultaneous, global
#> 
#> Adjusted impacts:
#>  outcome disease unit  raw published simultaneous global
#>    yield     LAM    % 4.81      3.16         2.87   2.87
#>    yield     SCK    % 8.40      7.51         7.82   7.82
#>    yield     MET    % 5.61      3.52         3.17   3.17
#>  culling     LAM   HR 1.74      1.55         1.53   1.53
#>  culling     SCK   HR 1.92      1.81         1.80   1.79
#>  culling     MET   HR 1.50      1.26         1.15   1.14
#> 
#> Totals:
#>  outcome       method raw_loss adjusted_loss
#>    yield    published  0.04966       0.04006
#>    yield simultaneous  0.04966       0.04015
#>    yield       global  0.04966       0.04015
```

## When something is wrong

All problems are listed with their table, row (counting from the first
row below the header) and column, so they can be fixed in one pass.
[`cm_check_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_check_inputs.md)
returns them without stopping:

``` r

cm_check_inputs(
  diseases = data.frame(id = c("LAM", "SCK", "SCK"), value = c(0.25, 1.4, 0.3)),
  associations = data.frame(disease1 = "LAM", disease2 = "MET", value = 2.0),
  impacts = data.frame(disease = c("LAM", "SCK"), outcome = "yield",
                       value = c("4.81", "8,40"), scale = "percent")
)
#> Found 4 problem(s) in the inputs:
#>   diseases, row 3, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 2, column 'value': Gives a probability of 1.4; it must be strictly between 0 and 1.
#>   associations, row 1, column 'disease2': Unknown disease 'MET' (not in the diseases table).
#>   impacts, row 2, column 'value': '8,40' is not a number.
```

[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
stops with the same list when there are errors. The installed example
with errors shows more cases:

``` r

cm_check_inputs(dir = system.file("extdata", "example_with_errors", package = "deconflate"))
#> Found 10 problem(s) in the inputs:
#>   diseases, row 4, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   associations, row 3, column 'disease2': Unknown disease 'CM' (not in the diseases table).
#>   associations, row 4: Duplicate pair SCK:LAM.
#>   associations, row 5, column 'value': Missing value for measure OR.
#>   impacts, row 2, column 'value': '8,40' is not a number.
#>   impacts, row 5, column 'value': Hazard ratios must be positive.
#>   impacts: Outcome 'yield' has no impact for: MET. Add a row with value 0 for no impact.
#>   uncertainty, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   uncertainty, row 3, column 'key': Key 'impact:yield:CM' does not match an input (expected prob:<disease>, assoc:<d1>:<d2> with a numeric measure, impact:<outcome>:<disease> or inter:<outcome>:<d1>:<d2>).
```

## The global dairy example

``` r

gd <- cm_read_inputs(dir = system.file("extdata", "global_dairy_2024", package = "deconflate"))
gd
#> <cm_inputs>
#> <cm_model>
#>   Diseases: 12 (CK, CM, DA, DYS, LAM, MET, MF, OC, PTB, RP, SCK, SCM)
#>   Disease pairs: 66 [independent (default): 28; specified: 38]
#>   Outcomes: yield, fertility, culling
#>   Uncertain inputs: 70 (use $sampler with cm_monte_carlo())
cmp <- compare_methods(gd$model, methods = c("published", "simultaneous"))
head(cmp$impacts, 12)
#>    outcome disease      scale         raw    published simultaneous
#> 1    yield      CK proportion 0.004321944 0.0002473297 -0.053709141
#> 2    yield      CM proportion 0.032499000 0.0133194993 -0.003380725
#> 3    yield      DA proportion 0.028369300 0.0079378291 -0.026316316
#> 4    yield     DYS proportion 0.049190880 0.0356881898  0.043096749
#> 5    yield     LAM proportion 0.048061000 0.0253240737  0.019792779
#> 6    yield     MET proportion 0.056130850 0.0284095146  0.027909168
#> 7    yield      MF proportion 0.005365131 0.0006897532 -0.015399737
#> 8    yield      OC proportion 0.037478390 0.0263859343  0.032174755
#> 9    yield     PTB proportion 0.043000000 0.0322917066  0.039410642
#> 10   yield      RP proportion 0.041986640 0.0225771112  0.024733053
#> 11   yield     SCK proportion 0.083964720 0.0710307185  0.082785938
#> 12   yield     SCM proportion 0.062931840 0.0558961905  0.065761493
```

## Tables typed in R

Any table can be a data frame instead of a file:

``` r

inp2 <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
  associations = data.frame(disease1 = c("d1", "d2"), disease2 = c("d2", "d3"),
                            value = c(2, 3), measure = "OR"),
  impacts = data.frame(disease = c("d1", "d2", "d3"), outcome = "yield",
                       value = c(2.5, 5, 7.5), scale = "percent")
)
deconflate(inp2$model)$adjusted
#>   outcome disease   raw   adjusted      change      scale units direction
#> 1   yield      d1 0.025 0.02143250 -0.14269981 proportion  <NA>  decrease
#> 2   yield      d2 0.050 0.03387080 -0.32258408 proportion  <NA>  decrease
#> 3   yield      d3 0.075 0.06934209 -0.07543874 proportion  <NA>  decrease
```

## Uncertainty

Each row of `uncertainty.csv` gives one input a distribution. The key
names the input:

- `prob:<disease>`: the disease’s `value` (e.g. a prevalence or
  incidence rate);
- `assoc:<d1>:<d2>`: the association’s `value` (e.g. an odds ratio);
- `impact:<outcome>:<disease>`: the impact, on the scale it was entered
  (e.g. percent);
- `inter:<outcome>:<d1>:<d2>`: an interaction.

`dist` and the parameters `p1`-`p4`:

| `dist` | `p1` | `p2` | `p3` | `p4` |
|----|----|----|----|----|
| `fixed` | value |  |  |  |
| `normal` | mean | sd | lower bound (optional) | upper bound (optional) |
| `lognormal` | meanlog | sdlog |  |  |
| `lognormal_ci` | estimate | lower CI | upper CI | level (default 0.95) |
| `beta` | shape1 | shape2 | min (default 0) | max (default 1) |
| `pert` | min | mode | max | lambda (default 4) |
| `pert_mean` | min | mean | max | lambda (default 4) |
| `uniform` | min | max |  |  |

``` r

read.csv(file.path(dir, "uncertainty.csv"))
#>                  key   dist    p1     p2 p3 p4                            note
#> 1           prob:LAM   beta 78.29 227.42 NA NA            beta(shape1, shape2)
#> 2      assoc:LAM:SCK normal  2.01   0.20  0 NA normal(mean, sd) truncated at 0
#> 3   impact:yield:LAM normal  4.81   0.87 NA NA       normal(mean, sd), percent
#> 4   impact:yield:SCK normal  8.40   1.19 NA NA       normal(mean, sd), percent
#> 5 impact:culling:LAM normal  1.74   0.17  0 NA normal(mean, sd) truncated at 0
mc <- cm_monte_carlo(inp$sampler, 200, method = c("published", "simultaneous"), seed = 1)
compare_methods(mc)
#> <cm_comparison> methods: published, simultaneous
#> Monte Carlo: 200 draws (0 rejected); statistic: mean
#> 
#> Adjusted impacts:
#>  outcome disease unit raw_mean published simultaneous
#>    yield     LAM    %     4.83      3.21         2.92
#>    yield     SCK    %     8.52      7.64         7.94
#>    yield     MET    %     5.61      3.52         3.13
#>  culling     LAM   HR     1.75      1.57         1.55
#>  culling     SCK   HR     1.92      1.81         1.80
#>  culling     MET   HR     1.50      1.26         1.14
```

[`cm_dist_table()`](https://rasmussenphilip.github.io/deconflate/reference/cm_dist_table.md)
builds the distributions alone, if you want to pass them to
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
yourself.
