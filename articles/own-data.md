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
#> Wrote 5 files to /tmp/RtmpHL3xwP/my-inputs
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
#> 6     MET culling  1.12 hazard_ratio              NA           NA
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
#> Warning: Outcome 'culling': adjusted hazard ratios cross 1 for MET. The raw
#> hazard ratios are smaller than the associated diseases alone would produce;
#> check whether these estimates were already adjusted for co-diseases or come
#> from populations with different comorbidity patterns.
res$adjusted[, c("outcome", "disease", "raw", "adjusted")]
#>   outcome disease    raw   adjusted
#> 1   yield     LAM 0.0481 0.02870596
#> 2   yield     SCK 0.0840 0.07819114
#> 3   yield     MET 0.0561 0.03169881
#> 4 culling     LAM 1.7400 1.62940621
#> 5 culling     SCK 1.9200 1.80906323
#> 6 culling     MET 1.1200 0.82814177
```

Compare the methods side by side:

``` r

compare_methods(inp$model)
#> <cm_comparison> methods: published, simultaneous, global
#> 
#> Adjusted impacts:
#>  outcome disease unit  raw published simultaneous global
#>    yield     LAM    % 4.81      3.16        2.870  2.870
#>    yield     SCK    % 8.40      7.51        7.820  7.820
#>    yield     MET    % 5.61      3.52        3.170  3.170
#>  culling     LAM   HR 1.74      1.60        1.630  1.630
#>  culling     SCK   HR 1.92      1.82        1.820  1.810
#>  culling     MET   HR 1.12      1.03        0.833  0.828
#> 
#> Sign changes (adjusted impact on the other side of zero, or of 1 for hazard ratios):
#>   simultaneous: culling MET
#>   global: culling MET
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
stops with the same list when there are errors.

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
#>    yield     LAM    %     4.83      3.21        2.920
#>    yield     SCK    %     8.52      7.64        7.940
#>    yield     MET    %     5.61      3.52        3.130
#>  culling     LAM   HR     1.75      1.61        1.650
#>  culling     SCK   HR     1.92      1.82        1.820
#>  culling     MET   HR     1.12      1.03        0.831
```

[`cm_dist_table()`](https://rasmussenphilip.github.io/deconflate/reference/cm_dist_table.md)
builds the distributions alone, if you want to pass them to
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
yourself.
