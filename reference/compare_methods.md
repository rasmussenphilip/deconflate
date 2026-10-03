# Compare adjustment methods

Runs several adjustment methods on the same inputs and tabulates the
results side by side:

- `"published"`: the simple proportional approximation used in Rasmussen
  et al. (2022, 2024) (eq. 16; for hazard ratios, `HR - 1` as in 2024);

- `"simultaneous"`: the exact solution of the system of equations;

- `"global"`: the iterative (maximum-entropy) model of disease
  combinations, which can also include interactions.

## Usage

``` r
compare_methods(x, ...)

# S3 method for class 'cm_model'
compare_methods(
  x,
  methods = c("published", "simultaneous", "global"),
  economics = NULL,
  ...
)

# S3 method for class 'cm_mc'
compare_methods(x, stat = c("mean", "median", "trimmed_mean"), ...)
```

## Arguments

- x:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  or a `cm_mc` object from
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  run with more than one method.

- ...:

  Passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  (models) or
  [`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md)
  (Monte Carlo runs).

- methods:

  Methods to compare (models only).

- economics:

  Optional economics list (see
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md);
  for hazard-ratio outcomes,
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)).
  Adds gaps and values per method to `totals`.

- stat:

  For Monte Carlo runs: the statistic to compare, `"mean"`, `"median"`
  or `"trimmed_mean"`.

## Value

A `cm_comparison` object with:

- `impacts`: one row per outcome and disease, with the raw impact and
  one column of adjusted impacts per method;

- `change`: the same layout with the relative change from the raw
  impact;

- `long`: the same information in long format, with sign-change flags
  (Monte Carlo: the summary statistics and stability flags per method);

- `totals`: per outcome and method, the expected loss before and after
  adjustment and, with `economics`, the gap and its value;

- `total_value`: total value per method (with `economics`);

- `diagnostics`:
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  diagnostics per method (models);

- `failed`: methods that could not be run, with the reason.

## Details

For a model, every method is run on the model's inputs. For a Monte
Carlo run made with several methods
(`cm_monte_carlo(..., method = c(...))`), the methods were applied to
identical draws, and the table compares their summaries.

## Examples

``` r
cmp <- compare_methods(example_uk_dairy_2022(), economics = uk_dairy_2022_economics())
cmp
#> <cm_comparison> methods: published, simultaneous, global
#> 
#> Adjusted impacts:
#>    outcome disease unit   raw published simultaneous    global
#>      yield      CO    %  0.00    0.0000    -6.06e-01 -6.06e-01
#>      yield      DA    %  4.04    2.3700     2.20e+00  2.20e+00
#>      yield     DYS    %  4.05    2.9400     2.95e+00  2.95e+00
#>      yield     FAS    %  7.33    7.3300     7.33e+00  7.33e+00
#>      yield     GIN    %  3.28    3.2800     3.28e+00  3.28e+00
#>      yield     LAM    %  5.54    4.8300     5.03e+00  5.03e+00
#>      yield     MAS    %  4.57    3.7600     4.08e+00  4.08e+00
#>      yield     MET    %  3.95    2.4100     1.86e+00  1.86e+00
#>      yield      MF    %  0.41    0.0999    -5.14e-01 -5.14e-01
#>      yield     NEO    %  4.20    4.2000     4.20e+00  4.20e+00
#>      yield     PTB    %  5.90    4.4300     4.13e+00  4.13e+00
#>      yield      RP    %  7.38    6.1000     6.58e+00  6.58e+00
#>      yield     SCK    %  3.05    1.9200     1.64e+00  1.64e+00
#>  fertility      CO    % 11.30   11.1000     1.15e+01  1.15e+01
#>  fertility      DA    %  0.00    0.0000    -3.94e-01 -3.94e-01
#>  fertility     DYS    %  6.96    6.0200     6.08e+00  6.08e+00
#>  fertility     FAS    %  0.00    0.0000    -2.78e-17 -7.81e-16
#>  fertility     GIN    %  1.20    1.2000     1.20e+00  1.20e+00
#>  fertility     LAM    % 12.50   11.9000     1.24e+01  1.24e+01
#>  fertility     MAS    %  0.00    0.0000    -9.41e-01 -9.41e-01
#>  fertility     MET    %  4.74    4.1000     4.62e+00  4.62e+00
#>  fertility      MF    %  0.00    0.0000     2.95e-01  2.95e-01
#>  fertility     NEO    %  7.21    7.2100     7.21e+00  7.21e+00
#>  fertility     PTB    %  5.79    3.8600     3.04e+00  3.04e+00
#>  fertility      RP    %  2.74    1.6800     1.17e+00  1.17e+00
#>  fertility     SCK    %  1.50    0.5500    -9.64e-01 -9.64e-01
#>    culling      CO    %  0.00    0.0000    -2.68e+00 -2.68e+00
#>    culling      DA    % 31.40   21.9000     2.25e+01  2.25e+01
#>    culling     DYS    % 14.20   11.9000     1.28e+01  1.28e+01
#>    culling     FAS    %  0.00    0.0000    -4.75e-16 -2.63e-15
#>    culling     GIN    %  0.00    0.0000    -2.61e-16  2.91e-15
#>    culling     LAM    % 25.60   22.6000     2.38e+01  2.38e+01
#>    culling     MAS    % 21.30   16.9000     1.78e+01  1.78e+01
#>    culling     MET    % 17.40   12.3000     1.26e+01  1.26e+01
#>    culling      MF    % 20.60   15.6000     1.60e+01  1.60e+01
#>    culling     NEO    %  9.89    9.8900     9.89e+00  9.89e+00
#>    culling     PTB    % 19.60   13.4000     1.15e+01  1.15e+01
#>    culling      RP    %  0.00    0.0000    -5.17e+00 -5.17e+00
#>    culling     SCK    % 15.70    9.6700     7.61e+00  7.61e+00
#> 
#> Sign changes (adjusted impact on the other side of zero, or of 1 for hazard ratios):
#>   simultaneous: yield MF, fertility SCK
#>   global: yield MF, fertility SCK
#> 
#> Totals:
#>    outcome       method raw_loss adjusted_loss observed disease_free     gap
#>    culling    published  0.24990       0.19760       27        22.55   4.454
#>  fertility    published  0.07573       0.06907      401       375.10  25.910
#>      yield    published  0.07168       0.06046     8737      9299.00 562.200
#>    culling simultaneous  0.24990       0.19410       27        22.61   4.389
#>  fertility simultaneous  0.07573       0.06448      401       376.70  24.290
#>      yield simultaneous  0.07168       0.05982     8737      9293.00 555.900
#>    culling       global  0.24990       0.19410       27        22.61   4.389
#>  fertility       global  0.07573       0.06448      401       376.70  24.290
#>      yield       global  0.07168       0.05982     8737      9293.00 555.900
#>   value
#>   59.48
#>  101.80
#>  169.90
#>   58.61
#>   95.43
#>  168.00
#>   58.61
#>   95.43
#>  168.00
#> 
#> Total value: published 402.25; simultaneous 393.13; global 393.13
cmp$totals
#>     outcome       method  raw_loss adjusted_loss observed disease_free
#> 1   culling    published 0.2498891    0.19756853       27     22.54568
#> 2 fertility    published 0.0757340    0.06906756      401    375.09323
#> 3     yield    published 0.0716780    0.06045754     8737   9299.20724
#> 4   culling simultaneous 0.2498891    0.19411785       27     22.61083
#> 5 fertility simultaneous 0.0757340    0.06448457      401    376.70814
#> 6     yield simultaneous 0.0716780    0.05982142     8737   9292.91542
#> 7   culling       global 0.2498891    0.19411785       27     22.61083
#> 8 fertility       global 0.0757340    0.06448457      401    376.70814
#> 9     yield       global 0.0716780    0.05982142     8737   9292.91542
#>          gap     value
#> 1   4.454317  59.48117
#> 2  25.906774 101.77735
#> 3 562.207239 169.89903
#> 4   4.389166  58.61117
#> 5  24.291862  95.43301
#> 6 555.915420 167.99764
#> 7   4.389166  58.61117
#> 8  24.291862  95.43301
#> 9 555.915420 167.99764
```
