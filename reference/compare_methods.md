# Compare adjustment methods

Runs several adjustment methods on the same inputs and tabulates the
results side by side:

- `"simultaneous"`: the exact solution of the system of equations;

- `"global"`: the iterative (maximum-entropy) model of disease
  combinations, which can also include interactions;

- `"published"`: the proportional approximation used in Rasmussen et al.
  (2022, 2024) (eq. 16 of the 2022 paper). It is not a method of
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md);
  it is kept here (and in
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  and the `reproduce_*()` functions) to compare earlier results with the
  exact solution. It applies to crude estimates only, can mask
  incompatible inputs, and is undefined when its denominator is zero.

## Usage

``` r
compare_methods(x, ...)

# S3 method for class 'cm_model'
compare_methods(x, methods = c("published", "simultaneous", "global"), ...)

# S3 method for class 'cm_analyses'
compare_methods(x, methods = c("published", "simultaneous", "global"), ...)

# S3 method for class 'cm_hr_model'
compare_methods(
  x,
  methods = c("published", "first_order", "snapshot"),
  overall_risk = NULL,
  ...
)

# S3 method for class 'cm_mc'
compare_methods(x, stat = c("mean", "median", "trimmed_mean"), ...)
```

## Arguments

- x:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md),
  [`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md),
  or a `cm_mc` object from
  [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  run with more than one method.

- ...:

  Passed to the adjustment (e.g. `joint`, `feasibility`, or arguments of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md));
  for Monte Carlo runs, to
  [`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md).

- methods:

  Methods to compare.

- overall_risk:

  For hazard-ratio models: optional overall risk, adding the
  attributable risk per method (see
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)).

- stat:

  For Monte Carlo runs: `"mean"`, `"median"` or `"trimmed_mean"`.

## Value

A `cm_comparison` object with `impacts` (raw and adjusted values, one
column per method), `change` (relative change from raw), `long` (long
format with sign-change flags, or the Monte Carlo summaries), `totals`
(naive and adjusted aggregate per method; attributable risk for
hazard-ratio models with `overall_risk`), `diagnostics`, `failed`
(methods that could not be run, gave an undefined (non-finite) result,
or whose attributable risk could not be computed, with reasons) and
`methods`. For a model, `undefined` keeps the results with non-finite
values for inspection; they are not among the estimates.

## Details

For a model, every method is run on the model's inputs; methods that
cannot be run (e.g. the published approximation for adjusted estimands)
are reported with the reason. For
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md),
each analysis is compared and the tables are stacked. For a hazard-ratio
model, the methods of
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)
(`"snapshot"`, `"first_order"`) are compared with the published approach
of Rasmussen et al. (2024) (`"published"`: `HR - 1` adjusted with eq.
16). For a Monte Carlo run made with several methods, the methods were
applied to identical draws, and their summaries are compared.

## Examples

``` r
compare_methods(example_supplement())
#> <cm_comparison> methods: published, simultaneous, global
#> Units: %
#> 
#> Adjusted values:
#>  disease raw published simultaneous global
#>       d1 2.5      2.07         2.14   2.14
#>       d2 5.0      3.70         3.39   3.39
#>       d3 7.5      6.75         6.93   6.93
#> 
#> Totals:
#>        method raw_sum adjusted_total
#>     published     2.5          2.111
#>  simultaneous     2.5          2.109
#>        global     2.5          2.109
```
