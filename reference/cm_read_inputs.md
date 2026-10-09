# Read model inputs from CSV files or data frames

Builds a model
([`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md))
from one table of each kind: the diseases, their associations, one
impact table and, optionally, interactions and three-way terms. Each
table can be a path to a CSV file (with any name) or a data frame typed
in R. Values with a distribution (see Uncertainty) are kept in the
model, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
uses them to compute intervals.
[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set to start from.

## Usage

``` r
cm_read_inputs(
  diseases,
  associations = NULL,
  impacts,
  interactions = NULL,
  three_way = NULL,
  adjusted_associations = c("error", "use_as_marginal")
)
```

## Arguments

- diseases, associations, impacts, interactions, three_way:

  Paths to CSV files (any names) or data frames. `diseases` and
  `impacts` are required; `associations` can be left out only for the
  sensitivity tools
  ([`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  needs at least one association).

- adjusted_associations:

  Passed to
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
with the distributions of uncertain values in `$distributions` and any
notes from the checks in attribute `"problems"`.

## Details

All tables are checked before anything is built, and every problem is
reported at once, with its table, row and column (row 1 is the first row
below the header).
[`cm_check_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_check_inputs.md)
runs the same checks without stopping.

For several outcomes (e.g. milk yield and calving interval), keep one
impact table per outcome and read and adjust each in turn: the disease
and association tables can be shared.

## Columns

Column names are not case-sensitive. Optional columns can be left out or
left empty; unrecognised columns are ignored (with a note). Every table
can also have a free-text `note` column.

**diseases**: `id`, `value` (required); `type` (`prevalence` (default),
`probability` or `incidence_rate`, converted with `1 - exp(-value)`),
`time_horizon`, `reference_population`, `source`. Ids must not contain
`|`, `;` or `:`, and `all` is reserved.

**associations**: `disease1`, `disease2`, `value` (required); `measure`
(`OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 \| disease2)) or
`phi`), `adjusted` (TRUE/FALSE: covariate-adjusted measures are rejected
unless `adjusted_associations = "use_as_marginal"`), `adjusted_for`,
`source`. Pairs without a row are unknown: the global model fills in
their association from the others. An odds ratio of 1 states that two
diseases are unrelated.

**impacts**: `disease`, `value` (required); `estimand`, `adjusted_for`,
`measure`, `source`, `label`, `units` (one label and one unit per
table). Every disease needs a row. Without a `measure` column the
impacts are additive, in any units (the same for every row; results come
back in those units): `estimand` is `crude` (default) or
`adjusted_linear` (with `adjusted_for`: ids separated by `;`, or `all`),
and a disease with no impact has the value 0. With a `measure` column
the impacts are event impacts (see
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)):
every row needs a measure (`HR`, `rate_ratio`, `RR`, `OR` or `RD`) and
an estimand (`snapshot_crude`, or `snapshot_stratified` with
`adjusted_for`), and they are adjusted with
`deconflate(..., event_model = TRUE)`.

**interactions** (additive impacts only): `disease1`, `disease2`,
`value` (required); `source`.

**three_way**: `disease1`, `disease2`, `disease3`, `ratio` (required);
`source`. See
[`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md).

## Uncertainty

Every table can have the columns `dist` and `p1`-`p4`. A row with a
`dist` has an uncertain value:
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
draws it from its distribution (`n_draws`), while the point value
(`value`, or `ratio` for three-way terms) gives the central estimate. A
row with an empty `dist` is a fixed point value, so point values and
uncertain values can be mixed freely. Values are drawn on the scale they
were entered on (e.g. an incidence rate, an odds ratio, a hazard ratio
or a percent impact). A point value outside the support of its
distribution is noted. Distributions and their parameters:

|  |  |  |  |  |
|----|----|----|----|----|
| `dist` | `p1` | `p2` | `p3` | `p4` |
| `fixed` | value |  |  |  |
| `normal` | mean | sd | lower bound (optional) | upper bound (optional) |
| `lognormal` | meanlog | sdlog |  |  |
| `lognormal_ci` | estimate | lower CI | upper CI | level (default 0.95) |
| `beta` | shape1 | shape2 | min (default 0) | max (default 1) |
| `pert` | min | mode | max | lambda (default 4) |
| `pert_mean` | min | mean | max | lambda (default 4) |
| `uniform` | min | max |  |  |

## Examples

``` r
dir <- file.path(tempdir(), "deconflate-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 6 files to /tmp/RtmpNxOcWp/deconflate-inputs
m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                    associations = file.path(dir, "associations.csv"),
                    impacts = file.path(dir, "yield.csv"))
m
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 (3 with an association, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Estimands: crude: 3
#>   Uncertain inputs (with a distribution): 4
deconflate(m, n_draws = 0)
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0

# Example files shipped with the package
ex <- system.file("extdata", "five_diseases", package = "deconflate")
five <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
                       associations = file.path(ex, "associations.csv"),
                       impacts = file.path(ex, "impacts_yield.csv"),
                       three_way = file.path(ex, "three_way.csv"))
five
#> <cm_model>
#>   Diseases: 5 (LAM, MAS, MET, SCK, RP)
#>   Disease pairs: 10 (10 with an association, 0 unknown)
#>   Three-way terms: 1
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Estimands: crude: 5
#>   Uncertain inputs (with a distribution): 16
err <- system.file("extdata", "example_with_errors", package = "deconflate")
cm_check_inputs(diseases = file.path(err, "diseases.csv"),
                associations = file.path(err, "associations.csv"),
                impacts = file.path(err, "impacts_yield.csv"))
#> Found 18 problem(s) in the inputs:
#>   diseases, row 5, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   diseases, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   diseases, row 4, column 'dist': Parameters p1-p4 are given without a distribution (dist).
#>   associations: Contingency-table counts (n11, n10, n01, n00) are no longer read: compute each table's odds ratio, n11 * n00 / (n10 * n01), and enter it in `value` with measure OR.
#>   associations, row 3, column 'disease2': Unknown disease 'CM' (not in the diseases table).
#>   associations, row 9, column 'disease2': Unknown disease 'MAS' (not in the diseases table).
#>   associations, row 4: Duplicate pair SCK:LAM.
#>   associations, row 7, column 'measure': measure = 'table' is no longer used: compute the odds ratio from the counts, n11 * n00 / (n10 * n01), and enter it with measure OR.
#>   associations, row 8, column 'measure': measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).
#>   associations, row 9, column 'measure': measure = 'unknown' is no longer used: leave the pair out. Pairs without a row are unknown, and the global model fills in their association.
#>   associations, row 5, column 'value': Missing value. Leave the pair out (no row) if its association is unknown; an odds ratio of 1 states that the diseases are unrelated.
#>   associations, row 6, column 'adjusted': A covariate-adjusted measure is not a marginal 2x2 association. Use the crude measure, or set adjusted_associations = 'use_as_marginal' to use it as an approximation.
#>   associations, row 2, column 'dist': PERT needs min <= mode <= max and min < max.
#>   impacts, row 2, column 'value': '8,40' is not a number.
#>   impacts, row 3, column 'adjusted_for': adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.
#>   impacts, column 'units': One units per table is allowed (found: % decrease, days).
#>   impacts, row 3, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
#>   impacts, row 4, column 'dist' (note): The point value 5.61 lies outside the distribution's support [6, 8]. 

# The same kind of tables typed in R
m2 <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
  associations = data.frame(disease1 = "d1", disease2 = "d2", value = 2),
  impacts = data.frame(disease = c("d1", "d2"), value = c(2.5, 5), units = "% of yield")
)
deconflate(m2)
#> <cm_result> method: simultaneous [% of yield]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    1.989 -0.20442
#>       d2 5.0    4.852 -0.02957
#> 
#> Raw sum: 1; adjusted total: 0.9267
#> Diagnostics: residual 0.00e+00, condition number 1.2, sign changes 0
#> 
#> Notes:
#> * No input has a distribution, so no draws were run: the results are point
#>   estimates.
```
