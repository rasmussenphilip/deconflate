# Read model inputs from CSV files or data frames

Builds a model
([`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md))
from three tables: the diseases, their associations and one impact
table. Each table can be a path to a CSV file (with any name) or a data
frame typed in R. Values with a distribution (see Uncertainty) are kept
in the model, and
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
uses them to compute intervals.
[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set to start from.

## Usage

``` r
cm_read_inputs(diseases, associations = NULL, impacts, ...)
```

## Arguments

- diseases, associations, impacts:

  Paths to CSV files (any names) or data frames. `diseases` and
  `impacts` are required; `associations` can be left out only for the
  sensitivity tools
  ([`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  needs at least one association).

- ...:

  Not used (arguments of earlier versions give a message).

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
with the distributions of uncertain values in `$distributions` (keyed by
table and row, e.g. `"associations:3"`) and any notes from the checks in
attribute `"problems"`. When a check finds an error, the error (class
`deconflate_input_problems`) holds every problem found in `$problems`,
as
[`cm_check_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_check_inputs.md)
returns them. The notes are also shown as a message of class
`deconflate_input_notes`, which holds them in `$problems`;
[`suppressMessages()`](https://rdrr.io/r/base/message.html) hides it.

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
left empty; unrecognised columns are ignored (with a note). In every
table, `disease` (and `associated_with`) take a disease id or a
combination of ids separated by `;` (e.g. `MAS;SCK`), and `among` takes
a subgroup (`SCK`, `no SCK`, `SCK;no MET`, or `no others`).

**diseases** (one row per disease; whole-population values only):
`disease`, `value` (required); `measure` (`prevalence` (default),
`probability` or `incidence_rate`, converted with `1 - exp(-value)`),
`label` (full name), `time_horizon`, `reference_population`, `source`,
`note`. Ids must not contain `|`, `;` or `:` or start with `no `; `all`
and `others` are reserved.

**associations** (one row per estimate; any row about two or more
diseases): `disease`, `value` (required); `associated_with`, `among`,
`adjusted_for`, `measure` (`OR` (default), `RR`, `RD`, `phi`,
`prevalence` or `OR_ratio`), `study`, `source`, `note`. See
[`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
for what each row means. Pairs that no row mentions are unknown: the
model fills in their association. An odds ratio of 1 states that two
diseases are unrelated. Several rows about the same association are
pooled.

**impacts** (one row per estimate): `disease`, `value` (required);
`among`, `adjusted_for`, `measure`, `label`, `units` (one label and one
unit per table), `study`, `source`, `note`. Every disease needs a row of
its own; a row for a combination (e.g. `SCK;SCM`) adds an interaction to
the model, and a combination adjusted for `all` is the interaction
itself. See
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md).
Without a `measure` column the impacts are additive, in any units (the
same for every row; results come back in those units). With a `measure`
column they are event impacts: every row needs a measure (`HR`,
`rate_ratio`, `RR`, `OR` or `RD`), and they are adjusted with
`deconflate(..., event_model = TRUE)`.

Files in the layout of earlier versions (columns `id` and `type`,
`disease1` and `disease2`, `estimand`; three-way and interaction files)
are reported with what to change.

## Uncertainty

Every table can have the columns `dist` and `p1`-`p4`. A row with a
`dist` has an uncertain value:
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
draws it from its distribution (`n_draws`), while `value` gives the
central estimate. A row with an empty `dist` is a fixed point value
(and, when a disease or association has several rows, an exact one), so
point values and uncertain values can be mixed freely. Values are drawn
on the scale they were entered on (e.g. an incidence rate, an odds
ratio, a hazard ratio or a percent impact). A point value outside the
support of its distribution is noted. Distributions and their
parameters:

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
#> Wrote 4 files to /tmp/RtmpxOqF10/deconflate-inputs
m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                    associations = file.path(dir, "associations.csv"),
                    impacts = file.path(dir, "yield.csv"))
m
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Association rows: 3 (OR: 3)
#>   Disease pairs: 3 (3 with a crude pair row, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Impact rows: 3 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 4
deconflate(m, n_draws = 0)
#> <cm_result> pairwise shortcut; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 0.00e+00, condition number 2, sign changes 0

# Example files shipped with the package
ex <- system.file("extdata", "five_diseases", package = "deconflate")
five <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
                       associations = file.path(ex, "associations.csv"),
                       impacts = file.path(ex, "impacts_yield.csv"))
five
#> <cm_model>
#>   Diseases: 5 (LAM, MAS, MET, SCK, RP)
#>   Association rows: 12 (OR: 7; OR_ratio: 1; RD: 1; RR: 1; phi: 1; prevalence: 1)
#>   Disease pairs: 10 (10 with a crude pair row, 0 unknown)
#>   Impacts: milk yield loss [% of yield] (additive)
#>   Impact rows: 5 (0 adjusted, 0 within subgroups)
#>   Uncertain inputs (with a distribution): 17
err <- system.file("extdata", "example_with_errors", package = "deconflate")
cm_check_inputs(diseases = file.path(err, "diseases.csv"),
                associations = file.path(err, "associations.csv"),
                impacts = file.path(err, "impacts_yield.csv"))
#> Found 17 problem(s) in the inputs:
#>   diseases, row 5, column 'disease': Duplicate disease 'SCK': the diseases table has one row per disease (rows about several diseases or subgroups belong in the associations table).
#>   diseases, row 6, column 'disease': Disease id 'LAM;SCK' must not contain '|', ';' or ':' (a combination of diseases belongs in the associations table).
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   diseases, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   diseases, row 4, column 'dist': Parameters p1-p4 are given without a distribution (dist).
#>   associations, row 6, column 'measure': measure = 'cond_prob' is now 'prevalence': P(disease | associated_with) is the prevalence of `disease` among animals with `associated_with`.
#>   associations, row 7, column 'measure': measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).
#>   associations, row 4, column 'value': Missing value. Leave the pair out (no row) if its association is unknown; an odds ratio of 1 states that the diseases are unrelated.
#>   associations, row 3, column 'associated_with': unknown disease in `associated_with`: CM.
#>   associations, row 8, column 'among': `among` names LAM, which the row is about.
#>   associations, row 9: an adjusted association must be an odds ratio between two single diseases (measure OR).
#>   associations, row 10, column 'among': OR_ratio needs one disease in `among` (present, e.g. 'SCK', or absent, 'no SCK').
#>   associations, row 5, column 'adjusted_for' (note): 'parity', 'herd' are not diseases of the model: treated as factors outside the model (ignored; the row is used as crude).
#>   associations, row 2, column 'dist': PERT needs min <= mode <= max and min < max.
#>   impacts, row 2, column 'value': '8,40' is not a number.
#>   impacts, row 5, column 'disease': unknown disease in `disease`: CM.
#>   impacts, column 'units': One unit per table is allowed (found: % decrease, days).
#>   impacts, row 3, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
#>   impacts, row 4, column 'dist' (note): The point value 5.61 lies outside the distribution's support [6, 8]. 

# The same kind of tables typed in R
m2 <- cm_read_inputs(
  diseases = data.frame(disease = c("d1", "d2"), value = c(0.10, 0.15)),
  associations = data.frame(disease = "d1", associated_with = "d2", value = 2),
  impacts = data.frame(disease = c("d1", "d2"), value = c(2.5, 5), units = "% of yield")
)
deconflate(m2)
#> <cm_result> pairwise shortcut [% of yield]
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
