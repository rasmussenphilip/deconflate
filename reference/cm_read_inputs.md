# Read model inputs from CSV files or data frames

Builds a population
([`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)),
one analysis per impact table
([`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)),
optionally a hazard-ratio model
([`cm_hr_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_hr_model.md))
and, if uncertainty tables are given, a Monte Carlo sampler. Each table
can be a path to a CSV file or a data frame typed in R; alternatively,
`dir` names a folder of CSV files (see Files).
[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set to start from.

## Usage

``` r
cm_read_inputs(
  diseases = NULL,
  associations = NULL,
  three_way = NULL,
  impacts = NULL,
  interactions = NULL,
  uncertainty = NULL,
  hazard_ratios = NULL,
  dir = NULL,
  missing_associations = c("independent", "unknown"),
  adjusted_associations = c("error", "use_as_marginal")
)
```

## Arguments

- diseases, associations, three_way, hazard_ratios:

  Paths to CSV files or data frames.

- impacts, interactions, uncertainty:

  A path or data frame, or a named list of them (one per analysis; see
  Files).

- dir:

  Optional folder with CSV files named as in Files. Arguments given
  explicitly take precedence over files.

- missing_associations, adjusted_associations:

  Passed to
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).

## Value

A `cm_inputs` list with `population`, `analyses` (a
[`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
object, or `NULL` without impact tables), `model` (the
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
when there is exactly one analysis), `hr_model` (or `NULL`), `sampler`
(a
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
for one analysis, a
[`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md)
for several, or `NULL`), `tables` (as read; columns of CSV files are
read as text and converted by the checks) and `problems` (notes only,
since errors stop).

## Details

All tables are checked before anything is built, and every problem is
reported at once, with its table, row and column (row 1 is the first row
below the header).
[`cm_check_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_check_inputs.md)
runs the same checks without stopping.

## Files

In a folder (`dir`), files are found by name:

- `diseases.csv` (required), `associations.csv`, `three_way.csv`;

- one impact table per analysis: `impacts_<analysis>.csv`, e.g.
  `impacts_yield.csv` and `impacts_fertility.csv`. A single
  `impacts.csv` is an analysis named `impacts`. All analyses share the
  diseases and associations;

- optional interactions per analysis: `interactions_<analysis>.csv` (or
  `interactions.csv` when there is one analysis);

- optional uncertainty: `uncertainty.csv` (shared: disease
  probabilities, associations, and impacts keyed by analysis) and/or
  `uncertainty_<analysis>.csv` (that analysis's impacts and
  interactions);

- optional `hazard_ratios.csv` (culling or mortality hazard ratios, for
  [`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md)).

In R, `impacts`, `interactions` and `uncertainty` can be one table or a
named list of tables, named after the analyses (`uncertainty` may also
have an element named `shared`).

## Columns

Column names are not case-sensitive. Optional columns can be left out or
left empty; unrecognised columns are ignored (with a note).

**diseases**: `id`, `value` (required); `type` (`prevalence` (default),
`probability` or `incidence_rate`, converted with `1 - exp(-value)`),
`time_horizon`, `reference_population`, `source`. Ids must not contain
`|`, `;` or `:`, and `all` is reserved.

**associations**: `disease1`, `disease2` (required); `value`, `measure`
(`OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 \| disease2)),
`phi`, `table`, `independent` or `unknown`), `n11`, `n10`, `n01`, `n00`
(counts, for `measure = table`; a zero cell gets 0.5 added to every
cell, which is reported), `adjusted` (TRUE/FALSE: covariate-adjusted
measures are rejected unless
`adjusted_associations = "use_as_marginal"`), `adjusted_for`, `source`.
Pairs that are not listed are independent (or unknown, see
`missing_associations`).

**three_way**: `disease1`, `disease2`, `disease3`, `ratio` (required);
`source`. See
[`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md).

**impacts** (one table per analysis): `disease`, `value` (required);
`estimand` (`crude` (default) or `adjusted_linear`), `adjusted_for` (ids
separated by `;`, or `all`; only with `adjusted_linear`), `source`,
`label`, `units` (one label and one unit per table). Every disease needs
a row (use 0 for no impact). Values are in any units, the same for every
row; results come back in those units.

**interactions**: `disease1`, `disease2`, `value` (required); `source`.

**hazard_ratios**: `disease`, `value` (required); `estimand` (`crude`
(default) or `adjusted`), `adjusted_for`, `source`. Every disease needs
a row (use 1 for no effect).

**uncertainty**: `key`, `dist` (required); `p1`-`p4` (parameters),
`note`. Keys name the input: `prob:<disease>`, `assoc:<d1>:<d2>`,
`impact:<analysis>:<disease>` and `inter:<analysis>:<d1>:<d2>`. In an
analysis's own file (`uncertainty_<analysis>.csv`), or when there is
only one analysis, the analysis can be left out: `impact:<disease>`,
`inter:<d1>:<d2>`. Values are drawn on the scale the input was entered
on (e.g. an incidence rate, or a percent impact). Distributions and
their parameters:

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
#> Wrote 8 files to /tmp/RtmpZ17cjz/deconflate-inputs
inp <- cm_read_inputs(dir = dir)
inp
#> <cm_inputs>
#> <cm_population>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
#>   Analyses: calving_interval, yield
#>   Hazard ratios: yes (use $hr_model with deconflate_hr())
#>   Uncertain inputs: batch sampler over 2 analyses (use $sampler with cm_monte_carlo())
deconflate(inp$analyses)
#> <cm_results> 2 analyses: calving_interval, yield
#> 
#> == calving_interval ==
#> <cm_result> method: simultaneous; calving interval increase [days]
#> 
#>  disease raw adjusted  change
#>      LAM  12    8.922 -0.2565
#>      SCK   4    1.935 -0.5162
#>      MET  18   14.044 -0.2198
#> 
#> Raw sum: 6.325; adjusted total: 4.373
#> Diagnostics: residual 0.00e+00, condition number 2, sign changes 0
#> 
#> == yield ==
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0
#> 

# Example files shipped with the package: the 2024 global dairy inputs,
# and a set with deliberate errors
gd <- cm_read_inputs(dir = system.file("extdata", "global_dairy_2024", package = "deconflate"))
gd$analyses
#> <cm_analyses> 2 analyses on one population
#>   Diseases: 12 (CK, CM, DA, DYS, LAM, MET, MF, OC, PTB, RP, SCK, SCM)
#>   Disease pairs: 66 [independent (default): 28; specified: 38]
#>   - fertility [% increase]
#>   - yield [% decrease]
cm_check_inputs(dir = system.file("extdata", "example_with_errors", package = "deconflate"))
#> Found 23 problem(s) in the inputs:
#>   interactions_culling: There is no impact table for analysis 'culling'.
#>   diseases, row 5, column 'id': Duplicate disease id 'SCK'.
#>   diseases, row 3, column 'value': Gives a probability of 1.1; it must be strictly between 0 and 1.
#>   associations, row 3, column 'disease2': Unknown disease 'CM' (not in the diseases table).
#>   associations, row 4: Duplicate pair SCK:LAM.
#>   associations, row 5, column 'value': Missing value for measure OR.
#>   associations, row 6, column 'adjusted': A covariate-adjusted measure is not a marginal 2x2 association. Use the crude measure, or set adjusted_associations = 'use_as_marginal' to use it as an approximation.
#>   hazard_ratios, row 3, column 'value': Hazard ratios must be positive.
#>   impacts_fertility, column 'units': One units per table is allowed (found: % increase, days).
#>   impacts_yield, row 2, column 'value': '8,40' is not a number.
#>   impacts_yield, row 3, column 'adjusted_for': adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.
#>   impacts_yield: No impact for: MET. Add a row with value 0 for no impact.
#>   uncertainty, row 2, column 'dist': 'normal' needs p1 (mean) and p2 (sd).
#>   uncertainty, row 3, column 'key': Key 'assoc:LAM:RP' does not match an association with a numeric measure (expected assoc:<d1>:<d2>).
#>   uncertainty, row 4, column 'key': Key 'impact:yield:CM': unknown disease 'CM'.
#>   uncertainty, row 5, column 'key': Key 'impact:culling:LAM': unknown analysis 'culling'.
#>   uncertainty, row 6, column 'key': Key 'impact:LAM': name the analysis (impact:<analysis>:<disease>), since there are several.
#>   uncertainty, row 7, column 'dist': PERT needs min <= mode <= max and min < max.
#>   uncertainty, row 9, column 'key': Key 'inter:yield:LAM:SCK' does not match an interaction of analysis 'yield'.
#>   uncertainty, row 10, column 'key': Key 'milk_price' must start with prob:, assoc:, impact: or inter:.
#>   uncertainty, row 11, column 'dist': Unknown distribution 'gamma' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).
#>   uncertainty_fertility, row 3, column 'key': Disease probabilities are shared by all analyses; put prob: keys in uncertainty.csv.
#>   uncertainty_fertility, row 2, column 'key': The same input is given more than once. 

# The same kind of tables typed in R
inp2 <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
  associations = data.frame(disease1 = "d1", disease2 = "d2", value = 2),
  impacts = list(yield = data.frame(disease = c("d1", "d2"), value = c(2.5, 5),
                                    units = "% of yield"))
)
deconflate(inp2$model)
#> <cm_result> method: simultaneous; yield [% of yield]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    1.989 -0.20442
#>       d2 5.0    4.852 -0.02957
#> 
#> Raw sum: 1; adjusted total: 0.9267
#> Diagnostics: residual 0.00e+00, condition number 1.2, sign changes 0
```
