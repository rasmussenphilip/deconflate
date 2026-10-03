# Read model inputs from CSV files or data frames

Builds a
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
(and, if an uncertainty table is given, a
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md))
from up to five tables. Each table can be a path to a CSV file or a data
frame typed in R; alternatively, `dir` names a folder containing
`diseases.csv`, `associations.csv`, `impacts.csv`, `interactions.csv`
and `uncertainty.csv` (only `diseases.csv` is required).
[`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)
writes an example set to start from.

## Usage

``` r
cm_read_inputs(
  diseases = NULL,
  associations = NULL,
  impacts = NULL,
  interactions = NULL,
  uncertainty = NULL,
  dir = NULL,
  missing_associations = c("independent", "unknown"),
  outcome_correlation = NULL
)
```

## Arguments

- diseases, associations, impacts, interactions, uncertainty:

  Paths to CSV files or data frames. `NULL` tables are taken from `dir`
  (if given) or left out.

- dir:

  Optional folder with the CSV files named as above.

- missing_associations:

  Passed to
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- outcome_correlation:

  Passed to
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

## Value

A `cm_inputs` list with `model`, `sampler` (or `NULL`), `tables` (the
tables as read) and `problems` (notes only, since errors stop).

## Details

All tables are checked before anything is built, and every problem is
reported at once, with its table, row and column (row 1 is the first row
below the header).
[`cm_check_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_check_inputs.md)
runs the same checks without stopping.

## Columns

Column names are not case-sensitive. Optional columns can be left out or
left empty; unrecognised columns are ignored (with a note).

**diseases**: `id`, `value` (required); `type` (`prevalence` (default),
`probability` or `incidence_rate`, converted with `1 - exp(-value)`),
`time_horizon`, `reference_population`, `source`.

**associations**: `disease1`, `disease2` (required); `value`, `measure`
(`OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 \| disease2)),
`phi`, `table`, `independent` or `unknown`), `n11`, `n10`, `n01`, `n00`
(counts, for `measure = table`), `adjusted` (TRUE/FALSE),
`adjusted_for`, `source`. Pairs that are not listed are independent (or
unknown, see `missing_associations`).

**impacts**: `disease`, `outcome`, `value` (required); `scale`
(`proportion` (default), `percent`, `absolute` or `hazard_ratio`),
`units` (required for `absolute`), `direction` (`decrease` (default) or
`increase`; ignored for hazard ratios), `adjusted_for` (ids separated by
`;`), `source`. Every disease needs one row per outcome (use 0 for no
impact).

**interactions**: `disease1`, `disease2`, `value`, `outcome` (required);
`scale` (`proportion` (default) or `percent`), `source`.

**uncertainty**: `key`, `dist` (required); `p1`-`p4` (parameters),
`note`. Keys name the input: `prob:<disease>`, `assoc:<d1>:<d2>`,
`impact:<outcome>:<disease>` or `inter:<outcome>:<d1>:<d2>`. Values are
drawn on the scale the input was entered on (e.g. an incidence rate, or
a percent impact). Distributions and their parameters:

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
#> Wrote 5 files to /tmp/RtmpP1FJSR/deconflate-inputs
inp <- cm_read_inputs(dir = dir)
inp
#> <cm_inputs>
#> <cm_model>
#>   Diseases: 3 (LAM, SCK, MET)
#>   Disease pairs: 3 [specified: 3]
#>   Outcomes: yield, culling
#>   Uncertain inputs: 5 (use $sampler with cm_monte_carlo())
deconflate(inp$model, method = "global")
#> Warning: Outcome 'culling': adjusted hazard ratios cross 1 for MET. The raw hazard ratios are smaller than the associated diseases alone would produce; check whether these estimates were already adjusted for co-diseases or come from populations with different comorbidity patterns.
#> <cm_result> method: global
#> 
#>  outcome disease    raw adjusted   change
#>    yield     LAM 0.0481  0.02871 -0.40320
#>    yield     SCK 0.0840  0.07819 -0.06915
#>    yield     MET 0.0561  0.03170 -0.43496
#>  culling     LAM 1.7400  1.62941 -0.06356
#>  culling     SCK 1.9200  1.80906 -0.05778
#>  culling     MET 1.1200  0.82814 -0.26059
#> 
#> Diagnostics:
#>  outcome max_reconstruction_residual n_sign_changes condition_number
#>    yield                    6.94e-18              0                2
#>  culling                    1.22e-13              1                2

# The same tables typed in R
inp2 <- cm_read_inputs(
  diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
  associations = data.frame(disease1 = "d1", disease2 = "d2", value = 2),
  impacts = data.frame(disease = c("d1", "d2"), outcome = "yield",
                       value = c(2.5, 5), scale = "percent")
)
deconflate(inp2$model)
#> <cm_result> method: simultaneous
#> 
#>  outcome disease   raw adjusted   change
#>    yield      d1 0.025  0.01989 -0.20442
#>    yield      d2 0.050  0.04852 -0.02957
#> 
#> Diagnostics:
#>  outcome max_reconstruction_residual n_sign_changes condition_number
#>    yield                    6.94e-18              0              1.2
```
