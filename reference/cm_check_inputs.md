# Check model input tables

Runs the checks of
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
and returns every problem found, without stopping.

## Usage

``` r
cm_check_inputs(
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

A data frame (class `cm_problems`) with `table`, `row`, `column`,
`severity` (`"error"` or `"note"`) and `problem`. No rows means the
inputs can be read.

## Examples

``` r
cm_check_inputs(
  diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 1.5)),
  impacts = data.frame(disease = c("d1", "d3"), outcome = "yield", value = c(2, "x"))
)
#> Found 4 problem(s) in the inputs:
#>   diseases, row 2, column 'value': Gives a probability of 1.5; it must be strictly between 0 and 1.
#>   impacts, row 2, column 'value': 'x' is not a number.
#>   impacts, row 2, column 'disease': Unknown disease 'd3' (not in the diseases table).
#>   impacts: Outcome 'yield' has no impact for: d2. Add a row with value 0 for no impact. 
```
