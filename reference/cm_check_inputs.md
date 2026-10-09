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

A data frame (class `cm_problems`) with `table`, `row`, `column`,
`severity` (`"error"` or `"note"`) and `problem`. No rows means the
inputs can be read.

## Examples

``` r
cm_check_inputs(
  diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 1.5)),
  impacts = data.frame(disease = c("d1", "d3"), value = c(2, "x"))
)
#> Found 4 problem(s) in the inputs:
#>   diseases, row 2, column 'value': Gives a probability of 1.5; it must be strictly between 0 and 1.
#>   impacts, row 2, column 'value': 'x' is not a number.
#>   impacts, row 2, column 'disease': Unknown disease 'd3' (not in the diseases table).
#>   impacts: No impact for: d2. Add a row for each (for no effect: 0). 
```
