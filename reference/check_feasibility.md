# Check whether the pairwise associations are jointly feasible

Each pairwise 2x2 table can be valid on its own while no distribution of
disease combinations satisfies all of them together. This function
checks joint feasibility before fitting the global model, and points to
the pairs that conflict.

## Usage

``` r
check_feasibility(
  model,
  method = c("auto", "triples", "lp"),
  max_lp_diseases = 14L,
  tol = 1e-09
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- method:

  `"auto"` (both checks if `lpSolve` is installed and the number of
  diseases is at most `max_lp_diseases`, otherwise triples only),
  `"triples"` or `"lp"`.

- max_lp_diseases:

  Largest number of diseases for the LP check.

- tol:

  Numerical tolerance.

## Value

A `cm_feasibility` object: `feasible` (`TRUE`, `FALSE`, or `NA` if only
the necessary triple check ran and passed), `triples` (data frame of
violating triples and their gaps), `lp` (`NULL` or a list with the
minimum total deviation and per-pair deviations) and `method`.

## Details

Two checks are available:

- `"triples"`: for every triple of diseases whose three pairs are all
  constrained, checks whether some probability of all three together
  makes every cell of the 2x2x2 table non-negative. This is exact for
  three diseases and a necessary condition in general. It needs no
  additional packages.

- `"lp"`: solves a linear programme over all 2^n combination
  probabilities that minimises the total absolute deviation from the
  pairwise joint probabilities, keeping marginals exact. A minimum of
  zero means the constraints are jointly feasible; otherwise the
  deviation needed for each pair identifies the conflicting
  associations. This is exact, and requires the `lpSolve` package.

## Examples

``` r
check_feasibility(example_supplement())
#> <cm_feasibility> jointly feasible [lp + triples]

# Three strongly linked diseases that cannot all be associated this way:
bad <- cm_model(
  cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
  cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05))
)
check_feasibility(bad, method = "triples")
#> <cm_feasibility> NOT jointly feasible [triples]
#> 
#> Conflicting triples:
#>  disease1 disease2 disease3   gap
#>         a        b        c 0.226
```
