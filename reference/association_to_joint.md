# Joint probability of a disease pair from any supported measure

Joint probability of a disease pair from any supported measure

## Usage

``` r
association_to_joint(measure, value, p1, p2)
```

## Arguments

- measure:

  One of the measures in
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
  other than `"unknown"` (`"table"` rows store the table's odds ratio in
  `value`).

- value:

  Association value.

- p1, p2:

  Marginal probabilities of `disease1` and `disease2`.

## Value

`P(d1 and d2)`. Errors with class `deconflate_infeasible` if the measure
is incompatible with the marginals.
