# Convert an odds ratio to a joint probability

Solves the quadratic of Rasmussen et al. (2022), eqs. 7-11, for the 2x2
table with marginal probabilities `p1` and `p2` and odds ratio `or`, and
returns `P(d1 and d2)`. A numerically stable form of the quadratic
formula is used, and the root is selected by feasibility rather than by
sign, so that odds ratios near or below 1 are handled correctly.

## Usage

``` r
or_to_joint(or, p1, p2)
```

## Arguments

- or:

  Odds ratio (positive).

- p1, p2:

  Marginal probabilities, strictly between 0 and 1.

## Value

The joint probability `P(d1 and d2)`.

## Examples

``` r
# Supplementary File, Rasmussen et al. (2022): P(1 | 2) = 0.163
or_to_joint(2, 0.10, 0.15) / 0.15
#> [1] 0.163196
```
