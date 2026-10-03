# Probabilities of disease combinations

Probabilities of disease combinations

## Usage

``` r
combination_probs(joint, min_prob = 0)
```

## Arguments

- joint:

  A
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result.

- min_prob:

  Drop combinations with probability below this value.

## Value

A data frame with one row per combination: a 0/1 column per disease, the
number of diseases present and the probability, sorted by decreasing
probability.
