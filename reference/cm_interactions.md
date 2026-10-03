# Describe pairwise impact interactions

An interaction `delta` is the additional impact when both diseases are
present, in the same units as the impact vector: positive values are
synergistic (more than the sum), negative values antagonistic.
Interactions cannot be inferred from associations and must come from
evidence or explicit scenarios. They require `method = "global"` in
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Usage

``` r
cm_interactions(disease1, disease2, value, source = NA_character_)
```

## Arguments

- disease1, disease2:

  Character vectors of disease ids.

- value:

  Numeric interaction values (same units as the impacts).

- source:

  Optional citation or scenario label.

## Value

A `cm_interactions` data frame.

## Examples

``` r
cm_interactions("d1", "d2", 0.5)
#>   disease1 disease2 value source
#> 1       d1       d2   0.5   <NA>
```
