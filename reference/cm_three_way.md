# Describe three-way disease associations (scenarios)

Pairwise associations do not determine how often three diseases occur
together. The global model
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md))
fills this in by maximum entropy, which assumes no three-way
association: the odds ratio of a pair is the same whether or not a third
disease is present (given the other diseases). A three-way term relaxes
this for a chosen triple, while all pairwise associations are still
matched.

## Usage

``` r
cm_three_way(disease1, disease2, disease3, ratio, source = NA_character_)
```

## Arguments

- disease1, disease2, disease3:

  Character vectors of disease ids.

- ratio:

  Positive ratio of conditional odds ratios (see Details).

- source:

  Optional citation or scenario label.

## Value

A `cm_three_way` data frame.

## Details

`ratio` is the ratio of conditional odds ratios,
`OR(d1, d2 | d3 present) / OR(d1, d2 | d3 absent)` (given the other
diseases), which is symmetric in the three diseases. `ratio = 1` is the
maximum-entropy assumption. Pairwise evidence cannot identify `ratio`,
so three-way terms are sensitivity scenarios unless there is direct
evidence.

When every pair in the triple has an association, a three-way term
leaves the pairwise tables unchanged, so additive results without
interactions do not change; it affects only results that depend on the
joint distribution: additive impacts with interactions, and event
impacts (`event_model = TRUE` in
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
When a pair is unknown (no association given), the global fit determines
its table, and a three-way term can change it, and with it the additive
results. A model with three-way terms is always adjusted with the global
method.

## Examples

``` r
cm_three_way("d1", "d2", "d3", ratio = 2)
#>   disease1 disease2 disease3 ratio source
#> 1       d1       d2       d3     2   <NA>
```
