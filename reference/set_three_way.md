# Set or replace one three-way association term

Set or replace one three-way association term

## Usage

``` r
set_three_way(model, disease1, disease2, disease3, ratio)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2, disease3:

  Disease ids.

- ratio:

  Ratio of conditional odds ratios (see
  [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md)).

## Value

The modified object. A distribution of the triple's ratio in the model's
`distributions` is dropped (the scenario value is fixed).
