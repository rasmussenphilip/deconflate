# Set or replace one association

Set or replace one association

## Usage

``` r
set_association(model, disease1, disease2, value, measure = "OR")
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2:

  Disease ids.

- value:

  Association value.

- measure:

  Measure (see
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)).

## Value

The modified object.
