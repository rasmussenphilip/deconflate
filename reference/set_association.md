# Set or replace one association in a model

Set or replace one association in a model

## Usage

``` r
set_association(model, disease1, disease2, value, measure = "OR")
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2:

  Disease ids.

- value:

  Association value.

- measure:

  Measure (see
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)).

## Value

The modified model.
