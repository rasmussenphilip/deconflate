# Set or replace one pairwise interaction

Set or replace one pairwise interaction

## Usage

``` r
set_interaction(model, disease1, disease2, value)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2:

  Disease ids.

- value:

  Interaction value (same units as the impacts).

## Value

The modified model. A distribution of the pair's interaction in the
model's `distributions` is dropped (the scenario value is fixed).
