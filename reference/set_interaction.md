# Set or replace one pairwise interaction in a model

Set or replace one pairwise interaction in a model

## Usage

``` r
set_interaction(model, disease1, disease2, value, outcome)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- disease1, disease2:

  Disease ids.

- value:

  Interaction value on the outcome's proportion scale.

- outcome:

  Outcome label.

## Value

The modified model.
