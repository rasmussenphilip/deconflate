# Describe pairwise impact interactions

An interaction `delta` is the additional impact when both diseases are
present, on the same scale as the outcome's impacts: positive values are
synergistic (more loss than the sum), negative values antagonistic.
Interactions cannot be inferred from associations and must come from
evidence or explicit scenarios. They require `method = "global"` in
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Usage

``` r
cm_interactions(
  disease1,
  disease2,
  value,
  outcome = "impact",
  scale = "proportion",
  source = NA_character_
)
```

## Arguments

- disease1, disease2:

  Character vectors of disease ids.

- value:

  Numeric interaction values.

- outcome:

  Outcome label(s) matching
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md).

- scale:

  `"proportion"` or `"percent"`; additive on that scale.

- source:

  Optional citation or scenario label.

## Value

A `cm_interactions` data frame.
