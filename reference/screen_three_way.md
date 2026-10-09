# Screen three-way association scenarios

Pairwise associations do not identify how often three diseases occur
together; the global model assumes no three-way association. This screen
sets a three-way term (see
[`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md))
for each triple in turn, refits the joint distribution (all pairwise
associations are still matched) and reports the change in the aggregate
and rankings.

## Usage

``` r
screen_three_way(
  model,
  ratios = c(0.5, 2),
  triples = NULL,
  event_model = FALSE,
  overall_risk = NULL
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- ratios:

  Ratios of conditional odds ratios to try.

- triples:

  Optional list of character vectors of three disease ids; default all
  triples (which can be many).

- event_model, overall_risk:

  As in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Value

A `cm_screen` data frame as in
[`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).

## Details

When the pairs of a triple are all constrained, a three-way term keeps
their tables fixed, so it changes additive results only through
interactions: without interactions such a scenario reproduces the
baseline. When a pair is unknown (unconstrained), the fitted pairwise
table changes with the three-way term, and additive results can change
even without interactions. Three-way terms also matter for event impacts
(`event_model = TRUE`).
