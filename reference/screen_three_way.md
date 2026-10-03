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
screen_three_way(model, ratios = c(0.5, 2), triples = NULL, valuation = NULL)
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

- valuation:

  Optional valuation list.

## Value

A `cm_screen` data frame as in
[`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md).

## Details

Three-way terms change additive results only through interactions: for a
model without interactions the global result depends on the pairs alone,
and every scenario reproduces the baseline. They matter for
interactions, and for hazard ratios
([`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md))
and
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md).
