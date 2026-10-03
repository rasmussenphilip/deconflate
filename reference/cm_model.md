# Combine inputs into a comorbidity model

Combine inputs into a comorbidity model

## Usage

``` r
cm_model(
  diseases,
  associations = NULL,
  impacts = NULL,
  interactions = NULL,
  missing_associations = c("independent", "unknown")
)
```

## Arguments

- diseases:

  A
  [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  object.

- associations:

  A
  [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
  object, or `NULL`.

- impacts:

  A
  [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  object, or `NULL` (e.g. when only the joint distribution or simulated
  impacts are needed).

- interactions:

  A
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  object, or `NULL`.

- missing_associations:

  How to treat disease pairs without a row in `associations`:
  `"independent"` imposes an odds ratio of 1 (as in Rasmussen et al.
  2022, 2024); `"unknown"` leaves them unconstrained (only usable with
  the global method). These are different assumptions.

## Value

A `cm_model` object.
