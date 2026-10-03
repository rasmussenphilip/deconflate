# Getting started with deconflate

``` r

library(deconflate)
```

## The problem

Single-disease impact estimates compare animals with and without a
disease. Animals with one disease are more (or less) likely to have
associated diseases, so each estimate partly reflects the impacts of
those other diseases. Adding such estimates double counts. `deconflate`
adjusts (“de-conflates”) the estimates so that they can be aggregated.

## Describing a system

A model combines disease probabilities, pairwise associations and raw
impacts. This is the three-disease example from the Supplementary File
of Rasmussen et al. (2022):

``` r

m <- cm_model(
  diseases = cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
  associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3), measure = "OR"),
  impacts = cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5),
                       outcome = "yield", scale = "percent")
)
m
#> <cm_model>
#>   Diseases: 3 (d1, d2, d3)
#>   Disease pairs: 3 [independent (default): 1; specified: 2]
#>   Outcomes: yield
```

Pairs without an association are independent by default
(`missing_associations = "independent"`). Use `"unknown"` to leave them
unconstrained instead; this is a different assumption (see
[`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md)).

## Adjusting impacts

``` r

compare_methods(m)
#>   outcome disease   raw  published simultaneous     global
#> 1   yield      d1 0.025 0.02065001   0.02143250 0.02143250
#> 2   yield      d2 0.050 0.03699294   0.03387080 0.03387080
#> 3   yield      d3 0.075 0.06748473   0.06934209 0.06934209
```

- `"simultaneous"` (the default) solves the additive impact equations
  exactly.
- `"published"` is the proportional approximation of Rasmussen et al.
  (2022), kept for reproduction and comparison.
- `"global"` fits the distribution of disease combinations and can
  include impact interactions.

``` r

res <- deconflate(m)
res
#> <cm_result> method: simultaneous
#> 
#>  outcome disease   raw adjusted   change
#>    yield      d1 0.025  0.02143 -0.14270
#>    yield      d2 0.050  0.03387 -0.32258
#>    yield      d3 0.075  0.06934 -0.07544
#> 
#> Diagnostics:
#>  outcome max_reconstruction_residual n_sign_changes condition_number
#>    yield                           0              0             1.53
```

The diagnostics report how well the adjusted impacts reconstruct the raw
ones, and flag adjusted impacts whose sign differs from the raw impact.
A sign change means that the raw impact is smaller than the associated
diseases alone would produce, so the inputs are inconsistent under the
additive model.

## Productivity gaps and losses

``` r

gaps <- productivity_gap(res, c(yield = 10000))
gaps$summary
#>   outcome observed disease_free      gap total_loss_fraction
#> 1   yield    10000     10215.47 215.4676          0.02109229
value_losses(gaps, unit_value = c(yield = 0.30))$by_disease
#>   outcome disease       gap     value
#> 1   yield      d1  21.89431  6.568292
#> 2   yield      d2  51.90090 15.570271
#> 3   yield      d3 141.67238 42.501715
```

## Checking the method with a simulated herd

[`simulate_raw_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/simulate_raw_impacts.md)
computes the raw impacts that single-disease studies would report in a
herd with known true impacts. Adjusting them should recover the truth:

``` r

truth <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)
m_sim <- m
m_sim$impacts <- simulate_raw_impacts(m, truth, outcome = "yield")
compare_methods(m_sim, methods = c("published", "simultaneous"))
#>   outcome disease        raw  published simultaneous
#> 1   yield      d1 0.02421306 0.01960287         0.02
#> 2   yield      d2 0.05406438 0.04185061         0.04
#> 3   yield      d3 0.06668175 0.05872788         0.06
```

## Where next

- [`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md):
  the 2022 and 2024 analyses.
- [`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md):
  the joint distribution, unknown associations, impact interactions and
  attribution.
- [`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md):
  Monte Carlo analysis, correlated outcomes, scenarios and sensitivity
  screening.
