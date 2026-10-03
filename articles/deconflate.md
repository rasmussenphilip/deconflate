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
#> <cm_comparison> methods: published, simultaneous, global
#> 
#> Adjusted impacts:
#>  outcome disease unit raw published simultaneous global
#>    yield      d1    % 2.5      2.07         2.14   2.14
#>    yield      d2    % 5.0      3.70         3.39   3.39
#>    yield      d3    % 7.5      6.75         6.93   6.93
#> 
#> Totals:
#>  outcome       method raw_loss adjusted_loss
#>    yield    published    0.025       0.02111
#>    yield simultaneous    0.025       0.02109
#>    yield       global    0.025       0.02109
```

- `"simultaneous"` (the default) solves the system of impact equations
  exactly.
- `"published"` is the simple proportional approximation used in
  Rasmussen et al. (2022, 2024), kept for reproduction and comparison.
- `"global"` is the iterative model: it fits the distribution of disease
  combinations and can include impact interactions.

[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
runs all three and tabulates them side by side; with `economics` it also
compares the gaps and their values.

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
#> <cm_comparison> methods: published, simultaneous
#> 
#> Adjusted impacts:
#>  outcome disease unit  raw published simultaneous
#>    yield      d1    % 2.42      1.96            2
#>    yield      d2    % 5.41      4.19            4
#>    yield      d3    % 6.67      5.87            6
#> 
#> Totals:
#>  outcome       method raw_loss adjusted_loss
#>    yield    published  0.02387       0.01998
#>    yield simultaneous  0.02387       0.02000
```

## Where next

- [`vignette("own-data")`](https://rasmussenphilip.github.io/deconflate/articles/own-data.md):
  reading your own inputs from CSV files.
- [`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md):
  culling hazard ratios and the culling attributable to disease.
- [`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md):
  the 2022 and 2024 analyses.
- [`vignette("interactions")`](https://rasmussenphilip.github.io/deconflate/articles/interactions.md):
  the joint distribution, unknown associations, impact interactions and
  attribution.
- [`vignette("uncertainty")`](https://rasmussenphilip.github.io/deconflate/articles/uncertainty.md):
  Monte Carlo analysis, method comparison, stability checks, correlated
  outcomes, scenarios and sensitivity screening.
