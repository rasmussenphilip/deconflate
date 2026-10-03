# Simulate the raw impacts a single-disease study would report

For validation: given a population whose disease combinations follow the
maximum-entropy distribution implied by `model` (see
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)),
and known true additive impacts (and optional interactions), computes
the crude difference in outcome between animals with and without each
disease. With `n = NULL` the exact expectation is returned. With `n`
set, `n` animals are sampled and the empirical crude differences are
returned, which adds sampling error.

## Usage

``` r
simulate_raw_impacts(
  model,
  true_impacts,
  interactions = NULL,
  outcome = "impact",
  n = NULL,
  joint = NULL
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  (its impacts, if any, are ignored).

- true_impacts:

  Named numeric vector of true proportional impacts, one per disease.

- interactions:

  Optional
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  with true interaction values (only rows for `outcome` are used).

- outcome:

  Outcome label for the returned impacts.

- n:

  Optional number of animals to sample.

- joint:

  Optional pre-computed
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result.

## Value

A
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
object with the simulated raw impacts.

## Details

Adjusting the returned raw impacts with
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
should recover `true_impacts` (exactly for `"simultaneous"`/`"global"`
without interactions, and for `"global"` with interactions, when
`n = NULL`).

## Examples

``` r
m <- example_supplement()
raw <- simulate_raw_impacts(m, c(d1 = 0.02, d2 = 0.04, d3 = 0.06))
raw$value
#> [1] 0.02421306 0.05406438 0.06668175
```
