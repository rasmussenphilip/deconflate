# Simulate the raw impacts that studies would report

For validation: given a population whose disease combinations follow the
maximum-entropy distribution of `model` (see
[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)),
and known true additive impacts (and optional interactions), computes
what a study would estimate for each disease: the crude difference
between animals with and without the disease, or the coefficient of an
additive regression adjusted for other diseases. With `n = NULL` the
exact population values are returned; with `n` set, `n` animals are
sampled, which adds sampling error.

## Usage

``` r
simulate_raw_impacts(
  model,
  true_impacts,
  interactions = NULL,
  estimand = "crude",
  adjusted_for = NA_character_,
  n = NULL,
  joint = NULL,
  units = NULL
)
```

## Arguments

- model:

  A
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  or
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  (its impacts are ignored).

- true_impacts:

  Named numeric vector of true impacts, one per disease (any units).

- interactions:

  Optional
  [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  with true interaction values (same units).

- estimand:

  `"crude"` or `"adjusted_linear"`, one per disease or recycled.

- adjusted_for:

  For `"adjusted_linear"`: adjustment sets (ids separated by `";"`, or
  `"all"`), one per disease or recycled.

- n:

  Optional number of animals to sample.

- joint:

  Optional pre-computed
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result. It must match the model and have converged; a joint fitted
  here that does not converge stops with class
  `deconflate_nonconvergence`.

- units:

  Optional units label for the returned impacts.

## Value

A
[`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
object with the simulated raw impacts.

## Details

Adjusting the returned raw impacts with
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
should recover `true_impacts` exactly (for `n = NULL`): with the
simultaneous or global method without interactions, and with the global
method with them.

## Examples

``` r
raw <- simulate_raw_impacts(example_supplement(), c(d1 = 2, d2 = 4, d3 = 6))
raw$value
#> [1] 2.421306 5.406438 6.668175
deconflate(cm_model(example_supplement(), raw))$adjusted$adjusted
#> [1] 2 4 6
```
