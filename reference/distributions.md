# Input distributions

Constructors for the distributions of uncertain inputs, drawn by
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
(`n_draws`): for the `distributions` of
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
and for `overall_risk`. (The input tables give the same distributions in
their `dist` and `p1`-`p4` columns; see
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md).)
Each returns a `cm_dist` object with a quantile function (used for
sampling, including Latin hypercube sampling), a distribution function
and a log-density.

## Usage

``` r
dist_fixed(value)

dist_normal(mean, sd, lower = -Inf, upper = Inf)

dist_lognormal(meanlog, sdlog)

dist_lognormal_ci(estimate, ci_lower, ci_upper, level = 0.95)

dist_beta(shape1, shape2, min = 0, max = 1)

dist_pert(min, mode, max, lambda = 4)

dist_pert_mean(min, mean, max, lambda = 4)

dist_uniform(min, max)

dist_mixture(..., weights = NULL)
```

## Arguments

- value:

  Constant value.

- mean, sd:

  Normal mean and standard deviation.

- lower, upper:

  Truncation bounds.

- meanlog, sdlog:

  Log-normal parameters (log scale).

- estimate, ci_lower, ci_upper, level:

  Point estimate and confidence interval for `dist_lognormal_ci()`.

- shape1, shape2:

  Beta shape parameters.

- min, max:

  Range.

- mode:

  Mode of the PERT distribution.

- lambda:

  PERT shape parameter (4 for the standard PERT).

- ...:

  `cm_dist` components for `dist_mixture()`.

- weights:

  Mixture weights (normalised internally).

## Value

A `cm_dist` object, with the quantile function `q`, distribution
function `p`, log-density `logd`, sampler `r`, `mean`, and its support:
`support` (a two-column matrix of disjoint intervals on which the
density is positive), its range `lower` and `upper`, and
`discrete = TRUE` for point masses. The support of a mixture is the
union of the supports of its components with positive weight, so it can
have gaps.

## Details

- `dist_fixed()`: a constant.

- `dist_normal()`: normal, optionally truncated to `[lower, upper]`
  (e.g. `lower = 0` for odds ratios, which conditions on positive
  values).

- `dist_lognormal()`: log-normal; `dist_lognormal_ci()` builds one from
  a point estimate and confidence interval (e.g. a published odds
  ratio).

- `dist_beta()`: beta, optionally rescaled to `[min, max]`.

- `dist_pert()`: PERT with `min`, `mode` and `max`. `dist_pert_mean()`
  builds one from a mean instead of a mode. Note: in Rasmussen et al.
  (2024), Tables 2-4, the central value of PERT distributions is the
  mode.

- `dist_uniform()`: uniform.

- `dist_mixture()`: a finite mixture of distributions (e.g. two pooled
  sources of evidence).

## Examples

``` r
d <- dist_pert(1.19, 3.30, 10.71)
d$mean
#> [1] 4.183333
d$r(5)
#> [1] 4.763111 2.570940 7.339168 6.517981 7.144411
dist_lognormal_ci(2.7, 1.5, 4.9)
#> <cm_dist> lognormal(meanlog = 0.9933, sdlog = 0.302), mean 2.826
```
