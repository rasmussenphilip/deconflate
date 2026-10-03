# Monte Carlo propagation of input uncertainty

Draws input sets with a sampler (usually from
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md))
and adjusts each with
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).
Because every draw is a complete
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
disease probabilities and associations are shared across all outcomes
within a draw, so the outcomes' uncertainty stays dependent.

## Usage

``` r
cm_monte_carlo(
  sampler,
  n_draws,
  method = "simultaneous",
  economics = NULL,
  seed = NULL,
  progress = FALSE,
  sampling = c("random", "lhs"),
  proposal = NULL,
  ...
)
```

## Arguments

- sampler:

  A function of the draw index returning a
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  typically from
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

- n_draws:

  Number of draws.

- method:

  Adjustment method(s) passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md):
  one or more of `"simultaneous"`, `"published"` and `"global"`.

- economics:

  Optional list with named elements `observed` and `unit_value` (and
  optionally `additional`), as in
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md).
  Elements may be numbers or `cm_dist` objects (drawn per draw). When
  given, productivity gaps and monetary losses are computed for every
  draw. Hazard-ratio outcomes cannot be valued here; use
  [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
  on individual results.

- seed:

  Optional random seed, for reproducibility.

- progress:

  Logical: print progress every 10% of draws?

- sampling:

  `"random"` (default) or `"lhs"` (Latin hypercube; needs a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)).

- proposal:

  Optional named list of `cm_dist` objects (importance sampling
  proposals), keyed as in `params` (e.g. `"impact:fertility:SCK"`).
  Needs a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

- ...:

  Passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Value

A `cm_mc` object with elements:

- `draws`: long data frame of raw and adjusted impacts by draw and
  method;

- `losses`: long data frame of gaps and values by draw, method, outcome
  and disease (if `economics` was given);

- `params`: one row per accepted draw with the sampled inputs, keyed
  `prob:<id>`, `assoc:<d1>:<d2>`, `impact:<outcome>:<disease>`,
  `inter:<outcome>:<d1>:<d2>`, `observed:<outcome>`,
  `unit_value:<outcome>`;

- `weights` (equal, or importance weights), `log_weights`, `ess`,
  `n_draws`, `n_rejected`, `rejections` (draw and reason),
  `sign_changes`, `method`, `sampling` and `proposal`.

## Details

Draws whose inputs are infeasible (an impossible probability or odds
ratio, an association incompatible with the sampled marginals, or a
jointly infeasible set for the global method) are rejected and counted.
Conditioning on feasibility changes the effective input distribution, so
the rejection rate is part of the result and should be reported.

## Several methods

With more than one `method`, every accepted draw is adjusted with each
method, so the methods are compared on identical inputs. A draw is
rejected if any method fails on it. Use
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
on the result for a side-by-side table.

## Stabilising the estimates

[`summary.cm_mc()`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md)
flags unstable estimates and suggests remedies (see
[`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)).
Two variance-reduction options are available here:

- `sampling = "lhs"`: Latin hypercube sampling of the inputs that are
  not drawn with an outcome correlation. It stratifies each input's
  distribution and usually reduces the Monte Carlo error of means.

- `proposal`: importance sampling. Named inputs are drawn from the given
  proposal distributions instead of their own, and each draw is weighted
  by the ratio of the input densities. A proposal must cover the whole
  range of the input's own distribution; a defensive mixture such as the
  one built by
  [`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md)
  does. Importance sampling reduces the error of means that exist; it
  cannot fix a mean that does not exist (see
  [`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)).

## Examples

``` r
s <- cm_sampler(example_supplement(),
                impacts = list("yield:d1" = dist_normal(2.5, 0.5)))
mc <- cm_monte_carlo(s, 100, method = c("published", "simultaneous"), seed = 1)
compare_methods(mc)
#> <cm_comparison> methods: published, simultaneous
#> Monte Carlo: 100 draws (0 rejected); statistic: mean
#> 
#> Adjusted impacts:
#>  outcome disease unit raw_mean published simultaneous
#>    yield      d1    %     2.53      2.09         2.17
#>    yield      d2    %     5.00      3.70         3.39
#>    yield      d3    %     7.50      6.75         6.93
```
