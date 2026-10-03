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
  ...
)
```

## Arguments

- sampler:

  A function of one argument (the draw index) returning a
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
  typically from
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

- n_draws:

  Number of draws.

- method:

  Adjustment method passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

- economics:

  Optional list with named elements `observed` and `unit_value` (and
  optionally `additional`), as in
  [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md).
  Elements may be numbers or `cm_dist` objects (drawn per draw). When
  given, productivity gaps and monetary losses are computed for every
  draw.

- seed:

  Optional random seed, for reproducibility.

- progress:

  Logical: print progress every 10% of draws?

- ...:

  Passed to
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).

## Value

A `cm_mc` object with elements:

- `draws`: long data frame of raw and adjusted impacts by draw;

- `losses`: long data frame of gaps and values by draw, outcome and
  disease (if `economics` was given);

- `params`: one row per accepted draw with the sampled inputs, keyed
  `prob:<id>`, `assoc:<d1>:<d2>`, `impact:<outcome>:<disease>`,
  `inter:<outcome>:<d1>:<d2>`, `observed:<outcome>`,
  `unit_value:<outcome>`;

- `weights` (initially equal), `n_draws`, `n_rejected`, `rejections`
  (draw and reason), `sign_changes` and `method`.

## Details

Draws whose inputs are infeasible (an impossible probability or odds
ratio, an association incompatible with the sampled marginals, or a
jointly infeasible set for the global method) are rejected and counted.
Conditioning on feasibility changes the effective input distribution, so
the rejection rate is part of the result and should be reported.
