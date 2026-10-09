# Monte Carlo propagation of input uncertainty

Draws input sets with a sampler (usually from
[`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md))
and adjusts each with
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).
Every draw is a complete
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md),
so the uncertainty of probabilities, associations and impacts is
propagated jointly.

## Usage

``` r
cm_monte_carlo(
  sampler,
  n_draws,
  method = "simultaneous",
  seed = NULL,
  progress = FALSE,
  sampling = c("random", "lhs"),
  lhs_replicates = 10L,
  proposal = NULL,
  ...
)
```

## Arguments

- sampler:

  A function of the draw index returning a
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  (typically from
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)),
  or a
  [`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md).

- n_draws:

  Number of draws.

- method:

  Adjustment method(s): `"simultaneous"`, `"global"`, or `"published"`
  (the approximation of Rasmussen et al. 2022, eq. 16, kept here to
  compare with and reproduce earlier analyses; see
  [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)).

- seed:

  Optional random seed, for reproducibility.

- progress:

  Logical: print progress every 10% of draws?

- sampling:

  `"random"` (default) or `"lhs"` (needs a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)).

- lhs_replicates:

  Number of independent Latin hypercube blocks.

- proposal:

  Optional named list of `cm_dist` objects (importance sampling
  proposals), keyed as in `params` (e.g. `"impact:SCK"`). Needs a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md).

- ...:

  Passed to the adjustment (e.g. `joint`, `feasibility`, or arguments of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  for the global method).

## Value

A `cm_mc` object (or `cm_mc_batch` for a batch sampler) with:

- `draws`: raw and adjusted impacts and contributions by draw, method
  and disease;

- `totals`: naive and adjusted aggregate by draw and method;

- `params`: one row per accepted draw with the sampled inputs;

- `weights` (normalised), `log_weights`, `ess` (Kish effective sample
  size), `block` (LHS block of each accepted draw) and `n_blocks` (the
  number of LHS blocks sampled, including blocks with no accepted draw);

- `n_draws`, `n_rejected`, `rejections` (draw, type, reason),
  `sign_changes`, `method`, `sampling`, `proposal`, `specs`, `label` and
  `units`.

## Details

Draws that fail are rejected and counted by type (see `rejections`):
infeasible inputs (an impossible probability or association, or jointly
infeasible pairs), singular or non-identifiable systems, numerical
non-convergence, unsupported combinations, and non-finite results.
Conditioning on acceptance changes the effective input distribution, so
the rejection rate is part of the result and should be reported.

## Several methods

With more than one `method`, every draw is adjusted with each method, so
the methods are compared on identical inputs. A draw is rejected if any
method fails on it. Use
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
on the result.

## Batch runs

With a
[`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md),
each draw of the shared disease and association inputs is used by every
analysis, with the same draw identifiers; each analysis draws its own
impacts. A draw that fails in one analysis is rejected for that analysis
only. Batch runs use simple random sampling.

## Variance reduction

- `sampling = "lhs"`: Latin hypercube sampling, in `lhs_replicates`
  independent blocks. LHS draws are not independent, so the Monte Carlo
  standard error is estimated from the spread of the block means.

- `proposal`: importance sampling. Named inputs are drawn from proposal
  distributions, and each draw is weighted by the ratio of the input's
  own density to the proposal density (self-normalised). Each proposal
  must cover the support of the input's own distribution, which is
  checked from the distributions' supports; point masses (fixed values)
  cannot be importance-sampled. Self-normalised estimates have a small
  finite-sample bias. A defensive mixture (see
  [`cm_suggest_proposal()`](https://rasmussenphilip.github.io/deconflate/reference/cm_suggest_proposal.md))
  guarantees support but not a finite variance or better precision:
  compare the standard errors.

## Examples

``` r
s <- cm_sampler(example_supplement(), impacts = list(d1 = dist_normal(2.5, 0.5)))
mc <- cm_monte_carlo(s, 100, method = c("published", "simultaneous"), seed = 1)
compare_methods(mc)
#> <cm_comparison> methods: published, simultaneous
#> Units: %
#> Monte Carlo: 100 draws (0 rejected); statistic: mean
#> 
#> Adjusted values:
#>  disease raw_mean published simultaneous
#>       d1     2.53      2.09         2.17
#>       d2     5.00      3.70         3.39
#>       d3     7.50      6.75         6.93
#> 
#> Totals:
#>        quantity       method  mean  q0.5 trimmed_mean     mcse stability
#>  adjusted_total    published 2.114 2.109        2.113 0.003894        ok
#>  adjusted_total simultaneous 2.112 2.107        2.111 0.003940        ok
#>         raw_sum          raw 2.503 2.498        2.502 0.004294        ok
```
