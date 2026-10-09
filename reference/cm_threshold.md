# Threshold search: where does a conclusion change?

Varies one input over a range, with everything else fixed at the model's
values, and finds the input values at which a conclusion changes:

- `"rank"`: two diseases swap places, ranked by their contribution to
  the total (including interaction shares; `by = "adjusted"` ranks by
  adjusted impact instead);

- `"sign"`: a disease's adjusted impact crosses zero (for event impacts,
  its adjusted hazard ratio crosses 1). This is a sign change implied by
  the model and its inputs, not by itself evidence of a protective
  effect;

- `"total"`: the total crosses `target`;

- `"change"`: the total departs from its value at the model's own input
  by the relative amount `target` (e.g. `0.1` for +10%), which answers
  questions such as "what association strength would increase the total
  by 10%?".

## Usage

``` r
cm_threshold(
  model,
  input,
  range,
  conclusion = c("rank", "sign", "total", "change"),
  target = NULL,
  diseases = NULL,
  by = c("contribution", "adjusted"),
  method = "auto",
  event_model = FALSE,
  overall_risk = NULL,
  n_grid = 101L,
  log_scale = NULL,
  tol = 1e-08,
  ...
)
```

## Arguments

- model:

  A
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

- input:

  The input to vary, keyed as in
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  (`distributions`): `"assoc:<d1>:<d2>"` (the pair's association, on its
  measure; a pair without an association is varied as an odds ratio, so
  unknown pairs can be explored), `"inter:<d1>:<d2>"` (an interaction),
  `"prob:<disease>"` (the disease's `value`, on the scale it was
  entered), `"impact:<disease>"` (a raw impact),
  `"three:<d1>:<d2>:<d3>"` (a three-way ratio), or, for event impacts,
  `"risk"` (the overall risk).

- range:

  Two numbers: the range of the input to search.

- conclusion:

  `"rank"`, `"sign"`, `"total"` or `"change"`.

- target:

  For `"total"`, the value of the aggregate (for event impacts, of the
  attributable risk); for `"change"`, the relative change from the
  baseline (e.g. `0.1`, or `-0.1`).

- diseases:

  Optional disease ids to restrict `"rank"` (pairs among them) and
  `"sign"`.

- by:

  For `"rank"`: `"contribution"` (default) or `"adjusted"`.

- method, event_model, overall_risk:

  As in
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  (a distribution of the overall risk is used at its mean).

- n_grid:

  Number of grid points.

- log_scale:

  Use a log-spaced grid? Default `TRUE` for associations measured as
  ratios and three-way ratios with a positive range.

- tol:

  Relative tolerance on the threshold. (The bracket is also narrowed to
  at most 1e-6 of the grid spacing, so that the continuity test is
  reliable.) Because of this argument, the `tol` of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  cannot be passed through `...`.

- ...:

  `joint`, or arguments of
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  (e.g. `backend = "sampled"`).

## Value

A `cm_threshold` object: `thresholds` (one row per crossing: conclusion,
item, status, threshold, bracket `lower`/`upper`, the compared quantity
just below and above, and a description), `details` (for each crossing,
the adjusted results at the two ends of its bracket), `scan` (each grid
point with its status, the aggregate and the condition number),
`regions` (stretches of unusable grid points), `summary` (per item:
number of thresholds, and whether part of the range could not be
evaluated), and the settings (`input`, `range`, `baseline` value of the
input, `method`, `conclusion`, `target`, `fixed`).

## Details

The total is the adjusted aggregate (in the units of the impacts) for
additive impacts, and the risk attributable to disease for event impacts
(`event_model = TRUE`), whose contributions are its Shapley allocation.

## Method

The response need not be monotonic, so the range is first scanned on a
grid. Every grid point is adjusted and classified: `"ok"`, or the reason
it could not be used (`"infeasible"` inputs, a `"singular"` system, an
`"undefined"` (non-finite) result, an `"unresolved"` joint fit, an
`"unsupported"` combination or another `"error"`); an input value that
is itself invalid (e.g. a non-positive odds ratio or a probability of 1)
is `"infeasible"`. Each pair of neighbouring usable grid points whose
conclusion differs brackets a crossing, which is refined by bisection. A
grid point where the compared quantity is exactly zero (or a stretch of
such points) is a crossing only if the quantity has opposite signs at
the usable points on either side; the threshold is then the first zero
point. A refined crossing is reported as a `"threshold"` only if the
quantity being compared is close to zero on both sides of the final
bracket (continuity); a jump across a pole (e.g. where the system
becomes singular) is reported as a `"discontinuity"`, never as a
threshold; so is a bracket whose bisection reaches a singular point. If
a bisection step lands on any other unusable point, the crossing is
`"unresolved"`. All crossings in the range are reported; stretches of
unusable grid points are listed in `regions`, so that "no crossing
found" can be told apart from "part of the range could not be
evaluated". A sign change between usable points separated by unusable
ones is not reported as a crossing; such stretches appear in `regions`.

Each point is adjusted as
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
would (point estimates; the method is chosen for the point's own
inputs). When the input is an impact, an interaction or the overall risk
(which do not change the joint distribution) and the adjustment uses the
joint distribution, it is fitted once and reused (or the `joint` passed
through `...` is used).

Thresholds are deterministic: they are for the model's input values.

## Examples

``` r
m <- example_supplement()
# At what odds ratio between d1 and d2 does the ranking change?
th <- cm_threshold(m, "assoc:d1:d2", c(0.2, 20), conclusion = "rank")
th
#> <cm_threshold> rank over assoc:d1:d2 in [0.2, 20] (baseline 2; method simultaneous)
#> Fixed: all inputs other than assoc:d1:d2 at the model's values.
#> 
#> No crossing found in the range.
# What odds ratio between d1 and d3 (an odds ratio of 1 in the example)
# would reduce the aggregate by 10%?
cm_threshold(m, "assoc:d1:d3", c(1, 50), conclusion = "change", target = -0.1)$thresholds
#>   conclusion  item    status threshold    lower    upper        below
#> 1     change total threshold  5.263907 5.263907 5.263907 4.931178e-10
#>           above                                  description
#> 1 -3.748852e-10 the aggregate falls below the baseline - 10%
```
