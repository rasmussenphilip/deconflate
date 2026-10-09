# deconflate: comorbidity adjustment of disease impact estimates

Single-disease impact estimates (e.g. the milk-yield reduction in cows
with lameness compared with cows without) are conflated with the impacts
of associated diseases, which are more (or less) common among affected
animals. Summing such estimates double counts. `deconflate` adjusts the
estimates so that they can be aggregated.

## Workflow

1.  Read the inputs with
    [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):
    the diseases, their associations, one impact table and, optionally,
    interactions and three-way terms (CSV files with any names, or data
    frames; see
    [`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)).
    Values can have distributions in their own rows. The same inputs can
    be built in R with
    [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md),
    [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md),
    [`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md),
    [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md),
    [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
    and
    [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).
    For several outcomes, use one impact table each and adjust each in
    turn.

2.  Adjust the impacts with
    [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).
    Additive impacts (in any units) are adjusted with the exact
    simultaneous solution when every pair has an association and there
    are no interactions or three-way terms, and with the global
    (maximum-entropy) model otherwise; pairs without an association are
    unknown and filled in by the global model. Event impacts (hazard
    ratios, rate ratios, risk ratios, odds ratios or risk differences of
    e.g. culling) use the snapshot hazard model (`event_model = TRUE`,
    with the overall risk), which also gives the risk attributable to
    disease. When inputs have distributions,
    [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
    also reports intervals (`n_draws`).

3.  Inspect each disease's contribution
    ([`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md),
    [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md),
    [`plot_burden()`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)),
    in the units of the impacts. Converting the results into other
    quantities (e.g. a productivity gap or a monetary value) is left to
    the user.
    [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
    sets the methods side by side, including the proportional
    approximation of Rasmussen et al. (2022, eq. 16), which is kept for
    comparison and reproduction only.

4.  Check joint feasibility with
    [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)
    (a screen runs by default).

5.  Explore the inputs with
    [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md),
    [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md)
    (also without any association estimates),
    [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
    and
    [`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md),
    and find where a conclusion (a ranking, a sign or the total) changes
    as one input varies with
    [`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md).

[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
and
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
recompute the published tables. See
[`vignette("deconflate")`](https://rasmussenphilip.github.io/deconflate/articles/deconflate.md)
to get started.

## Conventions

- Additive impacts are in any units, the same within an impact table;
  results come back in those units. Event impacts are adjusted to hazard
  ratios, and risks are proportions of animals over the period.

- `E[k, i]` (see
  [`excess_matrix()`](https://rasmussenphilip.github.io/deconflate/reference/excess_matrix.md))
  is the excess probability of disease `k` among animals with disease
  `i`: `P(k | i) - P(k | not i)`.

- Under additive impacts, the crude impact of disease `i` satisfies
  `m_raw[i] = m[i] + sum_k E[k, i] * m[k]` exactly.

- Errors have classes `deconflate_infeasible`, `deconflate_singular`,
  `deconflate_nonconvergence`, `deconflate_unsupported` and
  `deconflate_nonfinite` (all also `deconflate_error`); problems found
  by
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  have class `deconflate_input_problems`. Warning classes include
  `deconflate_sign_change`, `deconflate_nonconvergence` and
  `deconflate_nonfinite_warning` (all also `deconflate_warning`).

## References

Rasmussen P, Shaw APM, Munoz V, Bruce M, Torgerson PR (2022). Estimating
the burden of multiple endemic diseases and health conditions using
Bayes' Theorem: a conditional probability model applied to UK dairy
cattle. Preventive Veterinary Medicine 203:105617.

Rasmussen P, et al. (2024). Global losses due to dairy cattle diseases:
a comorbidity-adjusted economic analysis. Journal of Dairy Science
107:6945-6970.

## See also

Useful links:

- <https://github.com/rasmussenphilip/deconflate>

- <https://rasmussenphilip.github.io/deconflate/>

- Report bugs at <https://github.com/rasmussenphilip/deconflate/issues>

## Author

**Maintainer**: Philip Rasmussen <phr@sund.ku.dk>

Authors:

- Philip Rasmussen <phr@sund.ku.dk>
