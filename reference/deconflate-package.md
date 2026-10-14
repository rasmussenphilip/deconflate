# deconflate: comorbidity adjustment of disease impact estimates

Single-disease impact estimates (e.g. the milk-yield reduction in cows
with lameness compared with cows without) are conflated with the impacts
of associated diseases, which are more (or less) common among affected
animals. Summing such estimates double counts. `deconflate` adjusts the
estimates so that they can be aggregated.

## Workflow

1.  Read the inputs with
    [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):
    the diseases, the association rows and one impact table (CSV files
    with any names, or data frames; see
    [`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)).
    An association row can be about a pair, a combination or a subgroup,
    crude or adjusted; an impact row can be about a disease or a
    combination (an interaction), within a subgroup, crude or adjusted;
    several rows can describe the same thing. Values can have
    distributions in their own rows. The same inputs can be built in R
    with
    [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md),
    [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md),
    [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
    and
    [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).
    For several outcomes, use one impact table each and adjust each in
    turn.

2.  Adjust the impacts with
    [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md).
    The package has one model: the maximum-entropy distribution of
    disease combinations that matches every association row (pairs that
    no row mentions are unknown and filled in), and the impacts and
    interactions that reproduce every impact row (exactly, or by
    weighted least squares when there are more rows than unknowns). When
    the pairwise tables give the same answer, they are used directly.
    Event impacts (hazard ratios, rate ratios, risk ratios, odds ratios
    or risk differences of e.g. culling) use the snapshot hazard model
    (`event_model = TRUE`, with the overall risk), which also gives the
    risk attributable to disease. When inputs have distributions,
    [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
    also reports intervals (`n_draws`).

3.  Inspect each disease's contribution
    ([`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md),
    [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md),
    [`plot_burden()`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)),
    in the units of the impacts, and how well the rows fit together
    (`$estimates`, `$heterogeneity`). Converting the results into other
    quantities (e.g. a productivity gap or a monetary value) is left to
    the user.
    [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
    sets the model beside the proportional approximation of Rasmussen et
    al. (2022, eq. 16), which is kept for comparison and reproduction
    only.

4.  Check joint feasibility with
    [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md).

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

- Errors are condition objects of class `deconflate_error`. Failures
  that code may want to handle separately have a subclass as well:
  `deconflate_infeasible` (inputs that no population can have),
  `deconflate_singular` (an estimate that is not identifiable, or rows
  that do not determine the impacts), `deconflate_nonconvergence` (a fit
  that did not converge), `deconflate_unsupported` (a combination of
  inputs and options that is not supported; with
  `deconflate_unknown_pairs` when a calculation needs a crude pair row
  for every pair), `deconflate_nonfinite` (non-finite results) and
  `deconflate_invalid_input` (an argument, table, row or input value
  that is not valid). Problems found by
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  have class `deconflate_input_problems` (also
  `deconflate_invalid_input`).

- Warnings are condition objects of class `deconflate_warning`, with the
  subclasses `deconflate_sign_change` (an adjusted impact on the other
  side of zero, or an adjusted hazard ratio on the other side of 1, from
  its raw estimate), `deconflate_nonfinite_warning` (non-finite adjusted
  impacts), `deconflate_nonconvergence` (a joint distribution that did
  not converge) and `deconflate_incomplete_allocation` (a Shapley
  allocation that skipped combinations; see
  [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)).

- Messages are of class `deconflate_message`; the notes of
  [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
  are also `deconflate_input_notes`.

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
