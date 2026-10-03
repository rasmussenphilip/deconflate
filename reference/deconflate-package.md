# deconflate: comorbidity adjustment of disease impact estimates

Single-disease impact estimates (e.g. the milk-yield reduction in cows
with lameness compared with cows without) are conflated with the impacts
of associated diseases, which are more (or less) common among affected
animals. Summing such estimates double counts. `deconflate` adjusts the
estimates so that they can be aggregated.

## Workflow

1.  Describe the system with
    [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md),
    [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md),
    [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
    and, optionally,
    [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md),
    then combine them with
    [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).
    Or read the same tables from CSV files or data frames with
    [`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
    (see
    [`cm_template()`](https://rasmussenphilip.github.io/deconflate/reference/cm_template.md)).

2.  Adjust impacts with
    [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
    using one of three methods:

    - `"simultaneous"` (default): the exact solution of the additive
      impact equations, built from pairwise 2x2 tables.

    - `"published"`: the proportional approximation of Rasmussen et al.
      (2022), eqs. 15-16, for reproduction and comparison.

    - `"global"`: fits a maximum-entropy distribution of disease
      combinations with
      [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md),
      then solves the impact equations including pairwise interactions.

3.  Estimate, value and attribute productivity gaps with
    [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md),
    [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md),
    [`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md)
    and
    [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md).
    Culling hazard ratios can be adjusted directly
    (`scale = "hazard_ratio"`) and turned into culling attributable to
    disease with
    [`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md);
    [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
    tabulates the methods side by side.

4.  Check joint feasibility with
    [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md).

5.  Propagate uncertainty with
    [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
    and
    [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
    (check stability with
    [`cm_diagnose()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diagnose.md)),
    and explore scenarios with
    [`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md),
    [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md),
    [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md),
    [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
    and
    [`compare_scenarios()`](https://rasmussenphilip.github.io/deconflate/reference/compare_scenarios.md).

See
[`vignette("deconflate")`](https://rasmussenphilip.github.io/deconflate/articles/deconflate.md)
to get started.

## Conventions

- Impacts on the `"proportion"` scale are proportional changes relative
  to the disease-free value (e.g. `0.025` for a 2.5% yield reduction).

- `E[k, i]` (see
  [`excess_matrix()`](https://rasmussenphilip.github.io/deconflate/reference/excess_matrix.md))
  is the excess probability of disease `k` among animals with disease
  `i`: `P(k | i) - P(k | not i)`.

- Under additive impacts, the raw (crude) impact of disease `i`
  satisfies `m_raw[i] = m[i] + sum_k E[k, i] * m[k]` exactly.

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
