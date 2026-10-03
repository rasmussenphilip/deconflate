# Package index

## Describing a system

- [`cm_diseases()`](https://rasmussenphilip.github.io/deconflate/reference/cm_diseases.md)
  : Describe the diseases in a system
- [`cm_associations()`](https://rasmussenphilip.github.io/deconflate/reference/cm_associations.md)
  : Describe statistical associations between disease pairs
- [`cm_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/cm_impacts.md)
  : Describe raw (unadjusted) disease impact estimates
- [`cm_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/cm_interactions.md)
  : Describe pairwise impact interactions
- [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  : Combine inputs into a comorbidity model
- [`combine_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/combine_impacts.md)
  : Combine impact tables
- [`set_association()`](https://rasmussenphilip.github.io/deconflate/reference/set_association.md)
  : Set or replace one association in a model
- [`set_interaction()`](https://rasmussenphilip.github.io/deconflate/reference/set_interaction.md)
  : Set or replace one pairwise interaction in a model

## Adjustment

- [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  : Adjust raw impact estimates for comorbidity
- [`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
  : Compare adjustment methods
- [`pair_tables()`](https://rasmussenphilip.github.io/deconflate/reference/pair_tables.md)
  : Pairwise 2x2 tables for every disease pair
- [`excess_matrix()`](https://rasmussenphilip.github.io/deconflate/reference/excess_matrix.md)
  : Excess-probability matrix from pairwise tables
- [`or_to_joint()`](https://rasmussenphilip.github.io/deconflate/reference/or_to_joint.md)
  : Convert an odds ratio to a joint probability
- [`joint_to_or()`](https://rasmussenphilip.github.io/deconflate/reference/joint_to_or.md)
  : Odds ratio implied by a joint probability
- [`association_to_joint()`](https://rasmussenphilip.github.io/deconflate/reference/association_to_joint.md)
  : Joint probability of a disease pair from any supported measure

## Joint distribution and feasibility

- [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  : Fit the maximum-entropy distribution of disease combinations
- [`combination_probs()`](https://rasmussenphilip.github.io/deconflate/reference/combination_probs.md)
  : Probabilities of disease combinations
- [`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)
  : Check whether the pairwise associations are jointly feasible

## Gaps, losses and attribution

- [`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md)
  : Productivity gaps and their attribution
- [`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md)
  : Value productivity gaps in monetary terms
- [`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md)
  : Attribute the aggregate burden to diseases (Shapley allocation)
- [`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
  : Shapley attribution of a general loss function over disease
  combinations
- [`loss_additive()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md)
  [`loss_multiplicative()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md)
  : Loss functions for Shapley attribution
- [`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)
  : Contribution table
- [`summary(`*`<cm_result>`*`)`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_result.md)
  : Summarise an adjustment result

## Culling and hazard ratios

- [`hr_to_risk()`](https://rasmussenphilip.github.io/deconflate/reference/hr_to_risk.md)
  : Convert culling (or mortality) hazard ratios to period risks
- [`excess_to_hr()`](https://rasmussenphilip.github.io/deconflate/reference/excess_to_hr.md)
  : Convert an adjusted excess risk back to a hazard ratio
- [`hr_conversion()`](https://rasmussenphilip.github.io/deconflate/reference/hr_conversion.md)
  : Convert culling hazard ratios into excess-risk impacts
- [`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md)
  : Turn a hazard-ratio conversion into impacts
- [`adjusted_hr()`](https://rasmussenphilip.github.io/deconflate/reference/adjusted_hr.md)
  : Convert adjusted culling impacts back to hazard ratios

## Uncertainty and scenarios

- [`dist_fixed()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_normal()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_lognormal()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_lognormal_ci()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_beta()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_pert()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_pert_mean()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_uniform()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  [`dist_mixture()`](https://rasmussenphilip.github.io/deconflate/reference/distributions.md)
  : Input distributions for Monte Carlo analysis
- [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md)
  : Build a Monte Carlo sampler from input distributions
- [`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
  : Monte Carlo propagation of input uncertainty
- [`summary(`*`<cm_mc>`*`)`](https://rasmussenphilip.github.io/deconflate/reference/summary.cm_mc.md)
  : Summarise Monte Carlo results
- [`cm_reweight()`](https://rasmussenphilip.github.io/deconflate/reference/cm_reweight.md)
  : Reweight Monte Carlo draws (importance sampling)
- [`cm_scenario()`](https://rasmussenphilip.github.io/deconflate/reference/cm_scenario.md)
  : Scenario analysis by reweighting Monte Carlo draws

## Sensitivity

- [`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md)
  : One-at-a-time sensitivity analysis
- [`screen_associations()`](https://rasmussenphilip.github.io/deconflate/reference/screen_associations.md)
  : Screen disease pairs for influential associations
- [`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
  : Screen disease pairs for influential impact interactions
- [`compare_scenarios()`](https://rasmussenphilip.github.io/deconflate/reference/compare_scenarios.md)
  : Compare scenarios

## Plots

- [`plot(`*`<cm_result>`*`)`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)
  [`plot_burden()`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)
  [`plot(`*`<cm_mc>`*`)`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)
  [`plot(`*`<cm_screen>`*`)`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)
  [`plot(`*`<cm_oat>`*`)`](https://rasmussenphilip.github.io/deconflate/reference/plots.md)
  : Plots

## Examples and validation

- [`example_supplement()`](https://rasmussenphilip.github.io/deconflate/reference/example_supplement.md)
  : Worked example from the Supplementary File of Rasmussen et al.
  (2022)
- [`example_uk_dairy_2022()`](https://rasmussenphilip.github.io/deconflate/reference/example_uk_dairy_2022.md)
  : UK dairy example from Rasmussen et al. (2022)
- [`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md)
  : Economic inputs for the UK dairy example (Rasmussen et al. 2022,
  Table 1)
- [`example_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/example_global_dairy.md)
  : Global dairy inputs from Rasmussen et al. (2024), at their central
  values
- [`sampler_global_dairy()`](https://rasmussenphilip.github.io/deconflate/reference/sampler_global_dairy.md)
  : Monte Carlo sampler for the global dairy inputs (Rasmussen et al.
  2024)
- [`simulate_raw_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/simulate_raw_impacts.md)
  : Simulate the raw impacts a single-disease study would report
