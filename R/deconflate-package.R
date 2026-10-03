#' deconflate: comorbidity adjustment of disease impact estimates
#'
#' Single-disease impact estimates (e.g. the milk-yield reduction in cows with
#' lameness compared with cows without) are conflated with the impacts of
#' associated diseases, which are more (or less) common among affected
#' animals. Summing such estimates double counts. `deconflate` adjusts the
#' estimates so that they can be aggregated.
#'
#' @section Workflow:
#' 1. Describe the population with [cm_diseases()], [cm_associations()] and,
#'    optionally, three-way scenarios ([cm_three_way()]), combined with
#'    [cm_population()]. Add one impact vector ([cm_impacts()]) per analysis
#'    with [cm_model()], or several with [cm_analyses()]. Or read the same
#'    tables from CSV files or data frames with [cm_read_inputs()] (see
#'    [cm_template()]).
#' 2. Adjust impacts with [deconflate()], using one of three methods:
#'    * `"simultaneous"` (default): the exact solution of the additive impact
#'      equations, built from pairwise 2x2 tables. Crude estimates and
#'      coefficients from additive regressions with a recorded adjustment set
#'      (`estimand = "adjusted_linear"`) are supported.
#'    * `"published"`: the proportional approximation of Rasmussen et al.
#'      (2022), eq. 16, for crude estimates, for reproduction and comparison.
#'    * `"global"`: fits a maximum-entropy distribution of disease
#'      combinations with [fit_joint()], then solves the impact equations
#'      including pairwise interactions ([cm_interactions()]).
#' 3. Inspect each disease's contribution to the aggregate
#'    ([attribute_burden()], [contribution_table()]); optionally turn the
#'    aggregate into a productivity gap and a value with
#'    [productivity_gap()] and [value_losses()]. [compare_methods()]
#'    tabulates the methods side by side.
#' 4. Hazard ratios (e.g. of culling) combine multiplicatively and have their
#'    own adapter: [cm_hazard_ratios()], [cm_hr_model()], [deconflate_hr()]
#'    and [attributable_risk()].
#' 5. Check joint feasibility with [check_feasibility()] (a screen runs by
#'    default).
#' 6. Propagate uncertainty with [cm_sampler()] or [cm_batch_sampler()] and
#'    [cm_monte_carlo()] (check stability with [cm_diagnose()]), and explore
#'    scenarios with [cm_scenario()], [sensitivity_oat()],
#'    [screen_associations()], [screen_interactions()], [screen_three_way()]
#'    and [compare_scenarios()].
#'
#' [reproduce_rasmussen_2022()] and [reproduce_rasmussen_2024()] recompute
#' the published tables. See `vignette("deconflate")` to get started.
#'
#' @section Conventions:
#' * Impacts are additive and in any units, the same within an analysis;
#'   results come back in those units.
#' * `E[k, i]` (see [excess_matrix()]) is the excess probability of disease
#'   `k` among animals with disease `i`: `P(k | i) - P(k | not i)`.
#' * Under additive impacts, the crude impact of disease `i` satisfies
#'   `m_raw[i] = m[i] + sum_k E[k, i] * m[k]` exactly.
#' * Errors have classes `deconflate_infeasible`, `deconflate_singular`,
#'   `deconflate_nonconvergence`, `deconflate_unsupported` and
#'   `deconflate_nonfinite` (all also `deconflate_error`); problems found
#'   by [cm_read_inputs()] have class `deconflate_input_problems`. Warning classes
#'   include `deconflate_sign_change`, `deconflate_nonconvergence` and
#'   `deconflate_nonfinite_warning` (all also `deconflate_warning`).
#'
#' @references
#' Rasmussen P, Shaw APM, Munoz V, Bruce M, Torgerson PR (2022). Estimating the
#' burden of multiple endemic diseases and health conditions using Bayes'
#' Theorem: a conditional probability model applied to UK dairy cattle.
#' Preventive Veterinary Medicine 203:105617.
#'
#' Rasmussen P, et al. (2024). Global losses due to dairy cattle diseases: a
#' comorbidity-adjusted economic analysis. Journal of Dairy Science
#' 107:6945-6970.
#'
#' @keywords internal
"_PACKAGE"
