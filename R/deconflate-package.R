#' deconflate: comorbidity adjustment of disease impact estimates
#'
#' Single-disease impact estimates (e.g. the milk-yield reduction in cows with
#' lameness compared with cows without) are conflated with the impacts of
#' associated diseases, which are more (or less) common among affected
#' animals. Summing such estimates double counts. `deconflate` adjusts the
#' estimates so that they can be aggregated.
#'
#' @section Workflow:
#' 1. Describe the system with [cm_diseases()], [cm_associations()],
#'    [cm_impacts()] and, optionally, [cm_interactions()], then combine them
#'    with [cm_model()]. Or read the same tables from CSV files or data
#'    frames with [cm_read_inputs()] (see [cm_template()]).
#' 2. Adjust impacts with [deconflate()], using one of three methods:
#'    * `"simultaneous"` (default): the exact solution of the additive impact
#'      equations, built from pairwise 2x2 tables.
#'    * `"published"`: the proportional approximation of Rasmussen et al.
#'      (2022), eqs. 15-16, for reproduction and comparison.
#'    * `"global"`: fits a maximum-entropy distribution of disease combinations
#'      with [fit_joint()], then solves the impact equations including pairwise
#'      interactions.
#' 3. Estimate, value and attribute productivity gaps with
#'    [productivity_gap()], [value_losses()], [attribute_burden()] and
#'    [contribution_table()]. Culling hazard ratios can be adjusted directly
#'    (`scale = "hazard_ratio"`) and turned into culling attributable to
#'    disease with [attributable_risk()]; [compare_methods()] tabulates the
#'    methods side by side.
#' 4. Check joint feasibility with [check_feasibility()].
#' 5. Propagate uncertainty with [cm_sampler()] and [cm_monte_carlo()] (check
#'    stability with [cm_diagnose()]), and
#'    explore scenarios with [cm_scenario()], [sensitivity_oat()],
#'    [screen_associations()], [screen_interactions()] and
#'    [compare_scenarios()].
#'
#' See `vignette("deconflate")` to get started.
#'
#' @section Conventions:
#' * Impacts on the `"proportion"` scale are proportional changes relative to
#'   the disease-free value (e.g. `0.025` for a 2.5% yield reduction).
#' * `E[k, i]` (see [excess_matrix()]) is the excess probability of disease
#'   `k` among animals with disease `i`: `P(k | i) - P(k | not i)`.
#' * Under additive impacts, the raw (crude) impact of disease `i` satisfies
#'   `m_raw[i] = m[i] + sum_k E[k, i] * m[k]` exactly.
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
