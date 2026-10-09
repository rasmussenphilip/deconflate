#' deconflate: comorbidity adjustment of disease impact estimates
#'
#' Single-disease impact estimates (e.g. the milk-yield reduction in cows with
#' lameness compared with cows without) are conflated with the impacts of
#' associated diseases, which are more (or less) common among affected
#' animals. Summing such estimates double counts. `deconflate` adjusts the
#' estimates so that they can be aggregated.
#'
#' @section Workflow:
#' 1. Read the inputs with [cm_read_inputs()]: the diseases, their
#'    associations, one impact table and, optionally, interactions and
#'    three-way terms (CSV files with any names, or data frames; see
#'    [cm_template()]). Values can have distributions in their own rows. The
#'    same inputs can be built in R with [cm_diseases()],
#'    [cm_associations()], [cm_three_way()], [cm_impacts()],
#'    [cm_interactions()] and [cm_model()]. For several outcomes, use one
#'    impact table each and adjust each in turn.
#' 2. Adjust the impacts with [deconflate()]. Additive impacts (in any units)
#'    are adjusted with the exact simultaneous solution when every pair has
#'    an association and there are no interactions or three-way terms, and
#'    with the global (maximum-entropy) model otherwise; pairs without an
#'    association are unknown and filled in by the global model. Event
#'    impacts (hazard ratios, rate ratios, risk ratios, odds ratios or risk
#'    differences of e.g. culling) use the snapshot hazard model
#'    (`event_model = TRUE`, with the overall risk), which also gives the
#'    risk attributable to disease. When inputs have distributions,
#'    `deconflate()` also reports intervals (`n_draws`).
#' 3. Inspect each disease's contribution ([attribute_burden()],
#'    [contribution_table()], [plot_burden()]), in the units of the impacts.
#'    Converting the results into other quantities (e.g. a productivity gap
#'    or a monetary value) is left to the user. [compare_methods()] sets the
#'    methods side by side, including the proportional approximation of
#'    Rasmussen et al. (2022, eq. 16), which is kept for comparison and
#'    reproduction only.
#' 4. Check joint feasibility with [check_feasibility()] (a screen runs by
#'    default).
#' 5. Explore the inputs with [sensitivity_oat()], [screen_associations()]
#'    (also without any association estimates), [screen_interactions()] and
#'    [screen_three_way()], and find where a conclusion (a ranking, a sign
#'    or the total) changes as one input varies with [cm_threshold()].
#'
#' [reproduce_rasmussen_2022()] and [reproduce_rasmussen_2024()] recompute
#' the published tables. See `vignette("deconflate")` to get started.
#'
#' @section Conventions:
#' * Additive impacts are in any units, the same within an impact table;
#'   results come back in those units. Event impacts are adjusted to hazard
#'   ratios, and risks are proportions of animals over the period.
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
