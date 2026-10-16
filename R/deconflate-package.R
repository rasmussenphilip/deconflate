#' deconflate: comorbidity adjustment of disease impact estimates
#'
#' Single-disease impact estimates (e.g. the milk-yield reduction in cows with
#' lameness compared with cows without) are conflated with the impacts of
#' associated diseases, which are more (or less) common among affected
#' animals. Summing such estimates double counts. `deconflate` adjusts the
#' estimates so that they can be aggregated.
#'
#' @section Workflow:
#' 1. Read the inputs with [cm_read_inputs()]: the diseases, the association
#'    rows and one impact table (CSV files with any names, or data frames;
#'    see [cm_template()]). An association row can be about a pair, a
#'    combination or a subgroup, crude or adjusted; an impact row can be
#'    about a disease or a combination (an interaction), within a subgroup,
#'    crude or adjusted; several rows can describe the same thing. Values can
#'    have distributions in their own rows. The same inputs can be built in R
#'    with [cm_diseases()], [cm_associations()], [cm_impacts()] and
#'    [cm_model()]. For several outcomes, use one impact table each and
#'    adjust each in turn.
#' 2. Adjust the impacts with [deconflate()]. The package has one model: the
#'    maximum-entropy distribution of disease combinations that matches every
#'    association row (pairs that no row mentions are unknown and filled in),
#'    and the impacts and interactions that reproduce every impact row
#'    (exactly, or by weighted least squares when there are more rows than
#'    unknowns). When the pairwise tables give the same answer, they are used
#'    directly. Event impacts (hazard ratios, rate ratios, risk ratios, odds
#'    ratios or risk differences of e.g. culling) use the snapshot hazard
#'    model (`event_model = TRUE`, with the overall risk or the risk of
#'    disease-free animals), which also gives
#'    the risk attributable to disease. When inputs have distributions,
#'    `deconflate()` also reports intervals (`n_draws`).
#' 3. Inspect each disease's contribution ([contribution_table()],
#'    [summary()][summary.cm_result]), in the units of the impacts, and how
#'    well the rows fit together (`$estimates`, `$heterogeneity`). Converting
#'    the results into other quantities (e.g. a productivity gap or a
#'    monetary value) is left to the user.
#' 4. Check joint feasibility with [check_feasibility()].
#' 5. See which inputs matter, entered or not (unknown associations,
#'    interactions), with [cm_sensitivity()] (also without any association
#'    estimates); how often each disease ranks first over the draws in
#'    `$ranks`; where a conclusion (a ranking, a sign, or the loss of a set
#'    of diseases) changes as an input varies with [cm_threshold()]; and how
#'    likely it is over the draws of the other inputs with
#'    [cm_threshold_prob()].
#'
#' See `vignette("deconflate")` to get started.
#'
#' @section Conventions:
#' * Additive impacts are in any units, the same within an impact table;
#'   results come back in those units. Event impacts are adjusted to hazard
#'   ratios, and risks are proportions of animals over the period.
#' * `E[k, i]` is the excess probability of disease `k` among animals with
#'   disease `i`: `P(k | i) - P(k | not i)`.
#' * Under additive impacts, the crude impact of disease `i` satisfies
#'   `m_raw[i] = m[i] + sum_k E[k, i] * m[k]` exactly.
#' * Errors are condition objects of class `deconflate_error`. Failures that
#'   code may want to handle separately have a subclass as well:
#'   `deconflate_infeasible` (inputs that no population can have),
#'   `deconflate_singular` (an estimate that is not identifiable, or rows
#'   that do not determine the impacts), `deconflate_nonconvergence` (a fit
#'   that did not converge), `deconflate_unsupported` (a combination of
#'   inputs and options that is not supported), `deconflate_nonfinite`
#'   (non-finite results) and `deconflate_invalid_input` (an argument, table,
#'   row or input value that is not valid). Problems found by
#'   [cm_read_inputs()] have class `deconflate_input_problems`
#'   (also `deconflate_invalid_input`).
#' * Warnings are condition objects of class `deconflate_warning`, with the
#'   subclasses `deconflate_sign_change` (an adjusted impact on the other
#'   side of zero, or an adjusted hazard ratio on the other side of 1, from
#'   its raw estimate), `deconflate_nonfinite_warning` (non-finite adjusted
#'   impacts), `deconflate_nonconvergence` (a joint distribution that did not
#'   converge) and `deconflate_incomplete_allocation` (a Shapley allocation
#'   that skipped combinations).
#' * Messages are of class `deconflate_message`; the notes of
#'   [cm_read_inputs()] are also `deconflate_input_notes`.
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
