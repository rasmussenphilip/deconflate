#' Adjust raw impact estimates for comorbidity
#'
#' Raw impact estimates (comparisons of animals with and without a disease)
#' are conflations of the disease's own impact and the impacts of associated
#' diseases. `deconflate()` removes the part of each estimate that belongs to
#' other diseases, for one impact table at a time: additive impacts (the
#' default), or event impacts (`event_model = TRUE`). When inputs have
#' distributions, it also runs draws (`n_draws`) and reports intervals.
#'
#' @section The model:
#' The package has one model. The joint distribution of disease combinations
#' is the maximum-entropy distribution that matches the disease
#' probabilities and every association row ([fit_joint()]); pairs that no
#' row mentions are unknown and are filled in by the fit (listed in
#' `$unknown_pairs`). Every impact row is then an equation in the unknown
#' impacts: an impact for each disease and an interaction for each
#' combination that appears in the impact table (see [cm_impacts()] for what
#' each row means).
#'
#' When the impacts are additive and the inputs are only crude pair
#' associations (one for every pair) and impact rows about
#' single diseases, the pairwise 2x2 tables give exactly the same answer
#' (the simultaneous solution of Rasmussen et al. 2022, eq. 14, generalised
#' to adjusted estimates), so the package uses them ("pairwise shortcut",
#' recorded in the result), and fits the joint distribution once to confirm
#' that the associations can hold together (up to 20 diseases). The draws
#' then use the shortcut too.
#'
#' @section Additive impacts:
#' With impacts `b` and interactions `delta`, the outcome of an animal with
#' disease combination `c` is `y(c) = sum_d b_d c_d + sum_T delta_T all_T(c)`.
#' Each impact row is a linear function of `b` and `delta`: for a crude row,
#' the difference in mean outcome between animals with all the row's
#' diseases and animals with none of them (within `among`); for an adjusted
#' row, a regression coefficient. Results are in the units of the impacts.
#'
#' The expected aggregate impact per animal is
#' `sum_d p_d b_d + sum_T P(T) delta_T`. Each disease's contribution is its
#' own term plus an equal share of each interaction it is part of (the
#' closed-form Shapley value); contributions add up to the aggregate. Shares
#' are `NA` when the aggregate is zero.
#'
#' @section Several rows per item:
#' Impact rows are fitted together: with one row per unknown the equations
#' are solved exactly; with more rows, by weighted least squares with weights
#' 1/SD^2 (SDs from each row's distribution; a row without a distribution is
#' exact and is reproduced exactly). The result reports each row's fitted
#' value and standardised residual (`$estimates`) and the heterogeneity
#' statistic Q (`$heterogeneity`). Two exact rows that contradict each other
#' are an error.
#'
#' Association rows that describe the same quantity are pooled on the log
#' odds ratio scale of their 2x2 table (the logit of a prevalence, the log
#' of an `OR_ratio`, or, for rows of different kinds about the same
#' quantity, the logit of its probability), inverse-variance weighted, with random effects (DerSimonian-Laird) by
#' default and `pool = "fixed"` as an option; Q, I^2 and tau^2 are in
#' `$heterogeneity`. A row without a distribution is exact and wins. In the
#' draws, each row is drawn from its own distribution and pooled with the
#' central weights, and a between-study term keeps the pooled variance at
#' its random-effects value.
#'
#' @section Event impacts:
#' With `event_model = TRUE`, the impact table holds event impacts (hazard
#' ratios, rate ratios, risk ratios, odds ratios or risk differences; see
#' [cm_impacts()]), which are adjusted with the snapshot hazard model. Within
#' the period of the overall risk, an animal with disease combination `c` has
#' a constant hazard `h0 * exp(sum_d beta_d c_d + sum_T gamma_T all_T(c))`, so
#' its risk over the period is `R(c) = 1 - exp(-h0 * exp(...))`; `h0` is set
#' so that the population risk equals `overall_risk`, or so that the risk of
#' disease-free animals equals `baseline_risk` (give one of the two). The `beta`s and
#' `gamma`s (interactions on the log-hazard scale) are solved so that each
#' raw estimate is reproduced: a hazard or rate ratio by the ratio of the
#' average hazard multipliers of the two groups the row compares; a risk
#' ratio, odds ratio or risk difference by those of their average risks;
#' crude, or within strata of the row's adjustment set (combined with
#' Mantel-Haenszel-type weights). The adjusted hazard ratios `exp(beta)` are
#' each disease's own hazard multiplier. For each disease (and combination),
#' `$adjusted` also gives the risk of an animal with only that disease
#' (`risk_alone`) and its risk ratio, risk difference and odds ratio against
#' a disease-free animal (`rr_alone`, `rd_alone`, `or_alone`): the adjusted
#' result on the scale of risk-based estimates. `change` compares the
#' adjusted result with the raw estimate on the raw estimate's own measure.
#' The risk attributable to disease is the overall risk minus the
#' disease-free risk `1 - exp(-h0)`, allocated to diseases by Shapley values
#' over disease combinations (cell by cell) and split, as for additive
#' impacts, into each disease's own part and its share of the interactions.
#' `$conflation$A` holds the row equations at the adjusted values: how each
#' row's value (on its measure's scale) changes with each log hazard
#' multiplier (a local version of the additive row equations; for crude
#' hazard-ratio rows at hazard ratios of 1, they equal the additive ones).
#'
#' What the snapshot model is not: a Cox hazard ratio estimated over
#' follow-up is not, in general, the snapshot ratio, because animals with
#' high hazards leave first and the mixture of disease combinations among
#' survivors changes. The model is exact for its own comparisons and a
#' reasonable approximation when follow-up is short relative to the hazards
#' or the diseases are rare. Risk-based estimates (risk ratios, odds ratios,
#' risk differences) must refer to the same period as the overall or baseline
#' risk.
#'
#' @section Raw and adjusted totals:
#' The raw total is the result of the same model with the diseases unrelated
#' (the same disease values and impact rows, no association rows): for
#' additive impacts with one crude row per disease, it is `sum(p * raw)`;
#' for event impacts, the risk attributable to disease computed from the raw
#' estimates. The difference from the adjusted total is the double counting
#' removed. Above 20 diseases it is computed only for additive impacts with
#' single-disease rows without a subgroup (otherwise `NA`).
#'
#' @section Uncertainty:
#' Inputs with a distribution (the `dist` columns of the input tables, or
#' `distributions` in [cm_model()]; and `overall_risk` or `baseline_risk`
#' given as a distribution) are drawn `n_draws` times; each draw is adjusted in the same
#' way as the central estimate. The central estimate uses the point values
#' (`value`, and the mean of a risk given as a distribution); the draws give
#' 95% intervals (2.5% and 97.5% quantiles), means, Monte Carlo standard
#' errors and stability checks (`$draws$summary`). Draws whose values cannot
#' hold together are rejected and counted (`$draws$rejections`), not
#' replaced; notes in `$notes` report a high rejection share and limited
#' Monte Carlo precision. Without distributions, no draws are run.
#'
#' @section Diagnosing the inputs:
#' An adjusted impact on the other side of zero from its raw impact (for
#' event impacts, of 1), or one that has shrunk by more than 90%, usually
#' means that the inputs contradict each other under the model: a raw impact
#' smaller than what the associated diseases alone account for, an
#' association or prevalence that is too strong, or estimates that measure
#' different things (e.g. "disease alone vs healthy" entered as a crude
#' comparison). With `diagnose = TRUE`, `$problems` shows which inputs
#' drive such flagged impacts:
#' * for each input with a distribution, the change in each flagged impact
#'   when the input moves up by one SD (a shift of 1 in its normal score,
#'   `qnorm(F(x))`), and the shift that alone would remove the flag, when it
#'   is within 3 SD (from the local slope);
#' * for each input without a distribution, the value at which the flag
#'   disappears, searched from a quarter to three times its value (for RD
#'   and phi, within 0.5 of it); for additive impact rows it is solved
#'   exactly, and always shown for the flagged disease's own rows.
#'
#' The flag disappears when the adjusted impact reaches 10% of the raw impact
#' as entered (on the log scale for ratios; for risk difference rows, when it
#' crosses to the raw impact's side). With several rows for a disease, the
#' raw impact is their inverse-variance weighted mean (for event impacts, of
#' the ratio rows when there are any). Thresholds outside a quarter to three
#' times the value (shown only for a disease's own rows) are marked, and do
#' not count as a fix. The listed inputs are
#' those that alone remove a flag within 3 SD, and the four with the largest
#' effect on each flagged impact. Association pooling is fixed at the
#' central estimate, as in the draws.
#'
#' A separate table lists rows that disagree with other rows for the same
#' quantity (standardised residual of 2 or more), with the adjusted impact
#' that changes most and the total (`total_before`, `total_after`) without
#' each row.
#'
#' The check describes the central estimate; it does not change the results.
#' A flagged input is worth checking, not necessarily wrong.
#'
#' @section Errors and warnings:
#' Errors and warnings are condition objects with classes (see
#' [deconflate-package]). Those about particular inputs also hold them as
#' fields, so that code can use them without reading the message:
#' * `table` and `row` (or `rows`): the input rows concerned, e.g. an impact
#'   row that is not identifiable, or exact rows that contradict each other
#'   (with their `residuals` for impact rows);
#' * `terms`: diseases and interactions, e.g. those the impact rows do not
#'   determine (`deconflate_singular`), or those whose adjusted impact
#'   changes sign (the warning `deconflate_sign_change`);
#' * `triples`: the conflicting triples when the crude pair rows cannot all
#'   hold (`deconflate_infeasible`), as in [check_feasibility()];
#' * `residual`: how far a fit was from converging
#'   (`deconflate_nonconvergence`).
#'
#' @param model A [cm_model()] (e.g. from [cm_read_inputs()]) with at least
#'   one association.
#' @param event_model `FALSE` (default) for additive impacts; `TRUE` for
#'   event impacts (an impact table with a `measure` column).
#' @param overall_risk With `event_model = TRUE` (this or `baseline_risk` is
#'   required): the overall risk of the event in the population over the
#'   period, as a proportion (e.g. `0.25`), or a distribution (e.g.
#'   `dist_beta(250, 750)`).
#' @param baseline_risk With `event_model = TRUE`, instead of `overall_risk`:
#'   the risk of the event in disease-free animals over the period, as a
#'   proportion or a distribution. The overall risk then follows from the
#'   model.
#' @param n_draws Number of draws for the uncertainty (default 1000; 0 for
#'   point estimates only; otherwise at least 2). Used only when some input
#'   has a distribution.
#' @param seed Optional random seed for the draws and for a sampled joint
#'   distribution (one is chosen and stored in the result otherwise).
#' @param sampling `"random"` (default) or `"lhs"` (Latin hypercube, in
#'   `lhs_replicates` independent blocks; the Monte Carlo error is then
#'   estimated from the block means).
#' @param lhs_replicates Number of Latin hypercube blocks (at least 2; at
#'   most `n_draws / 2` are used).
#' @param pool `"random"` (default) or `"fixed"`: how several association
#'   rows about the same quantity are pooled.
#' @param joint Optional [fit_joint()] result. It is checked against the
#'   model (diseases, probabilities and association rows).
#' @param warn Logical: warn when adjusted impacts change sign or are not
#'   finite?
#' @param diagnose Logical (default `FALSE`): show which inputs drive
#'   problematic adjusted impacts, and which rows disagree with other rows
#'   (`$problems`, printed with the result; see the section "Checking the
#'   inputs"). It reruns the central estimate many times (twice per input
#'   with a distribution, a few times per input without one).
#' @param groups Optional groups of diseases to report together, as a named
#'   list of disease ids, e.g. `list(ketosis = c("SCK", "CK"), mastitis =
#'   c("SCM", "CM"))` for subclinical and clinical forms. The diseases stay
#'   separate in the model; a group's contribution is the sum of its
#'   diseases' contributions (an interaction between two diseases of a
#'   group is counted wholly in the group), with an interval from the sum
#'   within each draw. A disease can be in one group only, and a group
#'   cannot have the name of a disease. See the section "Groups of
#'   diseases".
#' @param ... Passed to [fit_joint()] (e.g. `backend = "sampled"`). With
#'   more than 20 diseases the sampled backend is used.
#'
#' @section Groups of diseases:
#' Forms of one disease (e.g. subclinical and clinical ketosis) have their
#' own prevalences, associations and impacts in the literature, so they are
#' best entered as separate diseases and added up afterwards with `groups`.
#' Contributions are additive, so the group totals and the diseases not in a
#' group add up to the aggregate (for event impacts, to the allocated
#' attributable risk). How the
#' forms overlap (e.g. clinical cases that also count as subclinical, or
#' each animal counted in its worst form only) is set by the association
#' row between them.
#'
#' @return For additive impacts, a `cm_result` with elements:
#'   * `adjusted`: one row per disease and per interaction: the raw value
#'     (of its row, when it has exactly one), the adjusted value, the
#'     relative change and the number of rows (with draws, also `lower` and
#'     `upper`);
#'   * `totals`: the raw aggregate (`raw_total`: the diseases unrelated; see
#'     "Raw and adjusted totals"), the adjusted aggregate and its
#'     interaction part;
#'   * `contributions`: per disease, its `contribution` to the aggregate
#'     (Shapley value), split into `main` and `interaction` parts, its
#'     `share` (with draws, also `lower` and `upper`), and its
#'     `raw_contribution` (the same with the diseases unrelated: they add up
#'     to the raw total);
#'   * `estimates`: every input row (impacts and associations) with its
#'     value, fitted value and standardised residual (and, for impact rows,
#'     the SD used as weight; 0 for an exact row);
#'   * `heterogeneity`: pooled association rows and the fit of the impact
#'     rows (Q, df, p, I^2, tau^2), or `NULL`;
#'   * `diagnostics`: the computation used (`"full model"` or
#'     `"pairwise shortcut"`), reconstruction residual, rank and condition
#'     number of the row equations, sign changes and the feasibility check;
#'   * `unknown_pairs`: pairs that no row mentions and the odds ratios the
#'     fit gave them;
#'   * `draws`: the uncertainty (summary, accepted draws, rejections, seed),
#'     or `NULL`;
#'   * `ranks`: with draws, how often each disease ranks 1st, 2nd, ... by its
#'     contribution (`rank_1`, `rank_2`, ...), and its mean rank;
#'   * `groups`: with `groups`, one row per group and per disease not in a
#'     group (the rows add up to the sum of all contributions): `group`, `diseases`,
#'     `raw_contribution`, `contribution`, `lower`, `upper` (from the draws,
#'     else `NA`), `main`, `interaction` and `share`; `group_ranks`: with
#'     draws, as `ranks`, for these rows;
#'   * `problems`: with `diagnose = TRUE`, a `cm_diagnosis` object:
#'     `flagged` (every disease with impact rows: raw, adjusted, the
#'     criterion and whether it is flagged), `inputs` (inputs with a
#'     distribution: effect per SD and SD needed, for each flagged disease;
#'     `listed` marks those printed), `thresholds` (inputs without a
#'     distribution), `no_single_fix`, `disagreeing` and `settings`;
#'   * `notes`, `conflation` (the row equations), `joint_pairs` (matrix of
#'     `P(j and k)`), `model`, `joint`, `method`, `outcome`, `units` and
#'     `settings` (what [cm_sensitivity()] and [cm_threshold()] need to
#'     rerun the model).
#'
#'   For event impacts, a `cm_event_result` with `adjusted` (adjusted hazard
#'   ratios of diseases and interactions, the change from the raw estimate,
#'   and the risk-scale results `risk_alone`, `rr_alone`, `rd_alone` and
#'   `or_alone`), `contributions` (as for additive impacts: each disease's
#'   share of the attributable risk, Shapley allocation, with its `main` and
#'   `interaction` parts and its `raw_contribution`, the diseases
#'   unrelated), `attributable` (`summary`: overall, disease-free,
#'   raw and adjusted attributable risk and the attributable fraction;
#'   `baseline_hazard`), `conflation` (the local row equations), `estimates`,
#'   `heterogeneity`, `diagnostics`, `unknown_pairs`, `draws`, `ranks` (by
#'   share of the attributable risk), `groups` and `group_ranks` (with
#'   `groups`), `problems` (with `diagnose = TRUE`), `notes`, `joint`, `model`, `method`,
#'   `outcome`, `units`, `overall_risk`, `risk_given` (`"overall"` or `"baseline"`) and
#'   `settings`. Risks are
#'   proportions of animals with the event during the period.
#' @export
#' @examples
#' res <- deconflate(example_supplement())
#' res$adjusted
#' res$totals
#'
#' # An estimate adjusted for one other disease, and a second estimate for d1
#' m <- example_supplement()
#' m$impacts <- cm_impacts(c("d1", "d1", "d2", "d3"), c(2.2, 2.6, 5, 7.5),
#'                         adjusted_for = c("d2", NA, NA, NA), units = "%")
#' m$distributions <- list("impacts:1" = dist_normal(2.2, 0.4),
#'                         "impacts:2" = dist_normal(2.6, 0.5))
#' deconflate(m, n_draws = 0)$estimates
#'
#' # Uncertain inputs
#' m2 <- cm_model(example_supplement(), example_supplement()$impacts,
#'                distributions = list("impact:d1" = dist_normal(2.5, 0.5),
#'                                     "assoc:d2:d3" = dist_lognormal_ci(3, 2, 4.5)))
#' deconflate(m2, n_draws = 200, seed = 1)
#'
#' # Event impacts (culling): hazard ratios and a risk ratio
#' cull <- cm_model(example_supplement(),
#'                  cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.2), measure = c("HR", "HR", "RR"),
#'                             outcome = "culling"))
#' deconflate(cull, event_model = TRUE, overall_risk = 0.25)
#'
#' # Errors hold the inputs they are about
#' bad <- cm_model(cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
#'                               cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05))),
#'                 cm_impacts(c("a", "b", "c"), c(1, 1, 1)))
#' e <- tryCatch(deconflate(bad), deconflate_infeasible = function(e) e)
#' e$triples
deconflate <- function(model, ...) UseMethod("deconflate")

#' @export
deconflate.default <- function(model, ...) {
  cm_abort("`model` must be created with cm_model() or cm_read_inputs().", class = "deconflate_invalid_input")
}

#' @rdname deconflate
#' @export
deconflate.cm_model <- function(model, event_model = FALSE, overall_risk = NULL, n_draws = 1000,
                                seed = NULL, sampling = c("random", "lhs"), lhs_replicates = 10L,
                                pool = c("random", "fixed"), joint = NULL, warn = TRUE,
                                diagnose = FALSE, baseline_risk = NULL, groups = NULL, ...) {
  dots <- list(...)
  if (!is.logical(diagnose) || length(diagnose) != 1L || is.na(diagnose)) {
    cm_abort("`diagnose` must be TRUE or FALSE.", class = "deconflate_invalid_input")
  }
  retired_args(dots)
  sampling <- cm_match_arg(sampling)
  pool <- cm_match_arg(pool)
  check_draw_args(n_draws, seed, sampling, lhs_replicates)
  n_draws <- as.integer(n_draws)
  check_event_model(model, event_model)
  groups <- check_groups(groups, model$diseases$disease)
  risk <- check_overall_risk(overall_risk, event_model, baseline_risk)
  a <- model$associations
  if (is.null(a) || !nrow(a)) {
    cm_abort("No association estimates were given, so there is nothing to de-conflate (the diseases would be treated as independent). To see how much associations could change the results, use cm_sensitivity(model) or cm_threshold().",
             class = "deconflate_unsupported")
  }
  plan <- plan_model(model, event_model, dots, has_joint = !is.null(joint), pool = pool,
                     risk_kind = risk$kind %||% "overall")
  # The seed also makes a sampled joint distribution reproducible.
  if (!is.null(seed) && identical(plan$fit_args$backend %||% dots$backend, "sampled") && is.null(dots$seed)) {
    plan$fit_args$seed <- seed
  }
  res <- run_point(model, plan, risk = risk$value, joint = joint, warn = warn, dots = dots, central = TRUE)
  res <- finish_result(res, model, plan)
  # What cm_sensitivity(), cm_threshold() and the diagnosis need to rerun it.
  res$settings <- list(plan = plan, risk = risk, dots = dots, groups = groups)
  if (diagnose) res$problems <- diagnose_inputs(model, plan, res, risk, dots)
  specs <- draw_specs(model, risk$dist, res$pooling$frozen)
  if (n_draws > 0L) {
    if (!length(specs)) {
      res$notes <- c(res$notes, "No input has a distribution, so no draws were run: the results are point estimates.")
    } else {
      res$draws <- deconflate_draws(model, specs, plan, res, risk$value, n_draws, seed, sampling,
                                    lhs_replicates, dots)
      res <- attach_intervals(res)
      res$ranks <- draw_ranks(res)
      res$notes <- c(res$notes, draw_notes(res$draws, event = inherits(res, "cm_event_result")))
    }
  } else {
    n_unc <- sum(!startsWith(names(specs) %||% character(0), "between:"))
    if (n_unc) {
      res$notes <- c(res$notes, sprintf("%d input%s a distribution, but n_draws = 0: the results are point estimates without intervals.",
                                        n_unc, if (n_unc == 1L) " has" else "s have"))
    }
  }
  if (!is.null(groups)) {
    gr <- group_results(res, groups)
    res$groups <- gr$table
    res$group_ranks <- gr$ranks
  }
  res
}

# Arguments of earlier versions.
retired_args <- function(dots) {
  nm <- names(dots) %||% character(0)
  if ("method" %in% nm) {
    cm_abort("The `method` argument has been removed: deconflate() always gives the model's answer, and uses the pairwise tables itself when they give the same answer (see ?deconflate). The published approximation of Rasmussen et al. (2022, 2024) was removed in version 0.6.",
             class = "deconflate_unsupported")
  }
  if ("feasibility" %in% nm) {
    cm_abort("The `feasibility` argument has been removed: the joint distribution is fitted for the central estimate, which checks that the associations can hold together (see check_feasibility() for the separate checks).",
             class = "deconflate_unsupported")
  }
  invisible(TRUE)
}

# Additive or event impacts, as the call says? A mismatch is an error.
check_event_model <- function(model, event_model) {
  if (!is.logical(event_model) || length(event_model) != 1L || is.na(event_model)) {
    cm_abort("`event_model` must be TRUE or FALSE.", class = "deconflate_invalid_input")
  }
  if (is.null(model$impacts)) cm_abort("The model has no impacts to adjust.", class = "deconflate_invalid_input")
  kind <- impact_kind(model$impacts)
  if (event_model && kind != "event") {
    cm_abort("event_model = TRUE needs event impacts: an impact table with a measure column (HR, rate_ratio, RR, OR or RD; see ?cm_impacts).",
             class = "deconflate_unsupported")
  }
  if (!event_model && kind == "event") {
    cm_abort("These are event impacts (the impact table has a measure column): use event_model = TRUE, with overall_risk or baseline_risk.",
             class = "deconflate_unsupported")
  }
  invisible(TRUE)
}

# Can the pairwise tables give the answer? Additive impacts, every
# association row a crude pair row, every pair covered, every impact row
# about one disease without a subgroup.
shortcut_ok <- function(model) {
  if (impact_kind(model$impacts) == "event") return(FALSE)
  specs <- association_specs(model)
  if (!length(specs) || !all(vapply(specs, `[[`, logical(1), "pair"))) return(FALSE)
  im <- impact_specs(model)
  if (any(vapply(im, function(s) length(s$g) > 1L || length(c(s$apos, s$aneg)) > 0L, logical(1)))) return(FALSE)
  ids <- model$diseases$disease
  if (length(ids) < 2L) return(TRUE)
  covered <- unique(vapply(specs, function(s) pair_key(s$g, s$h), character(1)))
  length(covered) == choose(length(ids), 2)
}

# How to adjust: the pairwise shortcut or the full model, the arguments of
# the joint fit, and notes.
plan_model <- function(model, event_model, dots = list(), has_joint = FALSE, pool = "random",
                       risk_kind = "overall") {
  n <- nrow(model$diseases)
  notes <- character(0)
  fit_args <- list()
  if (n > 20L && is.null(dots$backend)) fit_args$backend <- "sampled"
  shortcut <- !event_model && shortcut_ok(model)
  if (!shortcut && !is.null(fit_args$backend) && !has_joint) {
    notes <- c(notes, sprintf("With %d diseases, the joint distribution was fitted with the sampled backend (see ?fit_joint).", n))
  }
  if (shortcut && n > 20L && !has_joint) {
    notes <- c(notes, sprintf("With %d diseases, the joint distribution was not fitted: the pairwise tables give the answer, and the triple screen (a necessary condition) checked that the associations can hold together.", n))
  }
  list(event = event_model, shortcut = shortcut, fit_args = fit_args, notes = notes, pool = pool,
       confirm = shortcut && (n <= 20L || has_joint), risk_kind = risk_kind)
}

# "a", "a and b", "a, b and c".
join_and <- function(x) {
  if (length(x) <= 1L) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

# One adjustment of a model, as planned. `central`: also fit the joint
# distribution to confirm the shortcut (draws use the shortcut alone), and
# compute the details of event results (the local row equations).
# `frozen`: pooling of association rows and impact-row SDs fixed at the
# central estimate; `between`: between-study deviates of a draw;
# `allocate`: the Shapley allocation of the attributable risk (event
# impacts); `raw`: also the raw total (see raw_total()).
run_point <- function(model, plan, risk = NULL, joint = NULL, warn = TRUE, dots = list(),
                      joint_start = NULL, start = NULL, frozen = NULL, central = FALSE, allocate = TRUE,
                      raw = TRUE) {
  fj <- utils::modifyList(dots[intersect(names(dots), joint_arg_names)], plan$fit_args)
  between <- model$between
  pool <- plan$pool %||% "random"
  if (isTRUE(plan$shortcut)) {
    pt <- pair_tables_at(model, pool, frozen$associations, between)
    tri <- triple_screen(pt, stats::setNames(model$diseases$prob, model$diseases$disease))
    if (nrow(tri)) {
      cm_abort(sprintf("The associations cannot all hold together: no population has all of them (conflicting triples: %s). See check_feasibility().",
                       paste(utils::head(paste(tri$disease1, tri$disease2, tri$disease3, sep = "-"), 5),
                             collapse = ", ")), class = "deconflate_infeasible", triples = tri)
    }
    res <- adjust_additive(model, pt = pt, warn = warn, frozen_se = frozen$impacts)
    res$pooling <- attr(pt, "pooling")
    res$diagnostics$feasibility <- "triple screen passed (necessary condition only)"
    if (central && isTRUE(plan$confirm)) {
      jt <- tryCatch(do.call(event_joint, c(list(as_population(model), joint, pool = pool), fj)),
                     deconflate_nonconvergence = function(e) {
                       cm_abort(sprintf("The associations cannot all hold together: %s", conditionMessage(e)),
                                class = "deconflate_infeasible", residual = e$residual)
                     })
      res$joint <- jt
      res$diagnostics$feasibility <- sprintf("joint distribution fitted (max residual %.1e)", jt$max_residual)
    }
  } else {
    jt <- do.call(event_joint, c(list(as_population(model), joint, joint_start = joint_start, pool = pool,
                                      frozen = frozen$associations, between = between), fj))
    res <- if (isTRUE(plan$event)) {
      adjust_event(model, overall_risk = risk, joint = jt, warn = warn, start = start, frozen_se = frozen$impacts,
                   allocate = allocate, risk_kind = plan$risk_kind %||% "overall", details = central)
    } else {
      r <- adjust_additive(model, joint = jt, warn = warn, frozen_se = frozen$impacts)
      r$diagnostics$feasibility <- sprintf("joint distribution fitted (max residual %.1e)", jt$max_residual)
      r
    }
    res$pooling <- jt$pooling
  }
  if (raw) {
    rt <- raw_total(model, plan, risk, res$row_se %||% frozen$impacts)
    if (isTRUE(plan$event)) {
      res$attributable$summary$raw_attributable <- rt
    } else {
      res$totals$raw_total <- rt
    }
    if (central && !is.null(res$contributions)) {
      # Each disease's raw contribution (for the reports and plots).
      rc <- raw_total(model, plan, risk, res$row_se %||% frozen$impacts, by_disease = TRUE)
      res$contributions$raw_contribution <- unname(rc[match(res$contributions$disease, names(rc))])
    }
  }
  res
}

# The raw total: the same model with the diseases unrelated (the same
# disease values and impact rows, no association rows). For additive impacts
# with one crude row per disease it is sum(p * raw); for event impacts it is
# the attributable risk computed from the raw estimates. NA when it cannot
# be computed (more than 20 diseases with rows other than single-disease
# rows without a subgroup, or a failed fit). `by_disease = TRUE`: each
# disease's raw contribution instead (named; they add up to the raw total).
raw_total <- function(model, plan, risk, frozen_se = NULL, by_disease = FALSE) {
  ids <- model$diseases$disease
  n <- length(ids)
  P <- stats::setNames(model$diseases$prob, ids)
  m0 <- model
  m0$associations <- NULL
  m0$between <- NULL
  d <- m0$distributions
  if (length(d)) m0$distributions <- d[startsWith(names(d), "impacts:") | startsWith(names(d), "prob:")]
  event <- isTRUE(plan$event)
  specs <- impact_specs(model)
  simple <- all(vapply(specs, function(s) length(s$g) == 1L && !length(c(s$apos, s$aneg)), logical(1)))
  none <- if (by_disease) stats::setNames(rep(NA_real_, n), ids) else NA_real_
  if (n > 20L && (event || !simple)) return(none)
  if (!event && simple) {
    # Unrelated diseases and single-disease rows: each row is its disease's
    # own impact; several rows are fitted together (weighted mean).
    se <- frozen_se %||% impact_row_se(model, specs)
    term <- vapply(specs, `[[`, character(1), "term")
    y <- vapply(specs, `[[`, numeric(1), "value")
    b <- vapply(ids, function(dd) {
      k <- which(term == dd)
      if (!length(k)) return(NA_real_)
      ex <- k[!is.na(se[k]) & se[k] == 0]
      if (length(ex)) return(y[ex[1]])
      w <- 1 / ifelse(is.na(se[k]), 1, se[k])^2
      sum(w * y[k]) / sum(w)
    }, numeric(1))
    if (anyNA(b)) return(none)
    return(if (by_disease) P * b else sum(P * b))
  }
  cells <- disease_cells(n)
  colnames(cells) <- ids
  prob <- rep(1, nrow(cells))
  for (j in seq_len(n)) prob <- prob * ifelse(cells[, j] == 1, P[j], 1 - P[j])
  jt <- structure(list(cells = cells, prob = prob, diseases = ids, converged = TRUE, iterations = 0L,
                       max_residual = 0, backend = "exact"), class = "cm_joint")
  out <- tryCatch(withCallingHandlers({
    r <- if (event) {
      adjust_event(m0, overall_risk = risk, joint = jt, warn = FALSE, allocate = by_disease, frozen_se = frozen_se,
                   risk_kind = plan$risk_kind %||% "overall", details = FALSE)
    } else {
      adjust_additive(m0, joint = jt, warn = FALSE, frozen_se = frozen_se)
    }
    if (by_disease) {
      stats::setNames(r$contributions$contribution[match(ids, r$contributions$disease)], ids)
    } else if (event) {
      r$attributable$summary$attributable
    } else {
      r$totals$adjusted_total
    }
  }, deconflate_nonconvergence = function(w) {
    if (inherits(w, "warning")) invokeRestart("muffleWarning")
  }), deconflate_error = function(e) none)
  out
}

# Diagnostics of the central estimate: every input row, heterogeneity,
# unknown pairs and notes.
finish_result <- function(res, model, plan) {
  jt <- res$joint
  pooling <- res$pooling
  assoc <- NULL
  if (!is.null(jt) && !is.null(jt$fit)) {
    assoc <- association_diagnostics(jt$fit, model)
    pooling <- jt$pooling
  } else if (!is.null(pooling$groups) && length(pooling$groups)) {
    assoc <- pair_row_diagnostics(pooling)
  }
  res$pooling <- pooling
  im <- res$estimates
  cols <- c("table", "row", "item", "measure", "among", "adjusted_for", "study", "value", "sd", "fitted",
            "weight", "std_residual")
  fill <- function(d) {
    if (is.null(d)) return(NULL)
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[, cols]
  }
  if (!is.null(assoc)) {
    a <- model$associations
    k <- match(assoc$row, a$row)
    assoc$among <- a$among[k]
    assoc$adjusted_for <- a$adjusted_for[k]
    assoc$sd <- NA_real_
  }
  res$estimates <- rbind(fill(im), fill(assoc))
  rownames(res$estimates) <- NULL
  het <- if (!is.null(pooling$groups) && !is.null(pooling$info)) {
    association_heterogeneity(pooling$groups, pooling$info, pooling$pool %||% plan$pool)
  }
  if (!is.null(het)) het <- cbind(table = "associations", het, stringsAsFactors = FALSE)
  f <- res$fit
  if (!is.null(f) && is.finite(f$Q)) {
    het <- rbind(het, data.frame(table = "impacts", item = "all impact rows", rows = "all", k = nrow(im),
                                 scale = if (inherits(res, "cm_event_result")) "measure scale" else "units",
                                 pooled = NA_real_, method = "weighted least squares", Q = f$Q, df = f$df,
                                 p = f$p, I2 = if (f$Q > 0) max(0, (f$Q - f$df) / f$Q) else 0, tau2 = NA_real_,
                                 stringsAsFactors = FALSE))
  }
  res$heterogeneity <- het
  res$notes <- c(plan$notes,
                 adjustment_notes(association_specs(model), "association"),
                 adjustment_notes(impact_specs(model), "impact"))
  if (!is.null(het) && any(het$table == "associations")) {
    h <- het[het$table == "associations", , drop = FALSE]
    res$notes <- c(res$notes, sprintf("Several association rows were pooled for %d item%s (%s; see $heterogeneity).",
                                      nrow(h), if (nrow(h) > 1L) "s" else "",
                                      if (all(h$method == "exact row")) "rows without uncertainty are exact" else
                                        sprintf("%s effects", plan$pool)))
  }
  if (!is.null(f) && isTRUE(f$df > 0)) {
    res$notes <- c(res$notes, sprintf("The impact rows over-determine the impacts (%d more rows than unknowns): fitted by weighted least squares (Q = %.3g on %d df, p = %.3g; see $estimates). Rows without uncertainty are exact.",
                                      f$df, f$Q %||% NA_real_, f$df, f$p %||% NA_real_))
  }
  res$unknown_pairs <- unknown_pairs_table(model, res$joint)
  res
}

# Diagnostics of crude pair rows from the pairwise tables (shortcut without
# a joint fit).
pair_row_diagnostics <- function(pooling) {
  rows <- list()
  for (j in seq_along(pooling$groups)) {
    gr <- pooling$groups[[j]]
    info <- pooling$info[[j]]
    for (k in seq_along(gr$specs)) {
      s <- gr$specs[[k]]
      se <- info$se_rows[k]
      rows[[length(rows) + 1L]] <- data.frame(
        table = "associations", row = s$row, item = gr$label, measure = s$measure, study = s$study,
        value = s$value, fitted = NA_real_,
        weight = if (length(gr$specs) > 1L) info$weights[k] / sum(info$weights) else 1,
        std_residual = if (isTRUE(se > 0) && length(gr$specs) > 1L) (info$z[k] - info$estimate) / se else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out[order(out$row), , drop = FALSE]
}

# Pairs that no association row mentions, and the odds ratios the fit gave
# them.
unknown_pairs_table <- function(model, joint) {
  pt <- pair_tables(model)
  pt <- pt[pt$status == "unknown", , drop = FALSE]
  if (!nrow(pt) || is.null(joint)) {
    return(data.frame(disease1 = pt$disease1, disease2 = pt$disease2,
                      fitted_or = rep(NA_real_, nrow(pt)), stringsAsFactors = FALSE))
  }
  ids <- joint$diseases
  J <- crossprod(joint$cells * joint$prob, joint$cells)
  i1 <- match(pt$disease1, ids)
  i2 <- match(pt$disease2, ids)
  p11 <- J[cbind(i1, i2)]
  out <- data.frame(disease1 = pt$disease1, disease2 = pt$disease2,
                    fitted_or = unname(joint_to_or(p11, diag(J)[i1], diag(J)[i2])), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
