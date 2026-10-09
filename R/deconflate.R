#' Adjust raw impact estimates for comorbidity
#'
#' Raw impact estimates (comparisons of animals with and without a disease)
#' are treated as conflations of the disease's own impact and the impacts of
#' associated diseases. `deconflate()` removes the part of each estimate that
#' belongs to other diseases, for one impact table at a time: additive
#' impacts (the default), or event impacts (`event_model = TRUE`). When
#' inputs have distributions, it also runs draws (`n_draws`) and reports
#' intervals.
#'
#' @section Additive impacts:
#' For one vector of additive impacts `b` (in any units), the raw estimates
#' satisfy
#'
#' `raw = A %*% b + offset`,
#'
#' where `A` depends on the disease probabilities and associations (and on
#' each estimate's estimand), and `offset` holds the contribution of any
#' pairwise interactions. Results are in the units of the impacts supplied.
#'
#' For a crude estimate of disease `i`, `A[i, k] = P(k | i) - P(k | not i)`
#' (the excess probability; Rasmussen et al. 2022, eq. 14). For a coefficient
#' from an additive regression adjusted for the diseases `S`
#' (`estimand = "adjusted_linear"` in [cm_impacts()]), `A[i, k]` is the
#' coefficient of `D_i` in the population linear projection of `D_k` on
#' `D_i` and `D_S`, and `A[i, S] = 0`. These coefficients depend only on the
#' pairwise tables. They assume the probabilities and associations describe
#' the population the estimates come from. An adjustment set that is
#' collinear with the disease is not identifiable and is rejected.
#'
#' With interactions `delta`, the offset of disease `i` is the corresponding
#' regression coefficient of the interaction burden
#' `sum_{j<k} delta[j, k] D_j D_k`, which needs the joint distribution of
#' disease combinations.
#'
#' The expected aggregate impact per animal is
#' `sum_i p_i b_i + sum_{j<k} P(j and k) delta[j, k]`. Each disease's
#' contribution is its own term plus half of each interaction term it is
#' involved in (the closed-form Shapley value); contributions add up to the
#' aggregate. Shares are `NA` when the aggregate is zero.
#'
#' @section Methods:
#' * The simultaneous method solves `raw = A b` exactly from the pairwise
#'   2x2 tables. No joint distribution is needed, so it works for any number
#'   of diseases.
#' * The global method fits the maximum-entropy distribution of disease
#'   combinations ([fit_joint()]; the "iterative" model) and solves the same
#'   equations, including any interactions. Pairs without an association are
#'   unknown: the fit fills them in from the other associations.
#'
#' The two give the same results when every pair has an association, there
#' are no interactions and no three-way terms. `method = "auto"` (the
#' default) uses the simultaneous method then, and the global method
#' otherwise, with a note giving the reason; `"simultaneous"` behaves the
#' same way (the global method is used when it is needed), and `"global"`
#' always uses the global method. The global method's exact backend
#' enumerates 2^n combinations; with more than 20 diseases the sampled
#' backend is used (see [fit_joint()]).
#'
#' The proportional approximation of Rasmussen et al. (2022, eq. 16) is not a
#' method here: it is kept for comparison and reproduction in
#' [compare_methods()] and the `reproduce_*()` functions.
#'
#' @section Event impacts:
#' With `event_model = TRUE`, the impact table holds event impacts (hazard
#' ratios, rate ratios, risk ratios, odds ratios or risk differences; see
#' [cm_impacts()]), which are adjusted with the snapshot hazard model. Within
#' the period of the overall risk, an animal with disease combination `d` has
#' a constant hazard `h0 * exp(sum_i beta_i d_i)`, so its risk over the
#' period is `R(d) = 1 - exp(-h0 * exp(sum_i beta_i d_i))`; `h0` is set so
#' that the population risk, averaged over the joint distribution of disease
#' combinations, equals `overall_risk`. The `beta`s are solved so that each
#' raw estimate is reproduced:
#' * a hazard ratio or rate ratio by the ratio of the average hazard
#'   multipliers among animals with and without the disease;
#' * a risk ratio, odds ratio or risk difference by the ratio, odds ratio or
#'   difference of the average risks `R` among animals with and without the
#'   disease;
#'
#' at the start of the period (crude), or within strata of the estimate's
#' adjustment set, combined with Mantel-Haenszel-type weights (stratified).
#' The adjusted hazard ratios `exp(beta)` are each disease's own hazard
#' multiplier. The risk attributable to disease is `overall_risk` minus the
#' disease-free risk `1 - exp(-h0)`; an animal's risk cannot exceed 1, so it
#' is smaller than the sum of per-disease excess risks when diseases
#' co-occur. It is allocated to diseases by Shapley values over disease
#' combinations ([shapley_by_cell()]).
#'
#' What the snapshot model is not: a Cox hazard ratio estimated over
#' follow-up is not, in general, the snapshot ratio, because animals with
#' high hazards leave first and the mixture of disease combinations among
#' survivors changes. The model is exact for its own estimands and a
#' reasonable approximation when follow-up is short relative to the hazards
#' or the diseases are rare. Risk-based estimates (risk ratios, odds ratios,
#' risk differences) must refer to the same period as `overall_risk`.
#'
#' @section Uncertainty:
#' Inputs with a distribution (the `dist` columns of the input tables, or
#' `distributions` in [cm_model()]; and `overall_risk` given as a
#' distribution) are drawn `n_draws` times; each draw is adjusted in the same
#' way as the central estimate. The central estimate uses the point values
#' (`value`, and the mean of an `overall_risk` distribution); the draws give
#' 95% intervals (2.5% and 97.5% quantiles), means, Monte Carlo standard
#' errors and stability checks (`$draws$summary`). Draws whose values cannot
#' hold together (e.g. associations that no population can have at once) are
#' rejected and counted (`$draws$rejections`), not replaced; notes in
#' `$notes` report a high rejection share and limited Monte Carlo precision.
#' Without distributions, no draws are run.
#'
#' @section Feasibility:
#' Pairwise tables can each be valid while no population has all of them
#' (an invertible `A` does not mean the inputs are feasible). The
#' simultaneous method screens every triple of diseases (a necessary
#' condition; see [check_feasibility()]); `feasibility = "lp"` runs the exact
#' check. The global method fits the joint distribution, and stops if it
#' does not converge.
#'
#' @param model A [cm_model()] (e.g. from [cm_read_inputs()]) with at least
#'   one association.
#' @param method `"auto"` (default), `"simultaneous"` or `"global"` (see
#'   Methods). Not used for event impacts, which always use the snapshot
#'   model.
#' @param event_model `FALSE` (default) for additive impacts; `TRUE` for
#'   event impacts (an impact table with a `measure` column).
#' @param overall_risk With `event_model = TRUE` (required): the overall risk
#'   of the event in the population over the period, as a proportion (e.g.
#'   `0.25`), or a distribution (e.g. `dist_beta(250, 750)`).
#' @param n_draws Number of draws for the uncertainty (default 1000; 0 for
#'   point estimates only; otherwise at least 2). Used only when some input has a distribution.
#' @param seed Optional random seed for the draws and for a sampled joint
#'   distribution (one is chosen and stored
#'   in the result otherwise).
#' @param sampling `"random"` (default) or `"lhs"` (Latin hypercube, in
#'   `lhs_replicates` independent blocks; the Monte Carlo error is then
#'   estimated from the block means).
#' @param lhs_replicates Number of Latin hypercube blocks (at least 2; at
#'   most `n_draws / 2` are used).
#' @param joint Optional [fit_joint()] result for the global method or the
#'   event model. It is checked against the model (diseases, probabilities,
#'   associations and three-way terms).
#' @param warn Logical: warn when adjusted impacts change sign or are not
#'   finite?
#' @param feasibility For the simultaneous method: `"screen"` (default;
#'   triple screen), `"lp"` (exact check, needs `lpSolve`) or `"none"`.
#' @param ... Passed to [fit_joint()] (e.g. `backend = "sampled"`).
#'
#' @return For additive impacts, a `cm_result` with elements:
#'   * `adjusted`: raw and adjusted impacts, relative change and estimand
#'     (with draws, also `lower` and `upper`);
#'   * `totals`: the sum of `p_i * raw_i` (the naive aggregate), the adjusted
#'     aggregate and its interaction part;
#'   * `contributions`: per disease, main and interaction contributions,
#'     their total and share;
#'   * `diagnostics`: reconstruction residual, rank and condition number of
#'     `A`, sign changes and the feasibility check;
#'   * `unknown_pairs`: pairs without an association and the odds ratios the
#'     global fit gave them;
#'   * `draws`: the uncertainty (summary, accepted draws, rejections, seed),
#'     or `NULL`;
#'   * `notes`, `conflation` (`A` and `offset`), `interactions` (matrix of
#'     `delta`), `joint_pairs` (matrix of `P(j and k)`), `model`, `joint`,
#'     `method`, `label` and `units`.
#'
#'   For event impacts, a `cm_event_result` with `adjusted` (raw estimates
#'   and adjusted hazard ratios), `attributable` (`summary`: overall,
#'   disease-free and attributable risk and the attributable fraction;
#'   `by_disease`: the Shapley allocation; `baseline_hazard`), `diagnostics`,
#'   `unknown_pairs`, `draws`, `notes`, `joint`, `model`, `method` and
#'   `overall_risk`. Risks are proportions of animals with the event during
#'   the period.
#' @export
#' @examples
#' res <- deconflate(example_supplement())
#' res$adjusted
#' res$totals
#'
#' # A coefficient adjusted for one other disease
#' m <- example_supplement()
#' m$impacts <- cm_impacts(c("d1", "d2", "d3"), c(2.2, 5, 7.5),
#'                         estimand = c("adjusted_linear", "crude", "crude"),
#'                         adjusted_for = c("d2", NA, NA), units = "%")
#' deconflate(m)$adjusted
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
#'                             estimand = "snapshot_crude", label = "culling"))
#' deconflate(cull, event_model = TRUE, overall_risk = 0.25)
deconflate <- function(model, ...) UseMethod("deconflate")

#' @export
deconflate.default <- function(model, ...) {
  cm_abort("`model` must be created with cm_model() or cm_read_inputs().")
}

#' @rdname deconflate
#' @export
deconflate.cm_model <- function(model, method = c("auto", "simultaneous", "global"),
                                event_model = FALSE, overall_risk = NULL, n_draws = 1000,
                                seed = NULL, sampling = c("random", "lhs"), lhs_replicates = 10L,
                                joint = NULL, warn = TRUE,
                                feasibility = c("screen", "lp", "none"), ...) {
  given <- !missing(method)
  method <- public_method(method, c("auto", "simultaneous", "global"), "deconflate()")
  sampling <- match.arg(sampling)
  feasibility <- match.arg(feasibility)
  check_draw_args(n_draws, seed, sampling, lhs_replicates)
  n_draws <- as.integer(n_draws)
  check_event_model(model, event_model)
  risk <- check_overall_risk(overall_risk, event_model)
  a <- model$associations
  if (is.null(a) || !nrow(a)) {
    cm_abort("No association estimates were given, so there is nothing to de-conflate (the diseases would be treated as independent). To see how much associations could change the results, use screen_associations() or cm_threshold().",
             class = "deconflate_unsupported")
  }
  dots <- list(...)
  plan <- plan_method(model, method, event_model, given, dots, has_joint = !is.null(joint))
  # The seed also makes a sampled joint distribution reproducible.
  if (!is.null(seed) && identical(plan$fit_args$backend %||% dots$backend, "sampled") && is.null(dots$seed)) {
    plan$fit_args$seed <- seed
  }
  res <- run_point(model, plan, risk = risk$value, joint = joint, warn = warn,
                   feasibility = feasibility, dots = dots)
  res$notes <- plan$notes
  res$unknown_pairs <- unknown_pairs_table(model, res$joint)
  specs <- draw_specs(model, risk$dist)
  if (n_draws > 0L) {
    if (!length(specs)) {
      res$notes <- c(res$notes, "No input has a distribution, so no draws were run: the results are point estimates.")
    } else {
      res$draws <- deconflate_draws(model, specs, plan, res, risk$value, n_draws, seed, sampling,
                                    lhs_replicates, feasibility, dots)
      res <- attach_intervals(res)
      res$notes <- c(res$notes, draw_notes(res$draws))
    }
  }
  res
}

# Additive or event impacts, as the call says? A mismatch is an error.
check_event_model <- function(model, event_model) {
  if (!is.logical(event_model) || length(event_model) != 1L || is.na(event_model)) {
    cm_abort("`event_model` must be TRUE or FALSE.")
  }
  if (is.null(model$impacts)) cm_abort("The model has no impacts to adjust.")
  kind <- impact_kind(model$impacts)
  if (event_model && kind != "event") {
    cm_abort("event_model = TRUE needs event impacts: an impact table with a measure column (HR, rate_ratio, RR, OR or RD) and snapshot estimands (see ?cm_impacts).",
             class = "deconflate_unsupported")
  }
  if (!event_model && kind == "event") {
    cm_abort("These are event impacts (the impact table has a measure column): use event_model = TRUE, with overall_risk.",
             class = "deconflate_unsupported")
  }
  invisible(TRUE)
}

# Choose how to adjust: the method actually used, the arguments for the
# joint fit, and notes saying why.
plan_method <- function(model, method, event_model, given = TRUE, dots = list(), has_joint = FALSE) {
  n <- nrow(model$diseases)
  notes <- character(0)
  fit_args <- list()
  sampled_note <- if (has_joint) character(0) else
    sprintf("With %d diseases, the joint distribution was fitted with the sampled backend (see ?fit_joint).", n)
  if (n > 20L && is.null(dots$backend)) {
    fit_args$backend <- "sampled"
  }
  if (event_model) {
    if (given && method == "simultaneous") {
      notes <- c(notes, "Event impacts are always adjusted with the snapshot model, which uses the joint distribution of disease combinations; `method` was not used.")
    }
    if (!is.null(fit_args$backend)) notes <- c(notes, sampled_note)
    return(list(method = "snapshot", event = TRUE, notes = notes, fit_args = fit_args))
  }
  reasons <- global_reasons(model)
  use <- if (method == "global" || length(reasons)) "global" else "simultaneous"
  if (length(reasons) && method != "global") {
    notes <- c(notes, sprintf("%s because of %s.",
                              if (method == "simultaneous") "The global method was used instead of the simultaneous method" else "The global method was used",
                              join_and(reasons)))
  }
  if (use == "global" && !is.null(fit_args$backend)) notes <- c(notes, sampled_note)
  list(method = use, event = FALSE, notes = notes, fit_args = fit_args)
}

# "a", "a and b", "a, b and c".
join_and <- function(x) {
  if (length(x) <= 1L) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

# What needs the global method: interactions, three-way terms, unknown pairs.
global_reasons <- function(model) {
  out <- character(0)
  if (!is.null(model$interactions) && nrow(model$interactions)) out <- c(out, "interactions")
  if (!is.null(model$three_way) && nrow(model$three_way)) out <- c(out, "three-way terms")
  pt <- pair_tables(model)
  k <- sum(pt$status == "unknown")
  if (k) out <- c(out, sprintf("%d pair%s without an association (unknown)", k, if (k == 1L) "" else "s"))
  out
}

# One adjustment of a model, as planned.
run_point <- function(model, plan, risk = NULL, joint = NULL, warn = TRUE,
                      feasibility = "screen", dots = list(), joint_start = NULL, start = NULL) {
  fj <- utils::modifyList(dots[intersect(names(dots), joint_arg_names)], plan$fit_args)
  if (isTRUE(plan$event)) {
    do.call(adjust_event, c(list(model, method = "snapshot", overall_risk = risk, joint = joint,
                                 warn = warn, start = start, joint_start = joint_start), fj))
  } else {
    do.call(adjust_impacts, c(list(model, method = plan$method, joint = joint, warn = warn,
                                   feasibility = feasibility, joint_start = joint_start), fj))
  }
}

# Pairs without an association and the odds ratios the global fit gave them.
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

# The adjustment itself, including the published approximation (eq. 16 of
# Rasmussen et al. 2022), which is kept for compare_methods() and the
# reproduce_*() functions only.
adjust_impacts <- function(model, method = c("simultaneous", "published", "global"),
                           joint = NULL, warn = TRUE,
                           feasibility = c("screen", "lp", "none"), joint_start = NULL, ...) {
  method <- match.arg(method)
  feasibility <- match.arg(feasibility)
  imp <- model$impacts
  if (is.null(imp)) cm_abort("The model has no impacts to adjust.")
  if (impact_kind(imp) == "event") {
    cm_abort("These are event impacts (the impact table has a measure column): use event_model = TRUE, with overall_risk.",
             class = "deconflate_unsupported")
  }
  ids <- model$diseases$id
  if (!identical(imp$disease, ids)) {
    # Impacts replaced after cm_model() (e.g. `m$impacts <- ...`): validate
    # them and put them in disease order.
    model <- cm_model(model, imp, model$interactions)
    imp <- model$impacts
  }
  n <- length(ids)
  P <- stats::setNames(model$diseases$prob, ids)
  m_raw <- imp$value
  has_int <- !is.null(model$interactions) && nrow(model$interactions) > 0L
  if (method == "published" && any(imp$estimand != "crude")) {
    cm_abort("The published approximation is defined for crude estimates only; use method = 'simultaneous' or 'global' for adjusted estimands.",
             class = "deconflate_unsupported")
  }
  if (has_int && method != "global") {
    cm_abort("Interactions need the global method (they need probabilities of disease combinations).",
             class = "deconflate_unsupported")
  }

  cells <- NULL
  if (method == "global") {
    if (is.null(joint)) {
      joint <- withCallingHandlers(fit_joint(as_population(model), start = joint_start, ...),
                                   deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
    } else {
      validate_joint(joint, model)
    }
    if (!joint$converged) {
      cm_abort(if (identical(joint$backend, "sampled")) {
        sprintf("The sampled joint distribution is unresolved (residual %.2e); see its diagnostics, and try more samples or chains. This does not show that the inputs are infeasible.",
                joint$max_residual)
      } else {
        sprintf("The joint distribution did not converge (max residual %.2e). The pairwise associations may be jointly infeasible; see check_feasibility(method = 'lp').",
                joint$max_residual)
      }, class = "deconflate_nonconvergence")
    }
    cells <- joint$cells
    J <- crossprod(cells * joint$prob, cells)
    dimnames(J) <- list(ids, ids)
    feas <- sprintf("joint distribution fitted (max residual %.1e)", joint$max_residual)
  } else {
    pt <- pair_tables(model)
    if (anyNA(pt$p11)) {
      unk <- paste(pt$disease1, pt$disease2, sep = "-")[is.na(pt$p11)]
      cm_abort(sprintf(
        "No association for %s: the pairwise methods need every pair. The global method fills unknown pairs in (deconflate() uses it automatically); to treat a pair as unrelated, give it an odds ratio of 1.",
        paste(unk, collapse = ", ")), class = c("deconflate_unknown_pairs", "deconflate_unsupported"))
    }
    J <- diag(P, nrow = n)
    dimnames(J) <- list(ids, ids)
    J[cbind(pt$disease1, pt$disease2)] <- pt$p11
    J[cbind(pt$disease2, pt$disease1)] <- pt$p11
    feas <- switch(feasibility,
      screen = {
        tri <- triple_screen(pt, P)
        if (nrow(tri)) {
          cm_abort(sprintf("The pairwise associations are jointly infeasible: no population has all of them (conflicting triples: %s). See check_feasibility().",
                           paste(utils::head(paste(tri$disease1, tri$disease2, tri$disease3, sep = "-"), 5),
                                 collapse = ", ")), class = "deconflate_infeasible")
        }
        "triple screen passed (necessary condition only)"
      },
      lp = {
        cf <- check_feasibility(model, method = "lp")
        if (!isTRUE(cf$feasible)) {
          cm_abort("The pairwise associations are jointly infeasible (LP check); see check_feasibility().",
                   class = "deconflate_infeasible")
        }
        "LP check passed (jointly feasible)"
      },
      none = "unchecked"
    )
  }

  A <- projection_matrix(pair_covariance(J), imp, ids)
  offset <- rep(0, n)
  D <- matrix(0, n, n, dimnames = list(ids, ids))
  if (has_int) {
    int <- model$interactions
    D[cbind(int$disease1, int$disease2)] <- int$value
    D[cbind(int$disease2, int$disease1)] <- int$value
    # Interaction burden of each combination: sum_{j<k} delta_jk D_j D_k.
    burden <- rowSums((cells %*% D) * cells) / 2
    offset <- projection_offset(cells, joint$prob, burden, imp, ids)
  }

  if (method == "published") {
    conf <- as.vector((A - diag(n)) %*% m_raw)
    m_adj <- ifelse(m_raw == 0, 0, m_raw^2 / (m_raw + conf))
    rk <- qr(A)$rank
  } else {
    rk <- qr(A)$rank
    if (rk < n) {
      cm_abort(sprintf("The conflation matrix is singular (rank %d of %d): the impacts are not identifiable from these inputs.",
                       rk, n), class = "deconflate_singular")
    }
    m_adj <- tryCatch(as.vector(solve(A, m_raw - offset)),
                      error = function(e) cm_abort(sprintf("Could not solve the adjustment equations: %s", conditionMessage(e)),
                                                   class = "deconflate_singular"))
  }
  m_adj <- unname(m_adj)
  recon <- as.vector(A %*% m_adj) + offset
  resid <- if (all(is.finite(recon))) max(abs(recon - m_raw)) else NA_real_
  sign_change <- is.finite(m_adj) & abs(m_raw) > 1e-12 & abs(m_adj) > 1e-12 &
    sign(m_adj) != sign(m_raw)

  if (warn && any(sign_change)) {
    cm_warn(sprintf(
      "Adjusted impacts change sign for %s. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.",
      paste(ids[sign_change], collapse = ", ")), class = "deconflate_sign_change")
  }
  if (warn && any(!is.finite(m_adj))) {
    cm_warn(sprintf("Non-finite adjusted impacts for %s.", paste(ids[!is.finite(m_adj)], collapse = ", ")),
            class = "deconflate_nonfinite_warning")
  }

  main <- unname(P * m_adj)
  inter <- unname(0.5 * rowSums(D * J))
  contrib <- main + inter
  total <- sum(contrib)
  share <- if (is.finite(total) && abs(total) > 0) contrib / total else rep(NA_real_, n)

  structure(list(
    adjusted = data.frame(disease = ids, raw = m_raw, adjusted = m_adj,
                          change = ifelse(m_raw != 0, m_adj / m_raw - 1, NA_real_),
                          estimand = imp$estimand, adjusted_for = imp$adjusted_for,
                          stringsAsFactors = FALSE),
    totals = data.frame(raw_sum = sum(P * m_raw), adjusted_total = total,
                        interaction_total = sum(inter)),
    contributions = data.frame(disease = ids, main = main, interaction = inter,
                               total = contrib, share = share, stringsAsFactors = FALSE),
    diagnostics = data.frame(method = method, max_reconstruction_residual = resid,
                             rank = rk, condition_number = kappa(A, exact = TRUE),
                             n_sign_changes = sum(sign_change),
                             sign_changes = paste(ids[sign_change], collapse = ", "),
                             feasibility = feas, stringsAsFactors = FALSE),
    conflation = list(A = A, offset = stats::setNames(offset, ids)),
    interactions = D, joint_pairs = J, model = model, joint = joint, method = method,
    label = attr(imp, "label"), units = attr(imp, "units")
  ), class = "cm_result")
}

# Regression coefficient of D_i for the interaction burden g, given the
# estimand of each impact (crude: difference in means; adjusted_linear: the
# coefficient in a population regression on D_i and D_S).
projection_offset <- function(cells, prob, g, impacts, ids) {
  n <- length(ids)
  p <- colSums(cells * prob)
  J <- crossprod(cells * prob, cells)
  Sigma <- J - outer(p, p)
  diag(Sigma) <- p * (1 - p)
  cov_g <- as.vector(crossprod(cells, prob * g)) - p * sum(prob * g)
  vapply(seq_len(n), function(i) {
    S <- if (impacts$estimand[i] == "adjusted_linear") {
      match(resolve_adjusted_for(impacts$adjusted_for[i], ids[i], ids), ids)
    } else integer(0)
    X <- c(i, S)
    tryCatch(solve(Sigma[X, X, drop = FALSE], cov_g[X])[1],
             error = function(e) cm_abort(sprintf("The interaction offset of %s is not identifiable: %s",
                                                  ids[i], conditionMessage(e)),
                                          class = "deconflate_singular"))
  }, numeric(1))
}

#' Attribute the aggregate to diseases (Shapley allocation)
#'
#' The expected aggregate impact is
#' `sum_i p_i b_i + sum_{j < k} delta[j, k] P(j and k)`. Removing any disease
#' involved in a term removes that term, so each disease's Shapley value is
#' its own term plus an equal share of every interaction term it is involved
#' in. The shares add up to the aggregate. This is the `contributions` table
#' of [deconflate()].
#'
#' For event impacts, it is the Shapley allocation of the risk attributable
#' to disease (`$attributable$by_disease`).
#'
#' @param result A [deconflate()] result.
#' @return A data frame with, per disease, the main contribution, the share
#'   of interaction terms, the total and the share of the aggregate (`NA`
#'   when the aggregate is zero); for event impacts, the attributable risk
#'   and its share.
#' @export
attribute_burden <- function(result) {
  if (inherits(result, "cm_event_result")) {
    if (is.null(result$attributable$by_disease)) cm_abort("The result has no allocation of the attributable risk.")
    return(result$attributable$by_disease)
  }
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  result$contributions
}
