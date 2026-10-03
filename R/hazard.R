#' Describe raw culling (or mortality) hazard ratios
#'
#' Hazard ratios are not additive impacts, so they are adjusted by a
#' separate adapter, [deconflate_hr()], outside the additive engine.
#'
#' @param disease Character vector of disease ids (one row per disease; use
#'   1 for a disease with no effect).
#' @param value Positive hazard ratios.
#' @param estimand `"crude"` (unadjusted for other diseases) or `"adjusted"`
#'   (from a model that included the diseases in `adjusted_for`), one per row
#'   or recycled.
#' @param adjusted_for For `"adjusted"`: the diseases the estimate was
#'   adjusted for, separated by `";"`, or `"all"` (every other disease).
#' @param source Optional citation.
#' @return A `cm_hazard_ratios` data frame.
#' @export
#' @examples
#' cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3))
cm_hazard_ratios <- function(disease, value, estimand = "crude", adjusted_for = NA_character_,
                             source = NA_character_) {
  disease <- as.character(disease)
  n <- length(disease)
  if (length(value) != n) cm_abort("`value` must have one entry per disease.")
  check_numeric(value, "value")
  if (any(value <= 0)) cm_abort("Hazard ratios must be positive.")
  if (anyDuplicated(disease)) cm_abort("Each disease may have only one hazard ratio.")
  estimand <- as.character(recycle_arg(estimand, n, "estimand"))
  check_choices(estimand, c("crude", "adjusted"), "estimand")
  adjusted_for <- as.character(recycle_arg(adjusted_for, n, "adjusted_for"))
  has_adj <- vapply(adjusted_for, function(x) length(split_ids(x)) > 0, logical(1))
  if (any(estimand == "crude" & has_adj)) {
    cm_abort("`adjusted_for` is given for crude hazard ratios; set estimand = 'adjusted'.",
             class = "deconflate_unsupported")
  }
  if (any(estimand == "adjusted" & !has_adj)) {
    cm_abort("estimand = 'adjusted' needs `adjusted_for` (disease ids or 'all').")
  }
  out <- data.frame(disease = disease, value = value, estimand = estimand,
                    adjusted_for = adjusted_for,
                    source = as.character(recycle_arg(source, n, "source")),
                    stringsAsFactors = FALSE)
  class(out) <- c("cm_hazard_ratios", "data.frame")
  out
}

#' Combine a population with hazard ratios
#'
#' @param population A [cm_population()] (or [cm_model()], whose impacts are
#'   ignored).
#' @param hazard_ratios A [cm_hazard_ratios()] object with one value per
#'   disease.
#' @return A `cm_hr_model` object.
#' @export
cm_hr_model <- function(population, hazard_ratios) {
  check_population(population)
  if (!inherits(hazard_ratios, "cm_hazard_ratios")) {
    cm_abort("`hazard_ratios` must be created with cm_hazard_ratios().")
  }
  ids <- population$diseases$id
  unk <- setdiff(hazard_ratios$disease, ids)
  if (length(unk)) cm_abort(sprintf("Hazard ratios refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
  miss <- setdiff(ids, hazard_ratios$disease)
  if (length(miss)) {
    cm_abort(sprintf("No hazard ratio for: %s (use 1 for no effect).", paste(miss, collapse = ", ")))
  }
  for (r in seq_len(nrow(hazard_ratios))) {
    s <- split_ids(hazard_ratios$adjusted_for[r])
    if (length(s) == 1L && tolower(s) == "all") next
    unk <- setdiff(s, ids)
    if (length(unk)) cm_abort(sprintf("`adjusted_for` refers to unknown diseases: %s.", paste(unk, collapse = ", ")))
  }
  pop <- population
  pop$impacts <- NULL
  pop$interactions <- NULL
  class(pop) <- "cm_population"
  hr <- hazard_ratios[match(ids, hazard_ratios$disease), , drop = FALSE]
  rownames(hr) <- NULL
  structure(list(population = pop, hazard_ratios = hr), class = "cm_hr_model")
}

#' Adjust hazard ratios for comorbidity (hazard-ratio adapter)
#'
#' Adjusts raw hazard ratios of culling (or mortality) for the hazard ratios
#' of associated diseases. This is a separate adapter from the additive
#' engine ([deconflate()]), because hazard ratios combine multiplicatively.
#'
#' @section Methods:
#' * `"snapshot"` (default): a snapshot hazard-multiplier model. An animal's
#'   hazard is `h0 * exp(sum_i beta_i * D_i)`, and the raw hazard ratio of
#'   disease `i` is taken to be the ratio of the average hazard multiplier
#'   among animals with and without `i`, over the fitted distribution of
#'   disease combinations ([fit_joint()]) at the start of follow-up. The
#'   `beta`s are solved so that these ratios equal the raw hazard ratios. For
#'   an adjusted estimate, the ratio is computed within strata of its
#'   adjustment set and combined across strata with Mantel-Haenszel-type
#'   weights; with `adjusted_for = "all"` the hazard ratio is used as it is.
#' * `"first_order"`: the log-linear approximation,
#'   `log(HR_raw) = A beta`, with `A` as in [deconflate()] (pairwise tables
#'   only).
#' * `"published"`: Rasmussen et al. (2024): `HR - 1` adjusted with eq. 16
#'   (crude estimates only). For reproduction and comparison.
#'
#' @section What the snapshot model is not:
#' A Cox hazard ratio estimated over follow-up is not, in general, the
#' snapshot ratio: animals with high hazards leave first, so the mixture of
#' disease combinations among survivors changes over time, and the marginal
#' hazard ratio changes with it. The snapshot model therefore does not give
#' an exact de-conflation of published Cox coefficients. It is exact for its
#' own estimand (instantaneous marginal ratios at baseline), and a reasonable
#' approximation when follow-up is short relative to the hazards or the
#' diseases are rare. Results are hazard ratios; turn them into culling
#' attributable to disease with [attributable_risk()].
#'
#' @param model A [cm_hr_model()].
#' @param method `"snapshot"`, `"first_order"` or `"published"`.
#' @param joint Optional [fit_joint()] result, checked against the
#'   population. The snapshot method uses it; the other methods keep it in
#'   the result for [attributable_risk()].
#' @param warn Logical: warn when adjusted hazard ratios cross 1 or are not
#'   positive and finite?
#' @param ... Passed to [fit_joint()].
#' @return A `cm_hr_result` with `adjusted` (raw and adjusted hazard ratios),
#'   `diagnostics`, `joint`, `model` and `method`.
#' @export
#' @examples
#' hr <- cm_hr_model(example_supplement(), cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3)))
#' deconflate_hr(hr)$adjusted
#' deconflate_hr(hr, method = "first_order")$adjusted
deconflate_hr <- function(model, method = c("snapshot", "first_order", "published"),
                          joint = NULL, warn = TRUE, ...) {
  if (!inherits(model, "cm_hr_model")) cm_abort("`model` must be created with cm_hr_model().")
  method <- match.arg(method)
  pop <- model$population
  hr <- model$hazard_ratios
  ids <- pop$diseases$id
  n <- length(ids)
  P <- stats::setNames(pop$diseases$prob, ids)
  as_imp <- data.frame(estimand = ifelse(hr$estimand == "adjusted", "adjusted_linear", "crude"),
                       adjusted_for = hr$adjusted_for, stringsAsFactors = FALSE)
  if (method == "published" && any(hr$estimand != "crude")) {
    cm_abort("The published (2024) approach is defined for crude hazard ratios only.",
             class = "deconflate_unsupported")
  }
  if (method == "snapshot") {
    if (is.null(joint)) {
      joint <- withCallingHandlers(fit_joint(pop, ...),
                                   deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
    } else {
      validate_joint(joint, pop)
    }
    if (!joint$converged) {
      cm_abort("The joint distribution did not converge; see check_feasibility(method = 'lp').",
               class = "deconflate_nonconvergence")
    }
    J <- crossprod(joint$cells * joint$prob, joint$cells)
    feas <- sprintf("joint distribution fitted (max residual %.1e)", joint$max_residual)
  } else {
    # A supplied joint is not used here, but it is kept for
    # attributable_risk(), so it must belong to the population.
    if (!is.null(joint)) validate_joint(joint, pop)
    pt <- pair_tables(pop)
    if (anyNA(pt$p11)) {
      cm_abort("Unknown associations need method = 'snapshot'.", class = "deconflate_unsupported")
    }
    tri <- triple_screen(pt, P)
    if (nrow(tri)) cm_abort("The pairwise associations are jointly infeasible; see check_feasibility().",
                            class = "deconflate_infeasible")
    J <- diag(P, nrow = n)
    J[cbind(match(pt$disease1, ids), match(pt$disease2, ids))] <- pt$p11
    J[cbind(match(pt$disease2, ids), match(pt$disease1, ids))] <- pt$p11
    feas <- "triple screen passed (necessary condition only)"
  }
  dimnames(J) <- list(ids, ids)
  A <- projection_matrix(pair_covariance(J), as_imp, ids)
  b_raw <- log(hr$value)

  if (method == "published") {
    m <- hr$value - 1
    conf <- as.vector((A - diag(n)) %*% m)
    adj <- 1 + ifelse(m == 0, 0, m^2 / (m + conf))
    resid <- max(abs(as.vector(A %*% (adj - 1)) - m))
    cond <- kappa(A, exact = TRUE)
  } else {
    if (qr(A)$rank < n) cm_abort("The conflation matrix is singular.", class = "deconflate_singular")
    beta <- as.vector(solve(A, b_raw))
    if (method == "first_order") {
      adj <- exp(beta)
      resid <- max(abs(as.vector(A %*% beta) - b_raw))
      cond <- kappa(A, exact = TRUE)
    } else {
      sol <- solve_snapshot(b_raw, joint$cells, joint$prob, as_imp, ids, start = beta)
      adj <- exp(sol$beta)
      resid <- sol$residual
      cond <- kappa(sol$jacobian, exact = TRUE)
    }
  }
  adj <- unname(adj)
  sign_change <- is.finite(adj) & abs(hr$value - 1) > 1e-12 & abs(adj - 1) > 1e-12 &
    sign(adj - 1) != sign(hr$value - 1)
  if (warn && any(sign_change)) {
    cm_warn(sprintf("Adjusted hazard ratios cross 1 for %s. The raw hazard ratios are smaller than the associated diseases alone would produce; check the estimands and source populations.",
                    paste(ids[sign_change], collapse = ", ")), class = "deconflate_sign_change")
  }
  # The published approximation can divide by a quantity near or below zero.
  bad <- !is.finite(adj) | adj <= 0
  if (warn && any(bad)) {
    cm_warn(sprintf("Adjusted hazard ratios are not positive and finite for %s.",
                    paste(ids[bad], collapse = ", ")), class = "deconflate_nonfinite_warning")
  }
  structure(list(
    adjusted = data.frame(disease = ids, raw = hr$value, adjusted = adj,
                          change = adj / hr$value - 1, estimand = hr$estimand,
                          adjusted_for = hr$adjusted_for, stringsAsFactors = FALSE),
    diagnostics = data.frame(method = method, max_reconstruction_residual = resid,
                             condition_number = cond, n_sign_changes = sum(sign_change),
                             sign_changes = paste(ids[sign_change], collapse = ", "),
                             feasibility = feas, stringsAsFactors = FALSE),
    joint = joint, model = model, method = method
  ), class = "cm_hr_result")
}

# Solve the snapshot hazard-multiplier equations by damped Newton iterations
# with a numerical Jacobian. Equation i: the (stratified) log ratio of mean
# hazard multipliers among animals with and without disease i equals b_raw[i].
solve_snapshot <- function(b_raw, cells, prob, estimands, ids, start, tol = 1e-11,
                           max_iter = 200L) {
  n <- length(b_raw)
  strata <- lapply(seq_len(n), function(i) {
    S <- if (estimands$estimand[i] == "adjusted_linear") {
      match(resolve_adjusted_for(estimands$adjusted_for[i], ids[i], ids), ids)
    } else integer(0)
    if (!length(S)) return(rep(1L, nrow(cells)))
    as.integer(cells[, S, drop = FALSE] %*% 2^(seq_along(S) - 1)) + 1L
  })
  base <- lapply(seq_len(n), function(i) {
    x <- cells[, i]
    d1 <- rowsum(prob * x, strata[[i]], reorder = FALSE)[, 1]
    d0 <- rowsum(prob * (1 - x), strata[[i]], reorder = FALSE)[, 1]
    ok <- d1 > 0 & d0 > 0
    if (!any(ok)) {
      cm_abort(sprintf("The hazard ratio of %s is not identifiable: no stratum of its adjustment set contains animals with and without it.",
                       ids[i]), class = "deconflate_singular")
    }
    list(d1 = d1, d0 = d0, ok = ok, w = ifelse(ok, d1 * d0 / (d1 + d0), 0))
  })
  eval_F <- function(beta) {
    W <- prob * exp(as.vector(cells %*% beta))
    vapply(seq_len(n), function(i) {
      x <- cells[, i]
      bi <- base[[i]]
      n1 <- rowsum(W * x, strata[[i]], reorder = FALSE)[, 1]
      n0 <- rowsum(W * (1 - x), strata[[i]], reorder = FALSE)[, 1]
      r <- log(n1[bi$ok] / bi$d1[bi$ok]) - log(n0[bi$ok] / bi$d0[bi$ok])
      sum(bi$w[bi$ok] * r) / sum(bi$w[bi$ok]) - b_raw[i]
    }, numeric(1))
  }
  jac <- function(beta, f0) {
    h <- 1e-6
    vapply(seq_len(n), function(j) {
      e <- rep(0, n)
      e[j] <- h
      (eval_F(beta + e) - eval_F(beta - e)) / (2 * h)
    }, numeric(n))
  }
  size <- function(f) if (all(is.finite(f))) sum(f^2) else Inf
  beta <- start
  f <- eval_F(beta)
  for (it in seq_len(max_iter)) {
    if (max(abs(f)) < tol) break
    Jm <- jac(beta, f)
    step <- tryCatch(solve(Jm, f), error = function(e) NULL)
    if (is.null(step)) break
    lambda <- 1
    improved <- FALSE
    while (lambda >= 1e-8) {
      cand <- beta - lambda * step
      fc <- eval_F(cand)
      if (size(fc) < size(f)) {
        improved <- TRUE
        break
      }
      lambda <- lambda / 2
    }
    if (!improved) break
    beta <- cand
    f <- fc
  }
  resid <- max(abs(f))
  if (!is.finite(resid) || resid > 1e-8) {
    cm_abort(sprintf("Could not solve the snapshot hazard model (max residual %.2e).", resid),
             class = "deconflate_nonconvergence")
  }
  list(beta = unname(beta), residual = resid, jacobian = jac(beta, f))
}

#' Culling (or mortality) attributable to disease
#'
#' Converts adjusted hazard ratios into the part of an event's overall risk
#' (e.g. annual culling) that is attributable to disease, without counting
#' an animal with several diseases more than once, and allocates it to
#' diseases.
#'
#' Within the period, an animal with disease combination `x` has a constant
#' hazard `h0 * exp(sum_i beta[i] * x[i])`, with `beta = log(adjusted HR)`.
#' The baseline hazard `h0` is chosen so that the population risk, averaged
#' over the distribution of disease combinations ([fit_joint()]), equals
#' `overall_risk`. The disease-free risk is `1 - exp(-h0)`, and the
#' attributable risk is `overall_risk - (1 - exp(-h0))`. An animal's risk
#' cannot exceed 1, so the attributable risk is smaller than the sum of
#' per-disease excess risks when diseases co-occur.
#'
#' The attributable risk is allocated to diseases by Shapley values over all
#' disease combinations ([shapley_by_cell()]). Any combinations skipped with
#' `max_present` are reported as unallocated.
#'
#' This uses the same snapshot model as [deconflate_hr()] (constant hazards
#' within the period, multiplicative hazard ratios, no change in the mixture
#' of diseases over the period).
#'
#' @param result A [deconflate_hr()] result.
#' @param overall_risk Overall period risk of the event, as a proportion
#'   (e.g. `0.25` for an annual culling rate of 25%).
#' @param unit_value Optional value per animal removed (e.g. replacement
#'   cost less salvage value); adds `value` columns.
#' @param joint Optional [fit_joint()] result; by default the result's own
#'   joint distribution or a new fit. It is checked against the population.
#' @param allocate Logical: allocate the attributable risk to diseases?
#' @param max_present Passed to [shapley_by_cell()] (default: no skipping).
#' @return A `cm_attributable` list with `summary` (overall, disease-free and
#'   attributable risk, attributable fraction, value, and any unallocated
#'   part), `by_disease` and `baseline_hazard`.
#' @export
#' @examples
#' hr <- cm_hr_model(example_supplement(), cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3)))
#' attributable_risk(deconflate_hr(hr), overall_risk = 0.25, unit_value = 1300)
attributable_risk <- function(result, overall_risk, unit_value = NULL, joint = NULL,
                              allocate = TRUE, max_present = Inf) {
  if (!inherits(result, "cm_hr_result")) cm_abort("`result` must come from deconflate_hr().")
  check_numeric(overall_risk, "overall_risk")
  if (length(overall_risk) != 1L || overall_risk <= 0 || overall_risk >= 1) {
    cm_abort("`overall_risk` must be a single proportion between 0 and 1.")
  }
  a <- result$adjusted
  bad <- !is.finite(a$adjusted) | a$adjusted <= 0
  if (any(bad)) {
    cm_abort(sprintf("Adjusted hazard ratios must be positive and finite (check: %s).",
                     paste(a$disease[bad], collapse = ", ")), class = "deconflate_nonfinite")
  }
  pop <- result$model$population
  if (!is.null(joint)) {
    validate_joint(joint, pop)
  } else {
    joint <- result$joint %||% fit_joint(pop)
  }
  if (!isTRUE(joint$converged)) {
    cm_abort("The joint distribution did not converge.", class = "deconflate_nonconvergence")
  }
  ids <- joint$diseases
  beta <- stats::setNames(log(a$adjusted), a$disease)[ids]
  rel <- exp(as.vector(joint$cells %*% beta))
  risk_at <- function(h) sum(joint$prob * (1 - exp(-h * rel)))
  h0 <- stats::uniroot(function(h) risk_at(h) - overall_risk,
                       c(0, -log(1 - overall_risk)), extendInt = "upX", tol = 1e-14)$root
  r0 <- 1 - exp(-h0)
  summ <- data.frame(overall_risk = overall_risk, disease_free_risk = r0,
                     attributable = overall_risk - r0,
                     attributable_fraction = (overall_risk - r0) / overall_risk,
                     unallocated = NA_real_, stringsAsFactors = FALSE)
  if (!is.null(unit_value)) summ$value <- summ$attributable * unit_value
  by <- NULL
  if (allocate) {
    loss <- function(x) 1 - exp(-h0 * exp(sum(beta * x[ids]))) - r0
    sh <- withCallingHandlers(shapley_by_cell(joint, loss, max_present = max_present),
                              deconflate_incomplete_allocation = function(w) invokeRestart("muffleWarning"))
    by <- data.frame(disease = ids, hr_adjusted = unname(exp(beta)),
                     attributable = sh$shapley, share = sh$share, stringsAsFactors = FALSE)
    if (!is.null(unit_value)) by$value <- by$attributable * unit_value
    summ$unallocated <- summ$attributable - attr(sh, "allocated")
    if (attr(sh, "skipped_mass") > 0) {
      cm_warn(sprintf("%.3g of the attributable risk is not allocated (combinations with more than %g diseases were skipped).",
                      summ$unallocated, max_present), class = "deconflate_incomplete_allocation")
    }
  }
  structure(list(summary = summ, by_disease = by, baseline_hazard = h0),
            class = "cm_attributable")
}

#' @export
print.cm_attributable <- function(x, ...) {
  s <- x$summary
  cat("<cm_attributable> snapshot hazard-multiplier model\n")
  cat(sprintf("  Overall risk %.4g; disease-free risk %.4g; attributable %.4g (%.1f%% of the overall risk)\n",
              s$overall_risk, s$disease_free_risk, s$attributable, 100 * s$attributable_fraction))
  if (!is.null(s$value)) cat(sprintf("  Value: %.4g\n", s$value))
  if (!is.null(x$by_disease)) {
    cat("\n")
    print(x$by_disease, row.names = FALSE, digits = 4)
    if (isTRUE(abs(s$unallocated) > 1e-10)) cat(sprintf("\nUnallocated: %.3g\n", s$unallocated))
  }
  invisible(x)
}

#' @export
print.cm_hr_model <- function(x, ...) {
  cat("<cm_hr_model>\n")
  print(x$population)
  cat(sprintf("  Hazard ratios: %d (%d adjusted)\n", nrow(x$hazard_ratios),
              sum(x$hazard_ratios$estimand == "adjusted")))
  invisible(x)
}

#' @export
print.cm_hr_result <- function(x, ...) {
  cat(sprintf("<cm_hr_result> method: %s\n\n", x$method))
  print(x$adjusted[, c("disease", "raw", "adjusted", "change", "estimand")], row.names = FALSE, digits = 4)
  cat("\nDiagnostics:\n")
  print(x$diagnostics[, c("max_reconstruction_residual", "n_sign_changes", "condition_number", "feasibility")],
        row.names = FALSE, digits = 3)
  invisible(x)
}
