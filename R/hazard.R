# Event impacts: the snapshot hazard model -------------------------------------
#
# Event impacts (cm_impacts() with a `measure`) compare the risk of an event
# (e.g. death or culling) between animals with and without a disease. They
# are not additive, so deconflate(..., event_model = TRUE) adjusts them with
# a hazard model: within the period, an animal with disease combination d
# has a constant hazard h0 * exp(beta' d), so its period risk is
# R(d) = 1 - exp(-h0 exp(beta' d)). The overall risk r fixes h0:
# sum_d p(d) R(d) = r. The adjusted hazard ratios are exp(beta), each
# disease's own hazard multiplier.

# Adjust event impacts. `method`: "snapshot" (exact, any measure),
# "first_order" (log-linear approximation, hazard and rate ratios only) or
# "published" (HR - 1 with eq. 16 of Rasmussen et al. 2022; crude hazard
# ratios only, kept for compare_methods()). `overall_risk`: a proportion.
# `start`: optional starting values of beta (e.g. from the central fit).
adjust_event <- function(model, method = c("snapshot", "first_order", "published"),
                         overall_risk, joint = NULL, warn = TRUE, allocate = TRUE,
                         start = NULL, joint_start = NULL, ...) {
  method <- match.arg(method)
  imp <- model$impacts
  if (is.null(imp)) cm_abort("The model has no impacts to adjust.")
  if (impact_kind(imp) != "event") {
    cm_abort("event_model = TRUE needs event impacts: an impact table with a measure column (HR, rate_ratio, RR, OR or RD) and snapshot estimands.",
             class = "deconflate_unsupported")
  }
  if (!is.null(model$interactions) && nrow(model$interactions)) {
    cm_abort("Interactions apply to additive impacts only, not to event impacts.",
             class = "deconflate_unsupported")
  }
  ids <- model$diseases$id
  if (!identical(imp$disease, ids)) {
    model <- cm_model(model, imp, NULL)
    imp <- model$impacts
  }
  n <- length(ids)
  P <- stats::setNames(model$diseases$prob, ids)
  pop <- as_population(model)
  ratio_m <- c("HR", "rate_ratio")
  if (method != "snapshot" && !all(imp$measure %in% ratio_m)) {
    cm_abort(sprintf("The %s method needs hazard or rate ratios; risk ratios, odds ratios and risk differences need the snapshot model.",
                     if (method == "published") "published" else "first-order"),
             class = "deconflate_unsupported")
  }
  if (method == "published" && any(imp$estimand != "snapshot_crude")) {
    cm_abort("The published (2024) approach is defined for crude hazard ratios only.",
             class = "deconflate_unsupported")
  }
  as_imp <- data.frame(estimand = ifelse(imp$estimand == "snapshot_stratified", "adjusted_linear", "crude"),
                       adjusted_for = imp$adjusted_for, stringsAsFactors = FALSE)

  if (method == "snapshot") {
    joint <- event_joint(pop, joint, joint_start, ...)
    J <- crossprod(joint$cells * joint$prob, joint$cells)
    feas <- sprintf("joint distribution fitted (max residual %.1e)", joint$max_residual)
  } else {
    pt <- pair_tables(pop)
    if (anyNA(pt$p11)) {
      cm_abort(sprintf("The %s method needs an association for every pair (unknown: %s); the snapshot model fills unknown pairs in.",
                       if (method == "published") "published" else "first-order",
                       paste(utils::head(paste(pt$disease1, pt$disease2, sep = "-")[is.na(pt$p11)], 5),
                             collapse = ", ")),
               class = c("deconflate_unknown_pairs", "deconflate_unsupported"))
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

  v <- imp$value
  if (method == "published") {
    m <- v - 1
    conf <- as.vector((A - diag(n)) %*% m)
    adj <- 1 + ifelse(m == 0, 0, m^2 / (m + conf))
    resid <- max(abs(as.vector(A %*% (adj - 1)) - m))
    cond <- kappa(A, exact = TRUE)
  } else {
    if (qr(A)$rank < n) cm_abort("The conflation matrix is singular.", class = "deconflate_singular")
    # First-order (log-linear) values, also the start of the snapshot solver.
    b0 <- first_order_target(v, imp$measure, overall_risk)
    beta0 <- as.vector(solve(A, b0))
    if (method == "first_order") {
      adj <- exp(beta0)
      resid <- max(abs(as.vector(A %*% beta0) - b0))
      cond <- kappa(A, exact = TRUE)
    } else {
      st <- if (!is.null(start) && length(start) == n && all(is.finite(start))) start else beta0
      sol <- solve_event(v, imp$measure, joint$cells, joint$prob, as_imp, ids, overall_risk, start = st)
      adj <- exp(sol$beta)
      resid <- sol$residual
      cond <- kappa(sol$jacobian, exact = TRUE)
    }
  }
  adj <- unname(adj)
  raw_dir <- sign(event_scale(v, imp$measure))
  sign_change <- is.finite(adj) & adj > 0 & raw_dir != 0 & abs(adj - 1) > 1e-12 &
    sign(log(pmax(adj, 1e-300))) != raw_dir
  if (warn && any(sign_change)) {
    cm_warn(sprintf("Adjusted hazard ratios are on the other side of 1 from the raw estimates for %s. The raw estimates are weaker than the associated diseases alone would produce; check the estimands and source populations.",
                    paste(ids[sign_change], collapse = ", ")), class = "deconflate_sign_change")
  }
  # The published approximation can divide by a quantity near or below zero.
  bad <- !is.finite(adj) | adj <= 0
  if (warn && any(bad)) {
    cm_warn(sprintf("Adjusted hazard ratios are not positive and finite for %s.",
                    paste(ids[bad], collapse = ", ")), class = "deconflate_nonfinite_warning")
  }
  attributable <- NULL
  if (!any(bad)) {
    if (method != "snapshot") joint <- event_joint(pop, joint, joint_start, ...)
    attributable <- event_attributable(log(adj), joint, overall_risk, allocate = allocate)
  }
  structure(list(
    adjusted = data.frame(disease = ids, measure = imp$measure, raw = v, adjusted = adj,
                          change = ifelse(imp$measure %in% ratio_m, adj / v - 1, NA_real_),
                          estimand = imp$estimand, adjusted_for = imp$adjusted_for,
                          stringsAsFactors = FALSE),
    attributable = attributable,
    diagnostics = data.frame(method = method, max_reconstruction_residual = resid,
                             condition_number = cond, n_sign_changes = sum(sign_change),
                             sign_changes = paste(ids[sign_change], collapse = ", "),
                             feasibility = feas, stringsAsFactors = FALSE),
    joint = joint, model = model, method = method, overall_risk = overall_risk,
    label = attr(imp, "label"), units = attr(imp, "units")
  ), class = "cm_event_result")
}

# The population part of a model (for joint fits and their checks).
as_population <- function(model) {
  pop <- model
  pop$impacts <- NULL
  pop$interactions <- NULL
  pop$distributions <- NULL
  class(pop) <- "cm_population"
  pop
}

# A converged joint distribution for the event model: the one supplied
# (checked), or a new fit (optionally started from an earlier fit).
event_joint <- function(pop, joint, joint_start = NULL, ...) {
  if (is.null(joint)) {
    dots <- list(...)
    fj <- dots[intersect(names(dots), joint_arg_names)]
    joint <- withCallingHandlers(do.call(fit_joint, c(list(pop, start = joint_start), fj)),
                                 deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
  } else {
    validate_joint(joint, pop)
  }
  if (!isTRUE(joint$converged)) {
    cm_abort(if (identical(joint$backend, "sampled")) {
      "The sampled joint distribution is unresolved; see its diagnostics, and try more samples or chains."
    } else {
      "The joint distribution did not converge; the pairwise associations may be jointly infeasible (see check_feasibility(method = 'lp'))."
    }, class = "deconflate_nonconvergence")
  }
  joint
}

# Log-scale first-order targets of each row: log ratios, and for a risk
# difference the log of the implied risk ratio at the overall risk.
first_order_target <- function(value, measure, overall_risk) {
  out <- log(pmax(value, 1e-300))
  rd <- measure == "RD"
  out[rd] <- log(pmax(1 + value[rd] / overall_risk, 0.05))
  out
}

# Baseline cumulative hazard h0 such that the population risk over the period
# equals `overall_risk`, for hazard multipliers `mult` of the combinations.
baseline_hazard <- function(mult, prob, overall_risk) {
  keep <- prob > 0
  mult <- mult[keep]
  prob <- prob[keep]
  hi <- -log(1 - overall_risk) / min(mult)
  stats::uniroot(function(h) sum(prob * (1 - exp(-h * mult))) - overall_risk,
                 c(0, hi), extendInt = "upX", tol = 1e-14)$root
}

# Solve the snapshot equations by damped Newton iterations with a numerical
# Jacobian. Equation i: the model value of row i's measure (within strata of
# its adjustment set, combined with weights d1 d0 / (d1 + d0) on the
# measure's own scale: log for ratios, the difference for RD) equals its
# value. Hazard and rate ratios do not depend on h0; risk ratios, odds
# ratios and risk differences do, through the overall risk.
# An event estimate on its own scale: log for ratios, the difference for RD.
event_scale <- function(value, measure) {
  out <- value
  r <- measure != "RD"
  out[r] <- log(value[r])
  out
}

solve_event <- function(value, measure, cells, prob, estimands, ids, overall_risk, start,
                        tol = 1e-11, max_iter = 200L) {
  n <- length(value)
  target <- event_scale(value, measure)
  risk_based <- any(measure %in% c("RR", "OR", "RD"))
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
      cm_abort(sprintf("The estimate of %s is not identifiable: no stratum of its adjustment set contains animals with and without it.",
                       ids[i]), class = "deconflate_singular")
    }
    list(d1 = d1, d0 = d0, ok = ok, w = ifelse(ok, d1 * d0 / (d1 + d0), 0))
  })
  eval_F <- function(beta) {
    mult <- exp(as.vector(cells %*% beta))
    if (!all(is.finite(mult))) return(rep(Inf, n))
    R <- NULL
    if (risk_based) {
      h0 <- tryCatch(baseline_hazard(mult, prob, overall_risk), error = function(e) NA_real_)
      if (!is.finite(h0)) return(rep(Inf, n))
      R <- 1 - exp(-h0 * mult)
    }
    vapply(seq_len(n), function(i) {
      x <- cells[, i]
      bi <- base[[i]]
      val <- if (measure[i] %in% c("HR", "rate_ratio")) mult else R
      m1 <- rowsum(prob * x * val, strata[[i]], reorder = FALSE)[, 1][bi$ok] / bi$d1[bi$ok]
      m0 <- rowsum(prob * (1 - x) * val, strata[[i]], reorder = FALSE)[, 1][bi$ok] / bi$d0[bi$ok]
      e <- switch(measure[i],
        HR = , rate_ratio = , RR = log(m1) - log(m0),
        OR = log(m1) - log(1 - m1) - log(m0) + log(1 - m0),
        RD = m1 - m0)
      sum(bi$w[bi$ok] * e) / sum(bi$w[bi$ok]) - target[i]
    }, numeric(1))
  }
  jac <- function(beta) {
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
  if (!all(is.finite(f))) {
    beta <- rep(0, n)
    f <- eval_F(beta)
  }
  for (it in seq_len(max_iter)) {
    if (max(abs(f)) < tol) break
    Jm <- jac(beta)
    step <- tryCatch(solve(Jm, f), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) break
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
    cm_abort(sprintf("Could not solve the snapshot model (max residual %.2e): the estimates may be incompatible with the overall risk or with each other.", resid),
             class = "deconflate_nonconvergence")
  }
  list(beta = unname(beta), residual = resid, jacobian = jac(beta))
}

# The event risk attributable to disease, and its Shapley allocation over
# disease combinations. `beta`: log adjusted hazard ratios.
event_attributable <- function(beta, joint, overall_risk, allocate = TRUE, max_present = Inf) {
  ids <- joint$diseases
  rel <- exp(as.vector(joint$cells %*% beta))
  h0 <- baseline_hazard(rel, joint$prob, overall_risk)
  r0 <- 1 - exp(-h0)
  summ <- data.frame(overall_risk = overall_risk, disease_free_risk = r0,
                     attributable = overall_risk - r0,
                     attributable_fraction = (overall_risk - r0) / overall_risk,
                     unallocated = NA_real_, stringsAsFactors = FALSE)
  by <- NULL
  if (allocate) {
    b <- stats::setNames(beta, ids)
    loss <- function(x) 1 - exp(-h0 * exp(sum(b * x[ids]))) - r0
    sh <- withCallingHandlers(shapley_by_cell(joint, loss, max_present = max_present),
                              deconflate_incomplete_allocation = function(w) invokeRestart("muffleWarning"))
    by <- data.frame(disease = ids, attributable = sh$shapley, share = sh$share,
                     stringsAsFactors = FALSE)
    summ$unallocated <- summ$attributable - attr(sh, "allocated")
  }
  list(summary = summ, by_disease = by, baseline_hazard = h0)
}

#' @export
print.cm_event_result <- function(x, ...) {
  cat(sprintf("<cm_event_result> %s; method: %s", x$label %||% "event impacts", x$method))
  cat("\n\n")
  a <- x$adjusted
  tab <- data.frame(disease = a$disease, measure = a$measure, raw = a$raw,
                    adjusted_hr = a$adjusted, stringsAsFactors = FALSE)
  if (!is.null(a$lower)) {
    tab$lower <- a$lower
    tab$upper <- a$upper
  }
  if (any(a$estimand != "snapshot_crude")) tab$estimand <- a$estimand
  print(tab, row.names = FALSE, digits = 4)
  ar <- x$attributable
  if (!is.null(ar)) {
    s <- ar$summary
    cat(sprintf("\nOverall risk %.4g; disease-free risk %.4g; attributable to disease %.4g (%.1f%% of the overall risk)",
                s$overall_risk, s$disease_free_risk, s$attributable, 100 * s$attributable_fraction))
    if (!is.null(s$attributable_lower)) {
      cat(sprintf("; 95%% interval [%.4g, %.4g]", s$attributable_lower, s$attributable_upper))
    }
    cat("\n")
    if (!is.null(ar$by_disease)) {
      cat("\nAttributable risk by disease (Shapley allocation):\n")
      print(ar$by_disease, row.names = FALSE, digits = 4)
    }
  }
  d <- x$diagnostics
  cat(sprintf("\nDiagnostics: residual %.2e, condition number %.3g, sign changes %d%s\n",
              d$max_reconstruction_residual, d$condition_number, d$n_sign_changes,
              if (d$n_sign_changes) sprintf(" (%s)", d$sign_changes) else ""))
  print_unknown_pairs(x$unknown_pairs)
  print_draws_and_notes(x)
  invisible(x)
}
