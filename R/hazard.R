# Event impacts: the snapshot hazard model -------------------------------------
#
# Event impacts (cm_impacts() with a `measure`) compare the risk of an event
# (e.g. death or culling) between groups of animals. They are not additive,
# so deconflate(..., event_model = TRUE) adjusts them with a hazard model:
# within the period, an animal with disease combination c has a constant
# hazard h0 * exp(lp(c)), lp(c) = sum_d beta_d c_d + sum_T gamma_T all_T(c),
# so its period risk is R(c) = 1 - exp(-h0 exp(lp(c))). The overall risk r
# fixes h0: sum_c p(c) R(c) = r (or the baseline risk r0: 1 - exp(-h0) =
# r0). The adjusted hazard ratios are exp(beta), each disease's own hazard
# multiplier, and exp(gamma), the interactions.

# Adjust event impacts with the snapshot model on a joint distribution. The
# unknowns are the log hazard multipliers of each disease (beta) and of each
# combination in the impact rows (gamma, interactions on the log-hazard
# scale). Each row is one equation (see ?cm_impacts); with several rows they
# are fitted by weighted least squares on each measure's scale (rows without
# uncertainty exact). `risk`: the overall risk, or (risk_kind = "baseline")
# the risk in disease-free animals. `start`: optional starting values (e.g.
# the central fit's). `details`: also the local row equations
# (`$conflation`).
adjust_event <- function(model, overall_risk, joint, warn = TRUE, allocate = TRUE, start = NULL,
                         frozen_se = NULL, risk_kind = "overall", details = TRUE) {
  risk <- overall_risk
  imp <- model$impacts
  if (is.null(imp)) cm_abort("The model has no impacts to adjust.", class = "deconflate_invalid_input")
  if (impact_kind(imp) != "event") {
    cm_abort("event_model = TRUE needs event impacts: an impact table with a measure column (HR, rate_ratio, RR, OR or RD).",
             class = "deconflate_unsupported")
  }
  ids <- model$diseases$disease
  specs <- impact_specs(model)
  members <- model_terms(specs, ids)
  terms <- names(members)
  k <- length(terms)
  cells <- joint$cells
  prob <- joint$prob
  Fm <- term_matrix(cells, members)
  setups <- lapply(specs, event_row_setup, cells = cells, prob = prob, members = members)
  measure <- vapply(specs, `[[`, character(1), "measure")
  y <- vapply(specs, `[[`, numeric(1), "value")
  target <- event_scale(y, measure)
  se <- frozen_se %||% vapply(specs, function(s) {
    row_se(row_dist(model, "impacts", s$row), function(v) event_scale(v, s$measure))
  }, numeric(1))
  se[is.na(se)] <- 1
  exact <- se == 0
  risk_based <- any(measure %in% c("RR", "OR", "RD"))
  row_values <- function(theta) {
    mult <- exp(as.vector(Fm %*% theta))
    if (!all(is.finite(mult))) return(rep(Inf, length(specs)))
    R <- NULL
    if (risk_based) {
      h0 <- tryCatch(baseline_hazard(mult, prob, risk, risk_kind), error = function(e) NA_real_)
      if (!is.finite(h0)) return(rep(Inf, length(specs)))
      R <- 1 - exp(-h0 * mult)
    }
    vapply(seq_along(specs), function(r) {
      event_row_value(setups[[r]], measure[r], theta, mult, R, prob)
    }, numeric(1))
  }
  eval_F <- function(theta) row_values(theta) - target
  labels <- sprintf("row %s (%s)", vapply(specs, `[[`, numeric(1), "row"), vapply(specs, `[[`, character(1), "term"))
  st <- if (!is.null(start) && length(start) == k && all(is.finite(start))) start else rep(0, k)
  sol <- solve_event_rows(eval_F, st, se, exact, labels, terms)
  theta <- stats::setNames(sol$theta, terms)
  adj <- exp(theta)
  term_of_row <- vapply(specs, `[[`, character(1), "term")
  raw <- vapply(terms, function(t) {
    r <- which(term_of_row == t)
    if (length(r) == 1L) y[r] else NA_real_
  }, numeric(1), USE.NAMES = FALSE)
  meas_t <- vapply(terms, function(t) {
    r <- which(term_of_row == t)
    if (length(r) == 1L) measure[r] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  combo <- !(terms %in% ids)
  attributable <- event_attributable(theta, joint, risk, members, allocate = allocate, risk_kind = risk_kind)
  # Each term on the risk scale: an animal with exactly the term's diseases
  # against a disease-free animal.
  al <- alone_risks(theta, members, attributable$baseline_hazard)
  # The adjusted value on the measure of the term's row (HR for hazard and
  # rate ratios; risk ratio, odds ratio or risk difference alone vs
  # disease-free otherwise), to compare with the raw value.
  comparable <- ifelse(meas_t %in% c("HR", "rate_ratio"), adj,
                       ifelse(meas_t == "RR", al$rr, ifelse(meas_t == "OR", al$or,
                                                           ifelse(meas_t == "RD", al$rd, NA_real_))))
  raw_dir <- ifelse(is.na(raw), 0, sign(event_scale(ifelse(is.na(raw), 1, raw), ifelse(is.na(meas_t), "HR", meas_t))))
  adj_sc <- ifelse(is.na(comparable), 0,
                   event_scale(ifelse(is.na(comparable), 1, comparable), ifelse(is.na(meas_t), "HR", meas_t)))
  adj_dir <- ifelse(abs(adj_sc) > 1e-12, sign(adj_sc), 0)
  sign_change <- !combo & is.finite(comparable) & raw_dir != 0 & adj_dir != 0 & adj_dir != raw_dir
  if (warn && any(sign_change)) {
    cm_warn(sprintf("Adjusted results are on the other side of 1 (of 0 for risk differences) from the raw estimates for %s. The raw estimates are weaker than the associated diseases alone would produce; check the estimates and source populations. deconflate(diagnose = TRUE) shows which inputs drive this.",
                    paste(terms[sign_change], collapse = ", ")), class = "deconflate_sign_change",
            terms = terms[sign_change])
  }
  conflation <- NULL
  if (details) {
    # Local row equations: how each row's value (on its measure's scale)
    # changes with each log-hazard term, at the adjusted values.
    h <- 1e-6
    A <- vapply(seq_len(k), function(j) {
      e <- rep(0, k)
      e[j] <- h
      (row_values(theta + e) - row_values(theta - e)) / (2 * h)
    }, numeric(length(specs)))
    A <- matrix(A, nrow = length(specs), dimnames = list(paste0("impacts:", vapply(specs, `[[`, numeric(1), "row")), terms))
    conflation <- list(A = A, terms = members, local = TRUE)
  }
  estimates <- data.frame(table = "impacts", row = vapply(specs, `[[`, numeric(1), "row"),
                          item = term_of_row, measure = measure,
                          among = imp$among, adjusted_for = imp$adjusted_for, study = imp$study,
                          value = y, sd = ifelse(exact, 0, se),
                          fitted = event_unscale(sol$f + target, measure),
                          std_residual = ifelse(exact, NA_real_, sol$f / se), stringsAsFactors = FALSE)
  structure(list(
    adjusted = data.frame(disease = terms, kind = ifelse(combo, "interaction", "disease"),
                          measure = meas_t, raw = raw, adjusted = unname(adj),
                          change = ifelse(!combo & is.finite(comparable) & !is.na(raw) & raw != 0,
                                          comparable / raw - 1, NA_real_),
                          rows = vapply(terms, function(t) sum(term_of_row == t), integer(1), USE.NAMES = FALSE),
                          risk_alone = al$risk, rr_alone = al$rr, rd_alone = al$rd, or_alone = al$or,
                          stringsAsFactors = FALSE),
    attributable = attributable,
    estimates = estimates,
    conflation = conflation,
    fit = list(Q = sol$Q, df = sol$df, p = if (is.finite(sol$Q)) stats::pchisq(sol$Q, sol$df, lower.tail = FALSE) else NA_real_),
    diagnostics = data.frame(method = "snapshot", max_reconstruction_residual = sol$residual,
                             condition_number = sol$condition, n_sign_changes = sum(sign_change),
                             sign_changes = paste(terms[sign_change], collapse = ", "),
                             feasibility = sprintf("joint distribution fitted (max residual %.1e)", joint$max_residual),
                             stringsAsFactors = FALSE),
    joint = joint, model = model, method = "snapshot",
    overall_risk = attributable$summary$overall_risk, risk_given = risk_kind,
    label = attr(imp, "label"), units = attr(imp, "units"), row_se = se, terms = members
  ), class = "cm_event_result")
}

# Period risk of an animal with exactly the diseases of each term, and its
# risk ratio, risk difference and odds ratio against a disease-free animal
# (baseline cumulative hazard h0).
alone_risks <- function(theta, members, h0) {
  r0 <- 1 - exp(-h0)
  lp <- vapply(members, function(m) {
    sum(theta[vapply(members, function(s) all(s %in% m), logical(1))])
  }, numeric(1))
  r <- 1 - exp(-h0 * exp(lp))
  list(risk = unname(r), rr = unname(r / r0), rd = unname(r - r0),
       or = unname((r / (1 - r)) / (r0 / (1 - r0))), r0 = r0)
}

# Precomputed groups of cells for one event row: the comparison groups
# (all vs none of G; with an adjustment set, every subset of G within each
# stratum of the adjustment diseases), their probabilities, the sign of each
# group in the contrast and the weight of each stratum. HR and rate-ratio
# rows adjusted for all are the term's own value: within a subgroup, the
# sum of the log multipliers of the terms made of G and the subgroup's
# present diseases (as the coefficient of an additive regression within the
# subgroup is the sum of the corresponding terms).
event_row_setup <- function(s, cells, prob, members) {
  if (s$adj_all && s$measure %in% c("HR", "rate_ratio")) {
    j <- which(vapply(members, function(m) all(s$g %in% m) && all(m %in% c(s$g, s$apos)), logical(1)))
    return(list(own = TRUE, j = unname(j)))
  }
  A <- all_present(cells, s$apos) & none_present(cells, s$aneg) & prob > 0
  k <- length(s$g)
  full <- 2^k - 1
  sub <- as.integer(cells[, s$g, drop = FALSE] %*% 2^(seq_len(k) - 1))
  if (s$adjusted) {
    S <- s$adj
    strat <- if (length(S)) as.integer(cells[, S, drop = FALSE] %*% 2^(seq_along(S) - 1)) else rep(0L, nrow(cells))
    need <- 0:full
  } else {
    strat <- rep(0L, nrow(cells))
    need <- c(0L, full)
  }
  use <- A & sub %in% need
  code <- strat * (full + 1) + sub
  d <- rowsum(prob[use], code[use])
  codes <- as.numeric(rownames(d))
  d <- d[, 1]
  st <- codes %/% (full + 1)
  su <- codes %% (full + 1)
  ok_st <- unique(st)[vapply(unique(st), function(z) all(need %in% su[st == z & d > 0]), logical(1))]
  if (!length(ok_st)) {
    cm_abort(sprintf("%s: the estimate is not identifiable: no %s contains animals with %s and without %s.",
                     s$where, if (s$adjusted) "stratum of its adjustment set" else "part of the population",
                     if (k > 1L) "every combination of its diseases" else "the disease", if (k > 1L) "them" else "it"),
             class = "deconflate_singular", table = "impacts", row = s$row)
  }
  keep <- st %in% ok_st
  codes <- codes[keep]
  d <- d[keep]
  st <- st[keep]
  su <- su[keep]
  bits <- vapply(su, function(x) sum(bitwAnd(x, 2^(seq_len(k) - 1)) > 0), numeric(1))
  idx <- which(use & code %in% codes)
  g <- match(code[idx], codes)
  w <- vapply(ok_st, function(z) 1 / sum(1 / d[st == z]), numeric(1))
  # Interaction contrast (inclusion-exclusion over the subsets of G); for a
  # crude row, all of G against none of G.
  sgn <- if (s$adjusted) (-1)^(k - bits) else ifelse(su == full, 1, -1)
  list(own = FALSE, idx = idx, g = g, d = d, sign = sgn, stratum = match(st, ok_st), w = w)
}

# Value of an event row on its measure's scale.
event_row_value <- function(su, measure, theta, mult, R, prob) {
  if (isTRUE(su$own)) return(sum(theta[su$j]))
  val <- if (measure %in% c("HR", "rate_ratio")) mult else R
  m <- rowsum(prob[su$idx] * val[su$idx], su$g)[, 1] / su$d
  e <- switch(measure,
    HR = , rate_ratio = , RR = log(m),
    OR = log(m) - log(1 - m),
    RD = m)
  contr <- rowsum(su$sign * e, su$stratum)[, 1]
  sum(su$w * contr) / sum(su$w)
}

# Solve the event equations F(theta) = 0 (one per row) by Gauss-Newton:
# Newton for one row per unknown; weighted least squares with rows without
# uncertainty as exact constraints otherwise. Numerical Jacobian.
solve_event_rows <- function(eval_F, start, se, exact, labels, terms, tol = 1e-11, max_iter = 200L) {
  k <- length(start)
  m <- length(se)
  jac <- function(theta) {
    h <- 1e-6
    J <- vapply(seq_len(k), function(j) {
      e <- rep(0, k)
      e[j] <- h
      (eval_F(theta + e) - eval_F(theta - e)) / (2 * h)
    }, numeric(m))
    matrix(J, nrow = m)
  }
  w <- ifelse(exact, 0, 1 / se^2)
  mu <- 1e8 * max(c(w, 1))
  merit <- function(f) if (all(is.finite(f))) sum(w * f^2) + mu * sum(f[exact]^2) else Inf
  theta <- start
  f <- eval_F(theta)
  if (!all(is.finite(f))) {
    theta <- rep(0, k)
    f <- eval_F(theta)
  }
  J0 <- jac(theta)
  if (!all(is.finite(J0))) cm_abort("The snapshot model could not be evaluated at the start.", class = "deconflate_nonconvergence")
  identify_terms(J0, terms)
  J <- J0
  for (it in seq_len(max_iter)) {
    if (m == k && max(abs(f)) < tol) break
    step <- tryCatch(solve_rows(J, -f, ifelse(exact, 0, se), labels, terms, check = FALSE)$theta,
                     error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) break
    lambda <- 1
    improved <- FALSE
    m0 <- merit(f)
    while (lambda >= 1e-8) {
      cand <- theta + lambda * step
      fc <- eval_F(cand)
      if (merit(fc) < m0) {
        improved <- TRUE
        break
      }
      lambda <- lambda / 2
    }
    if (!improved) break
    small <- max(abs(lambda * step)) < 1e-12 || (m0 - merit(fc)) <= 1e-15 * max(1, m0)
    theta <- cand
    f <- fc
    if (small && (m != k || max(abs(f)) < tol)) break
    J <- jac(theta)
  }
  resid <- if (any(exact)) max(abs(f[exact])) else 0
  if (!is.finite(resid) || resid > 1e-8 || !all(is.finite(f))) {
    cm_abort(sprintf("Could not solve the snapshot model (max residual of the exact rows %.2e): the estimates may be incompatible with the overall risk or with each other.", resid),
             class = "deconflate_nonconvergence", residual = resid)
  }
  df <- m - k
  std <- ifelse(exact, NA_real_, f / se)
  list(theta = unname(theta), f = f, residual = resid, df = df,
       Q = if (df > 0 && any(!exact)) sum(std^2, na.rm = TRUE) else NA_real_,
       condition = kappa(jac(theta), exact = TRUE))
}

# From an event measure's scale back to the measure.
event_unscale <- function(x, measure) ifelse(measure == "RD", x, exp(x))

# The population part of a model (for joint fits and their checks).
as_population <- function(model) {
  pop <- model
  pop$impacts <- NULL
  d <- pop$distributions
  pop$distributions <- if (length(d)) d[!startsWith(names(d), "impacts:")]
  if (!length(pop$distributions)) pop$distributions <- NULL
  class(pop) <- "cm_population"
  pop
}

# A converged joint distribution: the one supplied (checked), or a new fit
# (optionally started from an earlier fit).
event_joint <- function(pop, joint, joint_start = NULL, pool = "random", frozen = NULL, between = NULL, ...) {
  if (is.null(joint)) {
    dots <- list(...)
    fj <- dots[intersect(names(dots), joint_arg_names)]
    joint <- withCallingHandlers(do.call(fit_joint_at, c(list(pop, start = joint_start, pool = pool,
                                                               frozen = frozen, between = between), fj)),
                                 deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
  } else {
    validate_joint(joint, pop, pool)
  }
  if (!isTRUE(joint$converged)) {
    cm_abort(if (identical(joint$backend, "sampled")) {
      sprintf("The sampled joint distribution is unresolved (residual %.2e); see its diagnostics, and try more samples or chains. This does not show that the inputs are infeasible.",
              joint$max_residual)
    } else {
      sprintf("The joint distribution did not converge (max residual %.2e): the association rows may not be able to hold together; see check_feasibility().",
              joint$max_residual)
    }, class = "deconflate_nonconvergence", residual = joint$max_residual)
  }
  joint
}

# Baseline cumulative hazard h0 such that the population risk over the period
# equals `overall_risk`, for hazard multipliers `mult` of the combinations
# (risk_kind = "baseline": such that the disease-free risk equals it).
baseline_hazard <- function(mult, prob, overall_risk, risk_kind = "overall") {
  keep <- prob > 0
  mult <- mult[keep]
  prob <- prob[keep]
  if (!all(is.finite(mult) & mult > 0)) {
    cm_abort("The hazard multipliers are not all finite and positive, so the baseline hazard is undefined.",
             class = "deconflate_nonfinite")
  }
  # A baseline risk fixes the hazard of disease-free animals directly.
  if (identical(risk_kind, "baseline")) return(-log(1 - overall_risk))
  hi <- -log(1 - overall_risk) / min(mult)
  tryCatch(stats::uniroot(function(h) sum(prob * (1 - exp(-h * mult))) - overall_risk,
                          c(0, hi), extendInt = "upX", tol = 1e-14)$root,
           error = function(e) {
             cm_abort(sprintf("Could not find the baseline hazard: %s", conditionMessage(e)),
                      class = "deconflate_nonconvergence")
           })
}

# An event estimate on its own scale: log for ratios, the difference for RD.
event_scale <- function(value, measure) {
  out <- value
  r <- measure != "RD"
  out[r] <- log(value[r])
  out
}

# The event risk attributable to disease, and its Shapley allocation over
# disease combinations. `theta`: log hazard multipliers of the terms
# (`members`: the diseases of each term); `risk`: the overall risk, or the
# baseline (disease-free) risk with risk_kind = "baseline". The allocation is
# split into each disease's own part (the allocation with the interaction
# terms set to zero, same baseline hazard) and its share of the interactions
# (the rest).
event_attributable <- function(theta, joint, overall_risk, members, allocate = TRUE, max_present = Inf,
                               risk_kind = "overall") {
  ids <- joint$diseases
  Fm <- term_matrix(joint$cells, members)
  rel <- exp(as.vector(Fm %*% theta))
  h0 <- baseline_hazard(rel, joint$prob, overall_risk, risk_kind)
  r0 <- 1 - exp(-h0)
  overall <- if (identical(risk_kind, "baseline")) sum(joint$prob * (1 - exp(-h0 * rel))) else overall_risk
  summ <- data.frame(overall_risk = overall, disease_free_risk = r0,
                     attributable = overall - r0,
                     attributable_fraction = (overall - r0) / overall,
                     unallocated = NA_real_, stringsAsFactors = FALSE)
  by <- NULL
  if (allocate) {
    shapley_for <- function(th) {
      th <- unname(th)
      loss <- function(x) {
        lp <- 0
        for (j in seq_along(members)) if (all(x[members[[j]]] == 1)) lp <- lp + th[j]
        1 - exp(-h0 * exp(lp)) - r0
      }
      withCallingHandlers(shapley_by_cell(joint, loss, max_present = max_present),
                          deconflate_incomplete_allocation = function(w) invokeRestart("muffleWarning"))
    }
    sh <- shapley_for(theta)
    combo <- lengths(members) > 1L
    main <- if (any(combo)) shapley_for(ifelse(combo, 0, theta))$shapley else sh$shapley
    by <- data.frame(disease = ids, main = main, interaction = sh$shapley - main,
                     attributable = sh$shapley, share = sh$share, stringsAsFactors = FALSE)
    summ$unallocated <- summ$attributable - attr(sh, "allocated")
  }
  list(summary = summ, by_disease = by, baseline_hazard = h0)
}

