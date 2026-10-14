# Additive impacts: the outcome model and its rows ------------------------------
#
# The outcome of an animal with disease combination c is
#   y(c) = sum_d b_d c_d + sum_T delta_T all_T(c),
# with an impact b_d for each disease and an interaction delta_T for each
# combination T that appears as the `disease` of an impact row. Every impact
# row is a linear function of these unknowns (see ?cm_impacts), so the rows
# give a linear system A theta = y, solved exactly (one row per unknown) or by
# weighted least squares (several rows; rows without uncertainty exact).

# Terms of the outcome model: each disease, then each combination in impact
# rows (in order of appearance), with their members.
model_terms <- function(specs, ids) {
  combos <- unique(unlist(lapply(specs, function(s) if (length(s$g) > 1L) s$term)))
  members <- c(as.list(ids), lapply(combos, split_ids))
  names(members) <- c(ids, combos)
  members
}

# 0/1 indicator of every term on the cells (all members present).
term_matrix <- function(cells, members) {
  out <- matrix(0, nrow(cells), length(members), dimnames = list(NULL, names(members)))
  for (j in seq_along(members)) {
    m <- members[[j]]
    out[, j] <- as.numeric(rowSums(cells[, m, drop = FALSE]) == length(m))
  }
  out
}

all_present <- function(cells, v) {
  if (length(v)) rowSums(cells[, v, drop = FALSE]) == length(v) else rep(TRUE, nrow(cells))
}
none_present <- function(cells, v) {
  if (length(v)) rowSums(cells[, v, drop = FALSE]) == 0 else rep(TRUE, nrow(cells))
}

# Proper subsets of a set with at least two members (the smaller
# combinations of a combination).
smaller_combinations <- function(g) {
  if (length(g) < 3L) return(list())
  out <- list()
  for (k in 2:(length(g) - 1L)) out <- c(out, utils::combn(g, k, simplify = FALSE))
  out
}

# Coefficients of one impact row on the terms, on a joint distribution.
row_coefficients <- function(s, cells, prob, Fm, members) {
  A <- all_present(cells, s$apos) & none_present(cells, s$aneg) & prob > 0
  XG <- all_present(cells, s$g)
  if (!s$adjusted) {
    w1 <- prob * (A & XG)
    w0 <- prob * (A & none_present(cells, s$g))
    if (sum(w1) <= 0 || sum(w0) <= 0) {
      cm_abort(sprintf("%s: the comparison is not identifiable (no animals with %s, or none without, %s).",
                       s$where, s$term, if (length(c(s$apos, s$aneg))) "in the subgroup" else "in the population"),
               class = "deconflate_singular", table = "impacts", row = s$row)
    }
    return(colSums(Fm * w1) / sum(w1) - colSums(Fm * w0) / sum(w0))
  }
  idx <- which(A)
  w <- prob[idx]
  cols <- list(XG = as.numeric(XG[idx]))
  if (length(s$g) > 1L) {
    for (d in s$g) cols[[d]] <- cells[idx, d]
    for (sub in smaller_combinations(s$g)) cols[[term_label(sub)]] <- as.numeric(all_present(cells[idx, , drop = FALSE], sub))
  }
  for (d in s$adj) cols[[d]] <- cells[idx, d]
  if (s$adj_all) {
    for (t in setdiff(names(members), c(s$term, names(cols)))) cols[[t]] <- Fm[idx, t]
  }
  Z <- cbind(1, do.call(cbind, cols))
  sw <- sqrt(w)
  q <- qr(Z * sw, tol = 1e-9)
  keep <- sort(q$pivot[seq_len(q$rank)])
  if (!(2L %in% keep)) {
    cm_abort(sprintf("%s: the estimate is not identifiable: %s is collinear with its adjustment set %s.",
                     s$where, s$term, if (length(c(s$apos, s$aneg))) "in the subgroup" else "in this population"),
             class = "deconflate_singular", table = "impacts", row = s$row)
  }
  B <- qr.coef(qr(Z[, keep, drop = FALSE] * sw, tol = 1e-12), Fm[idx, , drop = FALSE] * sw)
  B <- matrix(B, nrow = length(keep))
  out <- B[which(keep == 2L), ]
  out[!is.finite(out)] <- 0
  stats::setNames(out, names(members))
}

joint_coefficients <- function(specs, cells, prob, members) {
  Fm <- term_matrix(cells, members)
  A <- matrix(0, length(specs), length(members), dimnames = list(NULL, names(members)))
  for (r in seq_along(specs)) A[r, ] <- row_coefficients(specs[[r]], cells, prob, Fm, members)
  A
}

# Standard error of each impact row on its own scale (0 = exact).
impact_row_se <- function(model, specs, scale = function(v) v) {
  vapply(specs, function(s) row_se(row_dist(model, "impacts", s$row), scale), numeric(1))
}

# solve(A, y) for a square system, with a singular system as a classed error.
solve_exact <- function(A, y) {
  tryCatch(solve(A, y), error = function(e) {
    cm_abort(sprintf("Could not solve the adjustment equations: %s", conditionMessage(e)), class = "deconflate_singular")
  })
}

# Solve A theta = y: exactly (rows determine the unknowns exactly), or by
# weighted least squares with weights 1 / se^2, rows with se = 0 being exact
# constraints. `labels`: row labels for messages; `terms`: unknowns.
solve_rows <- function(A, y, se, labels, terms, check = TRUE) {
  m <- nrow(A)
  k <- ncol(A)
  se[is.na(se)] <- 1
  exact <- se == 0
  if (check) identify_terms(A, terms)
  tol_y <- 1e-8 * max(1, abs(y))
  conflict <- function(rows, resid) {
    if (!check) return(invisible(NULL))
    o <- order(-abs(resid))
    cm_abort(sprintf("The impact rows %s cannot all hold (largest residuals %s); rows without uncertainty are exact. Give them distributions (they are then fitted by weighted least squares), or check them.",
                     paste(labels[rows][o][seq_len(min(4, length(o)))], collapse = ", "),
                     paste(signif(resid[o][seq_len(min(4, length(o)))], 3), collapse = ", ")),
             class = "deconflate_infeasible", table = "impacts",
             rows = as.integer(sub("^impacts:", "", rownames(A)[rows][o])), residuals = resid[o])
  }
  if (m == k) {
    theta <- solve_exact(A, y)
  } else if (m < k) {
    theta <- qr.coef(qr(A), y)
    theta[!is.finite(theta)] <- 0
  } else if (all(exact)) {
    theta <- qr.coef(qr(A), y)
    res <- y - as.vector(A %*% theta)
    if (max(abs(res)) > tol_y) conflict(seq_len(m), res)
  } else if (!any(exact)) {
    sw <- 1 / se
    theta <- qr.coef(qr(A * sw), y * sw)
  } else {
    AE <- A[exact, , drop = FALSE]
    sv <- svd(AE, nu = nrow(AE), nv = k)
    d <- sv$d
    r <- sum(d > 1e-9 * max(d))
    Vr <- sv$v[, seq_len(r), drop = FALSE]
    theta0 <- as.vector(Vr %*% ((t(sv$u[, seq_len(r), drop = FALSE]) %*% y[exact]) / d[seq_len(r)]))
    resE <- y[exact] - as.vector(AE %*% theta0)
    if (max(abs(resE)) > tol_y) conflict(which(exact), resE)
    if (r < k) {
      N <- sv$v[, (r + 1):k, drop = FALSE]
      sw <- 1 / se[!exact]
      AU <- A[!exact, , drop = FALSE]
      phi <- qr.coef(qr((AU %*% N) * sw), (y[!exact] - as.vector(AU %*% theta0)) * sw)
      theta <- theta0 + as.vector(N %*% phi)
    } else {
      theta <- theta0
    }
  }
  theta <- unname(as.vector(theta))
  fitted <- as.vector(A %*% theta)
  resid <- y - fitted
  if (any(exact) && max(abs(resid[exact])) > tol_y) conflict(which(exact), resid[exact])
  std <- ifelse(exact, NA_real_, resid / se)
  df <- m - k
  Q <- if (df > 0 && any(!exact)) sum(std^2, na.rm = TRUE) else NA_real_
  list(theta = theta, fitted = fitted, resid = resid, std = std, Q = Q, df = df, exact = exact,
       condition = kappa(A, exact = TRUE), rank = qr(A, tol = 1e-9)$rank)
}

# Stop when the rows do not determine every unknown, naming those they do
# not determine.
identify_terms <- function(A, terms) {
  k <- ncol(A)
  r <- qr(A, tol = 1e-9)$rank
  if (r >= k) return(invisible(TRUE))
  sv <- svd(A, nv = k)
  N <- sv$v[, (r + 1):k, drop = FALSE]
  bad <- terms[apply(abs(N), 1, max) > 1e-8]
  cm_abort(sprintf("The impact rows do not determine %s separately (%d rows, %d unknowns, rank %d). Add a row that tells them apart (e.g. a row for a combination, or a row within a subgroup), or check for rows that repeat each other.",
                   join_and(bad), nrow(A), k, r), class = "deconflate_singular", terms = bad)
}

# Adjust additive impacts. `ctx`: either the pairwise tables (shortcut) or a
# joint distribution.
adjust_additive <- function(model, joint = NULL, pt = NULL, warn = TRUE, frozen_se = NULL) {
  imp <- model$impacts
  ids <- model$diseases$disease
  P <- stats::setNames(model$diseases$prob, ids)
  specs <- impact_specs(model)
  members <- model_terms(specs, ids)
  terms <- names(members)
  if (!is.null(joint)) {
    A <- joint_coefficients(specs, joint$cells, joint$prob, members)
    Fm <- term_matrix(joint$cells, members)
    pT <- colSums(Fm * joint$prob)
    J <- crossprod(joint$cells * joint$prob, joint$cells)
    dimnames(J) <- list(ids, ids)
    method <- "full model"
  } else {
    J <- pairwise_joint_matrix(pt, ids, P)
    A <- pairwise_coefficients(pair_covariance(J), specs, ids)
    pT <- P
    method <- "pairwise shortcut"
  }
  rownames(A) <- paste0("impacts:", vapply(specs, `[[`, numeric(1), "row"))
  y <- vapply(specs, `[[`, numeric(1), "value")
  se <- frozen_se %||% impact_row_se(model, specs)
  labels <- sprintf("row %s (%s)", vapply(specs, `[[`, numeric(1), "row"), vapply(specs, `[[`, character(1), "term"))
  sol <- solve_rows(A, y, se, labels, terms)
  theta <- stats::setNames(sol$theta, terms)
  n <- length(ids)
  b <- theta[ids]

  # Raw value of a term: its row when it has exactly one.
  term_rows <- lapply(terms, function(t) which(vapply(specs, `[[`, character(1), "term") == t))
  raw <- vapply(term_rows, function(r) if (length(r) == 1L) y[r] else NA_real_, numeric(1))
  combo <- !(terms %in% ids)
  sign_change <- !combo & is.finite(theta) & !is.na(raw) & abs(raw) > 1e-12 & abs(theta) > 1e-12 &
    sign(theta) != sign(raw)
  if (warn && any(sign_change)) {
    cm_warn(sprintf(
      "Adjusted impacts change sign for %s. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimates, and whether they come from populations with different comorbidity patterns.",
      paste(terms[sign_change], collapse = ", ")), class = "deconflate_sign_change", terms = terms[sign_change])
  }
  if (warn && any(!is.finite(theta))) {
    cm_warn(sprintf("Non-finite adjusted impacts for %s.", paste(terms[!is.finite(theta)], collapse = ", ")),
            class = "deconflate_nonfinite_warning", terms = terms[!is.finite(theta)])
  }

  # Contributions (closed-form Shapley values): each disease's own term and
  # an equal share of each interaction it is part of.
  main <- unname(P * b)
  inter <- numeric(n)
  for (j in which(combo)) {
    mem <- members[[j]]
    inter[match(mem, ids)] <- inter[match(mem, ids)] + theta[j] * pT[j] / length(mem)
  }
  contrib <- main + inter
  total <- sum(contrib)
  share <- if (is.finite(total) && abs(total) > 0) contrib / total else rep(NA_real_, n)
  raw_d <- raw[match(ids, terms)]

  estimates <- data.frame(table = "impacts", row = vapply(specs, `[[`, numeric(1), "row"),
                          item = vapply(specs, `[[`, character(1), "term"),
                          among = imp$among, adjusted_for = imp$adjusted_for, study = imp$study,
                          value = y, sd = ifelse(sol$exact, 0, se), fitted = sol$fitted,
                          std_residual = sol$std, stringsAsFactors = FALSE)
  structure(list(
    adjusted = data.frame(disease = terms, kind = ifelse(combo, "interaction", "disease"),
                          raw = raw, adjusted = unname(theta),
                          change = ifelse(!combo & !is.na(raw) & raw != 0, theta / raw - 1, NA_real_),
                          rows = lengths(term_rows), stringsAsFactors = FALSE),
    totals = data.frame(raw_sum = sum(P * raw_d), adjusted_total = total,
                        interaction_total = sum(inter)),
    contributions = data.frame(disease = ids, main = main, interaction = inter,
                               total = contrib, share = share, stringsAsFactors = FALSE),
    estimates = estimates,
    fit = list(Q = sol$Q, df = sol$df, p = if (is.finite(sol$Q)) stats::pchisq(sol$Q, sol$df, lower.tail = FALSE) else NA_real_),
    diagnostics = data.frame(method = method, max_reconstruction_residual = max(abs(sol$resid[sol$exact]), 0),
                             rank = sol$rank, condition_number = sol$condition,
                             n_sign_changes = sum(sign_change),
                             sign_changes = paste(terms[sign_change], collapse = ", "),
                             feasibility = NA_character_, stringsAsFactors = FALSE),
    conflation = list(A = A, terms = members), joint_pairs = J, model = model, joint = joint,
    method = method, label = attr(imp, "label"), units = attr(imp, "units"),
    row_se = se
  ), class = "cm_result")
}

# The published approximation of Rasmussen et al. (2022, eq. 16), kept for
# compare_methods() and the reproduce_*() functions: crude rows, one per
# disease, and a crude pair row for every pair.
adjust_published <- function(model, warn = TRUE, pool = "random", frozen = NULL) {
  imp <- model$impacts
  ids <- model$diseases$disease
  specs <- impact_specs(model)
  if (length(specs) != length(ids) || any(vapply(specs, function(s) length(s$g) > 1L || s$adjusted ||
                                                 length(c(s$apos, s$aneg)) > 0L, logical(1)))) {
    cm_abort("The published approximation is defined for one crude impact per disease (no combinations, subgroups or adjusted rows).",
             class = "deconflate_unsupported")
  }
  # (In draws: the central pooling weights and the draw's between-study
  # deviates, as for the model.)
  pt <- pair_tables_at(model, pool, frozen, model$between)
  if (anyNA(pt$p11)) {
    cm_abort(sprintf("The published approximation needs a crude pair row for every pair (missing: %s).",
                     paste(utils::head(paste(pt$disease1, pt$disease2, sep = "-")[is.na(pt$p11)], 5), collapse = ", ")),
             class = c("deconflate_unknown_pairs", "deconflate_unsupported"), pairs = missing_pairs(pt))
  }
  P <- stats::setNames(model$diseases$prob, ids)
  tri <- triple_screen(pt, P)
  if (nrow(tri)) {
    cm_abort("The pairwise associations are jointly infeasible; see check_feasibility().",
             class = "deconflate_infeasible", triples = tri)
  }
  J <- pairwise_joint_matrix(pt, ids, P)
  A <- pairwise_coefficients(pair_covariance(J), specs, ids)
  ord <- match(ids, vapply(specs, `[[`, character(1), "term"))
  A <- A[ord, , drop = FALSE]
  m_raw <- vapply(specs, `[[`, numeric(1), "value")[ord]
  n <- length(ids)
  conf <- as.vector((A - diag(n)) %*% m_raw)
  m_adj <- unname(ifelse(m_raw == 0, 0, m_raw^2 / (m_raw + conf)))
  sign_change <- is.finite(m_adj) & abs(m_raw) > 1e-12 & abs(m_adj) > 1e-12 & sign(m_adj) != sign(m_raw)
  main <- unname(P * m_adj)
  total <- sum(main)
  structure(list(
    adjusted = data.frame(disease = ids, kind = "disease", raw = m_raw, adjusted = m_adj,
                          change = ifelse(m_raw != 0, m_adj / m_raw - 1, NA_real_), rows = 1L,
                          stringsAsFactors = FALSE),
    totals = data.frame(raw_sum = sum(P * m_raw), adjusted_total = total, interaction_total = 0),
    contributions = data.frame(disease = ids, main = main, interaction = 0, total = main,
                               share = if (is.finite(total) && abs(total) > 0) main / total else NA_real_,
                               stringsAsFactors = FALSE),
    diagnostics = data.frame(method = "published", max_reconstruction_residual = {
      rec <- as.vector(A %*% m_adj); if (all(is.finite(rec))) max(abs(rec - m_raw)) else NA_real_
    }, rank = qr(A)$rank, condition_number = kappa(A, exact = TRUE),
    n_sign_changes = sum(sign_change), sign_changes = paste(ids[sign_change], collapse = ", "),
    feasibility = "triple screen passed (necessary condition only)", stringsAsFactors = FALSE),
    conflation = list(A = A), joint_pairs = J, model = model, joint = NULL, method = "published",
    label = attr(imp, "label"), units = attr(imp, "units")
  ), class = "cm_result")
}
