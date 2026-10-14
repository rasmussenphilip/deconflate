#' Convert an odds ratio to a joint probability
#'
#' Solves the quadratic of Rasmussen et al. (2022), eqs. 7-11, for the 2x2
#' table with marginal probabilities `p1` and `p2` and odds ratio `or`, and
#' returns `P(d1 and d2)`. A numerically stable form of the quadratic formula
#' is used, and the root is selected by feasibility rather than by sign, so
#' that odds ratios near or below 1 are handled correctly.
#'
#' @param or Odds ratio (positive).
#' @param p1,p2 Marginal probabilities, strictly between 0 and 1.
#' @return The joint probability `P(d1 and d2)`.
#' @noRd
or_to_joint <- function(or, p1, p2) {
  if (!is.finite(or) || or <= 0) cm_abort("Odds ratios must be positive and finite.", class = "deconflate_infeasible")
  if (abs(or - 1) < 1e-10) return(p1 * p2)
  a <- (or - 1) * p2
  b <- -((or - 1) * (p1 + p2) + 1)
  cc <- or * p1
  disc <- b^2 - 4 * a * cc
  q <- -0.5 * (b + if (b >= 0) sqrt(disc) else -sqrt(disc))
  roots <- c(q / a, cc / q) * p2          # candidate values of P(d1 and d2)
  lo <- max(0, p1 + p2 - 1)
  hi <- min(p1, p2)
  ok <- is.finite(roots) & roots >= lo - 1e-12 & roots <= hi + 1e-12
  if (!any(ok)) {
    cm_abort("No feasible 2x2 table for this odds ratio and these marginals.",
             class = "deconflate_infeasible")
  }
  min(max(roots[ok][1], lo), hi)
}

#' Odds ratio implied by a joint probability
#'
#' @param p11 Joint probability `P(d1 and d2)`.
#' @param p1,p2 Marginal probabilities.
#' @return The odds ratio.
#' @noRd
joint_to_or <- function(p11, p1, p2) {
  (p11 * (1 - p1 - p2 + p11)) / ((p1 - p11) * (p2 - p11))
}

#' Joint probability of a disease pair from any supported measure
#'
#' @param measure `"OR"`, `"RR"`, `"RD"`, `"phi"` or `"prevalence"` (see
#'   [cm_associations()]; the directional measures treat disease 1 as the
#'   outcome and disease 2 as the comparison).
#' @param value Association value.
#' @param p1,p2 Marginal probabilities of disease 1 and disease 2.
#' @return `P(d1 and d2)`. Errors with class `deconflate_infeasible` if the
#'   measure is incompatible with the marginals.
#' @noRd
association_to_joint <- function(measure, value, p1, p2) {
  p11 <- switch(measure,
    OR = return(or_to_joint(value, p1, p2)),
    RR = value * p1 / (1 - p2 + value * p2) * p2,
    RD = (value * (1 - p2) + p1) * p2,
    prevalence = value * p2,
    phi = p1 * p2 + value * sqrt(p1 * (1 - p1) * p2 * (1 - p2)),
    cm_abort(sprintf("Cannot convert measure '%s' to a joint probability.", measure), class = "deconflate_invalid_input")
  )
  lo <- max(0, p1 + p2 - 1)
  hi <- min(p1, p2)
  if (!is.finite(p11) || p11 < lo - 1e-12 || p11 > hi + 1e-12) {
    cm_abort(sprintf("%s = %g is incompatible with marginal probabilities %g and %g.",
                     measure, value, p1, p2), class = "deconflate_infeasible")
  }
  min(max(p11, lo), hi)
}

#' Pairwise 2x2 tables for every disease pair
#'
#' The tables that the crude pair rows of the associations give (one disease
#' in `disease`, one in `associated_with`, no `among`, not adjusted for
#' diseases of the model). Several rows of one pair are pooled (see
#' [deconflate()]).
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param pool `"random"` (default) or `"fixed"`: how several rows of one
#'   pair are pooled.
#' @return A data frame with one row per unordered pair: the measure (or
#'   `"pooled"`), its `status` (`"specified"`: a crude pair row;
#'   `"other rows"`: only other kinds of rows mention the pair, e.g. within a
#'   subgroup; `"unknown"`: no row mentions it), the association rows used,
#'   the joint probability `p11` (`NA` unless specified), the implied odds
#'   ratio, and the excess probabilities
#'   `ep_2_given_1 = P(d2 | d1) - P(d2 | not d1)` and `ep_1_given_2`.
#' @noRd
pair_tables <- function(model, pool = c("random", "fixed")) {
  check_population(model)
  pair_tables_at(model, cm_match_arg(pool))
}

# pair_tables() with the pooling of the draws (frozen weights, between-study
# deviates).
pair_tables_at <- function(model, pool = "random", frozen = NULL, between = NULL) {
  ids <- model$diseases$disease
  p <- stats::setNames(model$diseases$prob, ids)
  n <- length(ids)
  if (n < 2L) {
    return(data.frame(disease1 = character(0), disease2 = character(0),
                      measure = character(0), status = character(0), rows = character(0),
                      p11 = numeric(0), or = numeric(0),
                      ep_2_given_1 = numeric(0), ep_1_given_2 = numeric(0)))
  }
  pairs <- t(utils::combn(ids, 2))
  d1 <- pairs[, 1]
  d2 <- pairs[, 2]
  np <- length(d1)
  key <- pair_key(d1, d2)
  measure <- rep(NA_character_, np)
  status <- rep("unknown", np)
  rows <- rep(NA_character_, np)
  p11 <- rep(NA_real_, np)
  specs <- association_specs(model)
  # Pairs mentioned by any row.
  for (s in specs) {
    mem <- c(s$g, s$h)
    if (length(mem) >= 2L) {
      cb <- utils::combn(mem, 2)
      status[key %in% pair_key(cb[1, ], cb[2, ])] <- "other rows"
    }
  }
  groups <- Filter(function(g) g$pair, association_groups(specs))
  info <- NULL
  if (length(groups)) {
    ctx <- list(P = p)
    info <- group_targets(groups, ctx, model, pool, frozen, between)
    check_exact_rows(groups, info)
    for (j in seq_along(groups)) {
      s <- groups[[j]]$specs[[1]]
      k <- match(pair_key(s$g, s$h), key)
      status[k] <- "specified"
      measure[k] <- if (length(groups[[j]]$specs) > 1L) "pooled" else s$measure
      rows[k] <- paste(groups[[j]]$rows, collapse = ", ")
      p11[k] <- info[[j]]$target
    }
  }
  p1 <- unname(p[d1])
  p2 <- unname(p[d2])
  out <- data.frame(
    disease1 = d1, disease2 = d2, measure = measure, status = status, rows = rows,
    p11 = p11,
    or = ifelse(is.na(p11), NA_real_, joint_to_or(p11, p1, p2)),
    ep_2_given_1 = p11 / p1 - (p2 - p11) / (1 - p1),
    ep_1_given_2 = p11 / p2 - (p1 - p11) / (1 - p2),
    stringsAsFactors = FALSE
  )
  attr(out, "pooling") <- list(groups = groups, info = info, pool = pool,
                               frozen = freeze_pooling(groups, info))
  out
}

# n x n matrix of joint probabilities P(j and k) from pairwise tables
# (diagonal = marginal probabilities).
pairwise_joint_matrix <- function(pt, ids, P) {
  J <- diag(P, nrow = length(ids))
  dimnames(J) <- list(ids, ids)
  J[cbind(pt$disease1, pt$disease2)] <- pt$p11
  J[cbind(pt$disease2, pt$disease1)] <- pt$p11
  J
}

# Covariance matrix of the disease indicators from the pairwise tables:
# Sigma[j, k] = P(j and k) - p_j p_k, Sigma[j, j] = p_j (1 - p_j).
pair_covariance <- function(J) {
  p <- diag(J)
  S <- J - outer(p, p)
  diag(S) <- p * (1 - p)
  S
}

# Coefficients of impact rows about one disease (no subgroup) on the
# diseases' own impacts, from the covariance matrix of the disease
# indicators: for a crude row of disease i, the excess probabilities
# P(k | i) - P(k | not i); for a row adjusted for diseases S, the coefficient
# of D_i in the population linear projection of D_k on (D_i, D_S), with 0 for
# the diseases in S.
pairwise_coefficients <- function(Sigma, specs, ids) {
  n <- length(ids)
  A <- matrix(0, length(specs), n, dimnames = list(NULL, ids))
  for (r in seq_along(specs)) {
    s <- specs[[r]]
    i <- match(s$g, ids)
    S <- match(if (s$adj_all) setdiff(ids, s$g) else s$adj, ids)
    X <- c(i, S)
    sxx <- Sigma[X, X, drop = FALSE]
    rc <- if (!all(is.finite(sxx))) NA_real_ else if (length(X) == 1L) 1 else rcond(sxx)
    if (!is.finite(rc) || rc < 1e-10 || sxx[1, 1] <= 0) {
      cm_abort(sprintf("%s: the estimate is not identifiable: its adjustment set is (nearly) collinear with the disease in this population.",
                       s$where), class = "deconflate_singular", table = "impacts", row = s$row)
    }
    A[r, i] <- 1
    others <- setdiff(seq_len(n), X)
    if (length(others)) A[r, others] <- solve(sxx, Sigma[X, others, drop = FALSE])[1, ]
  }
  A
}
