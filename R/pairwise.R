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
#' @export
#' @examples
#' # Supplementary File, Rasmussen et al. (2022): P(1 | 2) = 0.163
#' or_to_joint(2, 0.10, 0.15) / 0.15
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
#' @export
joint_to_or <- function(p11, p1, p2) {
  (p11 * (1 - p1 - p2 + p11)) / ((p1 - p11) * (p2 - p11))
}

#' Joint probability of a disease pair from any supported measure
#'
#' @param measure One of the measures in [cm_associations()] other than
#'   `"unknown"` (`"table"` rows store the table's odds ratio in `value`).
#' @param value Association value.
#' @param p1,p2 Marginal probabilities of `disease1` and `disease2`.
#' @return `P(d1 and d2)`. Errors with class `deconflate_infeasible` if the
#'   measure is incompatible with the marginals.
#' @export
association_to_joint <- function(measure, value, p1, p2) {
  p11 <- switch(measure,
    independent = p1 * p2,
    OR = ,
    table = return(or_to_joint(value, p1, p2)),
    RR = value * p1 / (1 - p2 + value * p2) * p2,
    RD = (value * (1 - p2) + p1) * p2,
    cond_prob = value * p2,
    phi = p1 * p2 + value * sqrt(p1 * (1 - p1) * p2 * (1 - p2)),
    cm_abort(sprintf("Cannot convert measure '%s' to a joint probability.", measure))
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
#' @param model A [cm_population()] or [cm_model()].
#' @return A data frame with one row per unordered pair: the measure used,
#'   whether it was specified or defaulted, the joint probability `p11`
#'   (`NA` for unknown pairs), the implied odds ratio, and the excess
#'   probabilities `ep_2_given_1 = P(d2 | d1) - P(d2 | not d1)` and
#'   `ep_1_given_2`.
#' @export
pair_tables <- function(model) {
  check_population(model)
  ids <- model$diseases$id
  p <- stats::setNames(model$diseases$prob, ids)
  n <- length(ids)
  if (n < 2L) {
    return(data.frame(disease1 = character(0), disease2 = character(0),
                      measure = character(0), status = character(0),
                      p11 = numeric(0), or = numeric(0),
                      ep_2_given_1 = numeric(0), ep_1_given_2 = numeric(0)))
  }
  pairs <- t(utils::combn(ids, 2))
  d1 <- pairs[, 1]
  d2 <- pairs[, 2]
  np <- length(d1)
  assoc <- model$associations
  akey <- if (is.null(assoc)) character(0) else pair_key(assoc$disease1, assoc$disease2)
  hit <- match(pair_key(d1, d2), akey)

  default <- if (model$missing_associations == "independent") "independent" else "unknown"
  measure <- rep(default, np)
  status <- rep(paste0(default, " (default)"), np)
  p11 <- if (default == "independent") p[d1] * p[d2] else rep(NA_real_, np)
  spec <- which(!is.na(hit))
  if (length(spec)) {
    a <- assoc[hit[spec], , drop = FALSE]
    measure[spec] <- a$measure
    status[spec] <- "specified"
    # Directional measures are defined with the row's own orientation.
    p11[spec] <- vapply(seq_along(spec), function(j) {
      if (a$measure[j] == "unknown") return(NA_real_)
      association_to_joint(a$measure[j], a$value[j], p[[a$disease1[j]]], p[[a$disease2[j]]])
    }, numeric(1))
  }
  p1 <- unname(p[d1])
  p2 <- unname(p[d2])
  p11 <- unname(p11)
  data.frame(
    disease1 = d1, disease2 = d2, measure = measure, status = status,
    p11 = p11,
    or = ifelse(is.na(p11), NA_real_, joint_to_or(p11, p1, p2)),
    ep_2_given_1 = p11 / p1 - (p2 - p11) / (1 - p1),
    ep_1_given_2 = p11 / p2 - (p1 - p11) / (1 - p2),
    stringsAsFactors = FALSE
  )
}

#' Excess-probability matrix from pairwise tables
#'
#' @param model A [cm_population()] or [cm_model()].
#' @return An n x n matrix `E` with `E[k, i] = P(k | i) - P(k | not i)`, the
#'   excess probability of disease `k` among animals with disease `i`
#'   (`ep_ki` in Rasmussen et al. 2022, eq. 14). The diagonal is 0. Errors if
#'   any pair is `"unknown"`: pairwise methods need every pair specified.
#' @export
excess_matrix <- function(model) {
  pt <- pair_tables(model)
  ids <- model$diseases$id
  if (anyNA(pt$p11)) {
    unk <- paste(pt$disease1, pt$disease2, sep = "-")[is.na(pt$p11)]
    cm_abort(sprintf(
      "Unknown associations for %s. Specify them, set missing_associations = 'independent', or use method = 'global'.",
      paste(unk, collapse = ", ")), class = c("deconflate_unknown_pairs", "deconflate_unsupported"))
  }
  E <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (r in seq_len(nrow(pt))) {
    E[pt$disease2[r], pt$disease1[r]] <- pt$ep_2_given_1[r]
    E[pt$disease1[r], pt$disease2[r]] <- pt$ep_1_given_2[r]
  }
  E
}

# n x n matrix of joint probabilities P(j and k) from pairwise tables
# (diagonal = marginal probabilities).
pairwise_joint_matrix <- function(model) {
  pt <- pair_tables(model)
  ids <- model$diseases$id
  J <- diag(model$diseases$prob, nrow = length(ids))
  dimnames(J) <- list(ids, ids)
  for (r in seq_len(nrow(pt))) {
    J[pt$disease1[r], pt$disease2[r]] <- pt$p11[r]
    J[pt$disease2[r], pt$disease1[r]] <- pt$p11[r]
  }
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

# Conflation matrix A (raw = A %*% adjusted) from the covariance matrix of the
# disease indicators, for impacts that are crude differences or coefficients
# of an additive regression adjusted for a set of diseases. Row i holds, for
# each omitted disease k, the coefficient of D_i in the population linear
# projection of D_k on (D_i, D_S). For crude estimates (S empty) this is the
# excess probability P(k | i) - P(k | not i).
projection_matrix <- function(Sigma, impacts, ids) {
  n <- length(ids)
  A <- diag(n)
  dimnames(A) <- list(ids, ids)
  for (i in seq_len(n)) {
    S <- if (impacts$estimand[i] == "adjusted_linear") {
      match(resolve_adjusted_for(impacts$adjusted_for[i], ids[i], ids), ids)
    } else integer(0)
    X <- c(i, S)
    sxx <- Sigma[X, X, drop = FALSE]
    rc <- if (!all(is.finite(sxx))) NA_real_ else if (length(X) == 1L) 1 else rcond(sxx)
    if (!is.finite(rc) || rc < 1e-10 || sxx[1, 1] <= 0) {
      cm_abort(sprintf("The estimand of %s is not identifiable: its adjustment set is (nearly) collinear with it in this population.",
                       ids[i]), class = "deconflate_singular")
    }
    others <- setdiff(seq_len(n), X)
    if (length(others)) {
      g <- solve(sxx, Sigma[X, others, drop = FALSE])
      A[i, others] <- g[1, ]
    }
    if (length(S)) A[i, S] <- 0
  }
  A
}
