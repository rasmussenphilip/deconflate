#' Fit the maximum-entropy distribution of disease combinations
#'
#' Starting from independence (the product of marginal probabilities; slide 12
#' of the 2023 seminar), iterative proportional fitting (IPF) repeatedly
#' rescales the 2^n combination probabilities so that every specified
#' pairwise 2x2 table is matched. From an independence start, IPF converges to
#' the maximum-entropy distribution subject to these constraints: a
#' log-linear model with all two-way terms and no explicit three-way or
#' higher log-linear interaction terms. Three-disease co-occurrence still
#' departs from independence through the pairwise terms.
#'
#' Three-way terms from [cm_three_way()] are put into the starting
#' distribution. IPF keeps every log-linear term of its starting distribution
#' that the constraints do not fix, so the fitted distribution matches all
#' pairwise tables and has the requested three-way terms.
#'
#' Pairs whose measure is `"unknown"` are not constrained. Pairs set to
#' `"independent"` are constrained to an odds ratio of 1.
#'
#' The constraint system is underdetermined for three or more diseases
#' (2^n - 1 free probabilities, n(n+1)/2 constraints), which is why the
#' maximum-entropy criterion (an assumption, not something the pairwise
#' evidence identifies) is needed to pick one solution. If the pairwise
#' tables are jointly infeasible, IPF does not converge; non-convergence
#' after `max_iter` sweeps suggests, but does not prove, infeasibility (see
#' [check_feasibility()]).
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param tol Convergence tolerance: maximum absolute deviation of any
#'   constrained cell probability, recomputed on the returned distribution.
#' @param max_iter Maximum number of full IPF sweeps.
#' @param max_diseases Safety limit on the number of diseases (the table has
#'   2^n cells).
#'
#' @return A `cm_joint` object: `cells` (2^n x n 0/1 matrix), `prob`
#'   (probability of each combination), `converged`, `iterations`,
#'   `max_residual` and the `targets` it was fitted to (used by
#'   [deconflate()] to check that a supplied joint matches the model).
#' @export
fit_joint <- function(model, tol = 1e-10, max_iter = 10000L, max_diseases = 20L) {
  check_population(model)
  ids <- model$diseases$id
  p <- model$diseases$prob
  n <- length(ids)
  if (n > max_diseases) {
    cm_abort(sprintf("%d diseases give 2^%d cells; increase `max_diseases` if intended.", n, n))
  }
  cells <- disease_cells(n)
  colnames(cells) <- ids
  logp <- as.vector(cells %*% log(p) + (1 - cells) %*% log(1 - p))
  tw <- model$three_way
  if (!is.null(tw) && nrow(tw)) {
    for (r in seq_len(nrow(tw))) {
      k <- match(c(tw$disease1[r], tw$disease2[r], tw$disease3[r]), ids)
      logp <- logp + log(tw$ratio[r]) * cells[, k[1]] * cells[, k[2]] * cells[, k[3]]
    }
  }
  prob <- exp(logp - max(logp))
  prob <- prob / sum(prob)

  pt <- pair_tables(model)
  cons <- pt[!is.na(pt$p11), , drop = FALSE]
  idx <- match(c(cons$disease1, cons$disease2), ids)
  i1 <- idx[seq_len(nrow(cons))]
  i2 <- idx[nrow(cons) + seq_len(nrow(cons))]
  # Group index 1..4 = (0,0), (0,1), (1,0), (1,1) for (d1, d2).
  groups <- lapply(seq_len(nrow(cons)), function(r) 2L * cells[, i1[r]] + cells[, i2[r]] + 1L)
  targets <- lapply(seq_len(nrow(cons)), function(r) {
    p1 <- p[i1[r]]
    p2 <- p[i2[r]]
    p11 <- cons$p11[r]
    c(1 - p1 - p2 + p11, p2 - p11, p1 - p11, p11)
  })
  margin_resid <- function(pr) max(abs(colSums(cells * pr) - p))
  pair_resid <- function(pr) {
    if (!length(groups)) return(0)
    max(vapply(seq_along(groups), function(r) {
      max(abs(vapply(1:4, function(j) sum(pr[groups[[r]] == j]), numeric(1)) - targets[[r]]))
    }, numeric(1)))
  }

  iter <- 0L
  for (it in seq_len(max_iter)) {
    iter <- it
    for (r in seq_along(groups)) {
      g <- groups[[r]]
      cur <- vapply(1:4, function(j) sum(prob[g == j]), numeric(1))
      fac <- ifelse(cur > 0, targets[[r]] / cur, 0)
      prob <- prob * fac[g]
    }
    # Marginals of every disease (needed for diseases without constrained
    # pairs, e.g. when all their pairs are unknown).
    for (i in seq_len(n)) {
      cur <- sum(prob[cells[, i] == 1])
      if (cur > 0 && cur < 1) {
        prob <- prob * ifelse(cells[, i] == 1, p[i] / cur, (1 - p[i]) / (1 - cur))
      }
    }
    if (max(pair_resid(prob), margin_resid(prob)) < tol) break
  }
  resid <- max(pair_resid(prob), margin_resid(prob))
  converged <- is.finite(resid) && resid < tol
  if (!converged) {
    cm_warn(sprintf(
      "IPF did not converge after %d sweeps (max residual %.2e). The pairwise tables may be jointly infeasible; see check_feasibility().",
      iter, resid), class = "deconflate_nonconvergence")
  }
  structure(list(cells = cells, prob = prob, diseases = ids,
                 converged = converged, iterations = iter, max_residual = resid,
                 constrained_pairs = cons[, c("disease1", "disease2", "measure", "status")],
                 targets = list(prob = p, pairs = cons[, c("disease1", "disease2", "p11")],
                                three_way = if (is.null(tw)) NULL else as.data.frame(tw))),
            class = "cm_joint")
}

# Check that a joint distribution belongs to a population: same diseases in
# the same order, matching marginals, pairwise tables and three-way terms,
# and converged. Impact-only changes do not invalidate a joint.
validate_joint <- function(joint, model, tol = 1e-7) {
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().")
  ids <- model$diseases$id
  stale <- function(why) {
    cm_abort(sprintf("The supplied joint distribution does not match the model (%s); refit it with fit_joint().", why))
  }
  if (!identical(joint$diseases, ids)) stale("different diseases or order")
  if (!isTRUE(joint$converged)) {
    cm_abort("The supplied joint distribution did not converge; its probabilities cannot be used.",
             class = "deconflate_nonconvergence")
  }
  if (max(abs(colSums(joint$cells * joint$prob) - model$diseases$prob)) > tol) {
    stale("disease probabilities")
  }
  pt <- pair_tables(model)
  cons <- pt[!is.na(pt$p11), , drop = FALSE]
  tg <- joint$targets$pairs
  if (is.null(tg)) stale("no recorded targets")
  key_m <- pair_key(cons$disease1, cons$disease2)
  key_j <- pair_key(tg$disease1, tg$disease2)
  if (!setequal(key_m, key_j)) stale("different constrained pairs")
  if (max(abs(cons$p11 - tg$p11[match(key_m, key_j)]), 0) > tol) stale("pairwise associations")
  tw_m <- model$three_way
  tw_j <- joint$targets$three_way
  sig <- function(tw) {
    if (is.null(tw) || !nrow(tw)) return(character(0))
    sort(paste(apply(tw[, c("disease1", "disease2", "disease3")], 1,
                     function(x) paste(sort(x), collapse = "|")),
               signif(tw$ratio, 10)))
  }
  if (!identical(sig(tw_m), sig(tw_j))) stale("three-way terms")
  invisible(joint)
}

#' Probabilities of disease combinations
#'
#' @param joint A [fit_joint()] result.
#' @param min_prob Drop combinations with probability below this value.
#' @return A data frame with one row per combination: a 0/1 column per
#'   disease, the number of diseases present and the probability, sorted by
#'   decreasing probability.
#' @export
combination_probs <- function(joint, min_prob = 0) {
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().")
  out <- as.data.frame(joint$cells)
  out$n_diseases <- rowSums(joint$cells)
  out$prob <- joint$prob
  out <- out[out$prob >= min_prob, , drop = FALSE]
  out <- out[order(-out$prob), , drop = FALSE]
  rownames(out) <- NULL
  out
}
