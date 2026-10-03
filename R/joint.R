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
#' Pairs whose measure is `"unknown"` are not constrained. Pairs set to
#' `"independent"` are constrained to an odds ratio of 1.
#'
#' The constraint system is underdetermined for three or more diseases
#' (2^n - 1 free probabilities, n(n+1)/2 constraints), which is why the
#' maximum-entropy criterion is needed to pick one solution. If the pairwise
#' tables are jointly infeasible, IPF does not converge and the result is
#' flagged.
#'
#' @param model A [cm_model()].
#' @param tol Convergence tolerance: maximum absolute deviation of any
#'   constrained cell probability.
#' @param max_iter Maximum number of full IPF sweeps.
#' @param max_diseases Safety limit on the number of diseases (the table has
#'   2^n cells).
#'
#' @return A `cm_joint` object: `cells` (2^n x n 0/1 matrix), `prob`
#'   (probability of each combination), `converged`, `iterations` and
#'   `max_residual`.
#' @export
fit_joint <- function(model, tol = 1e-10, max_iter = 10000L, max_diseases = 20L) {
  check_model(model)
  ids <- model$diseases$id
  p <- model$diseases$prob
  n <- length(ids)
  if (n > max_diseases) {
    cm_abort(sprintf("%d diseases give 2^%d cells; increase `max_diseases` if intended.", n, n))
  }
  cells <- disease_cells(n)
  colnames(cells) <- ids
  prob <- as.vector(exp(cells %*% log(p) + (1 - cells) %*% log(1 - p)))

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

  converged <- nrow(cons) == 0L
  resid <- 0
  iter <- 0L
  while (!converged && iter < max_iter) {
    iter <- iter + 1L
    resid <- 0
    for (r in seq_along(groups)) {
      g <- groups[[r]]
      cur <- vapply(1:4, function(j) sum(prob[g == j]), numeric(1))
      resid <- max(resid, abs(cur - targets[[r]]))
      fac <- ifelse(cur > 0, targets[[r]] / cur, 0)
      prob <- prob * fac[g]
    }
    converged <- resid < tol
  }
  if (!converged) {
    cm_warn(sprintf(
      "IPF did not converge after %d sweeps (max residual %.2e). The pairwise tables may be jointly infeasible.",
      iter, resid), class = "deconflate_nonconvergence")
  }
  structure(list(cells = cells, prob = prob, diseases = ids,
                 converged = converged, iterations = iter, max_residual = resid,
                 constrained_pairs = cons[, c("disease1", "disease2", "measure", "status")]),
            class = "cm_joint")
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
