#' Check whether the pairwise associations are jointly feasible
#'
#' Each pairwise 2x2 table can be valid on its own while no distribution of
#' disease combinations satisfies all of them together. This function checks
#' joint feasibility before fitting the global model, and points to the
#' pairs that conflict.
#'
#' Two checks are available:
#' * `"triples"`: for every triple of diseases whose three pairs are all
#'   constrained, checks whether some probability of all three together
#'   makes every cell of the 2x2x2 table non-negative. This is exact for
#'   three diseases and a necessary condition in general. It needs no
#'   additional packages.
#' * `"lp"`: solves a linear programme over all 2^n combination
#'   probabilities that minimises the total absolute deviation from the
#'   pairwise joint probabilities, keeping marginals exact. A minimum of zero
#'   means the constraints are jointly feasible; otherwise the deviation
#'   needed for each pair identifies the conflicting associations. This is
#'   exact, and requires the `lpSolve` package.
#'
#' @param model A [cm_model()].
#' @param method `"auto"` (both checks if `lpSolve` is installed and the
#'   number of diseases is at most `max_lp_diseases`, otherwise triples
#'   only), `"triples"` or `"lp"`.
#' @param max_lp_diseases Largest number of diseases for the LP check.
#' @param tol Numerical tolerance.
#' @return A `cm_feasibility` object: `feasible` (`TRUE`, `FALSE`, or `NA`
#'   if only the necessary triple check ran and passed), `triples` (data
#'   frame of violating triples and their gaps), `lp` (`NULL` or a list with
#'   the minimum total deviation and per-pair deviations) and `method`.
#' @export
#' @examples
#' check_feasibility(example_supplement())
#'
#' # Three strongly linked diseases that cannot all be associated this way:
#' bad <- cm_model(
#'   cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
#'   cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05))
#' )
#' check_feasibility(bad, method = "triples")
check_feasibility <- function(model, method = c("auto", "triples", "lp"),
                              max_lp_diseases = 14L, tol = 1e-9) {
  check_model(model)
  method <- match.arg(method)
  ids <- model$diseases$id
  P <- stats::setNames(model$diseases$prob, ids)
  pt <- pair_tables(model)
  J <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  for (r in seq_len(nrow(pt))) {
    J[pt$disease1[r], pt$disease2[r]] <- J[pt$disease2[r], pt$disease1[r]] <- pt$p11[r]
  }

  # Triple screen.
  tri <- list()
  if (length(ids) >= 3L) {
    combos <- utils::combn(length(ids), 3)
    for (j in seq_len(ncol(combos))) {
      i3 <- combos[, j]
      a <- i3[1]; b <- i3[2]; c3 <- i3[3]
      pab <- J[a, b]; pac <- J[a, c3]; pbc <- J[b, c3]
      if (anyNA(c(pab, pac, pbc))) next
      lower <- max(0, pab + pac - P[[a]], pab + pbc - P[[b]], pac + pbc - P[[c3]])
      upper <- min(pab, pac, pbc, 1 - P[[a]] - P[[b]] - P[[c3]] + pab + pac + pbc)
      if (lower > upper + tol) {
        tri[[length(tri) + 1L]] <- data.frame(
          disease1 = ids[a], disease2 = ids[b], disease3 = ids[c3],
          gap = lower - upper, stringsAsFactors = FALSE)
      }
    }
  }
  triples <- if (length(tri)) do.call(rbind, tri) else
    data.frame(disease1 = character(0), disease2 = character(0),
               disease3 = character(0), gap = numeric(0))

  run_lp <- method == "lp" ||
    (method == "auto" && length(ids) <= max_lp_diseases &&
       requireNamespace("lpSolve", quietly = TRUE))
  lp_res <- NULL
  if (run_lp) {
    if (!requireNamespace("lpSolve", quietly = TRUE)) {
      cm_abort("method = 'lp' requires the lpSolve package: install.packages('lpSolve').")
    }
    if (length(ids) > max_lp_diseases) {
      cm_abort(sprintf("The LP check is limited to %d diseases (increase `max_lp_diseases`).", max_lp_diseases))
    }
    lp_res <- feasibility_lp(model, pt, P, tol)
  }

  feasible <- if (nrow(triples)) FALSE else if (!is.null(lp_res)) lp_res$feasible else NA
  structure(list(feasible = feasible, triples = triples, lp = lp_res,
                 method = if (run_lp) "lp + triples" else "triples"),
            class = "cm_feasibility")
}

feasibility_lp <- function(model, pt, P, tol) {
  n <- length(P)
  cells <- disease_cells(n)
  cons <- pt[!is.na(pt$p11), , drop = FALSE]
  k <- nrow(cons)
  if (k == 0L) {
    return(list(feasible = TRUE, objective = 0, status = 0L, pairs = NULL))
  }
  nc <- nrow(cells)
  idx <- match(c(cons$disease1, cons$disease2), names(P))
  # Rows: total, marginals, pairs. Columns: cells, then (s+, s-) per pair.
  A <- matrix(0, 1 + n + k, nc + 2 * k)
  A[1, seq_len(nc)] <- 1
  for (i in seq_len(n)) A[1 + i, seq_len(nc)] <- cells[, i]
  for (r in seq_len(k)) {
    A[1 + n + r, seq_len(nc)] <- cells[, idx[r]] * cells[, idx[k + r]]
    A[1 + n + r, nc + 2 * r - 1] <- 1
    A[1 + n + r, nc + 2 * r] <- -1
  }
  rhs <- c(1, unname(P), cons$p11)
  obj <- c(rep(0, nc), rep(1, 2 * k))
  sol <- lpSolve::lp("min", obj, A, rep("=", nrow(A)), rhs)
  if (sol$status != 0) {
    return(list(feasible = FALSE, objective = NA_real_, status = sol$status,
                pairs = NULL))
  }
  s <- sol$solution[nc + seq_len(2 * k)]
  adj <- s[seq(2, 2 * k, 2)] - s[seq(1, 2 * k, 2)]   # achieved minus target
  pairs <- data.frame(disease1 = cons$disease1, disease2 = cons$disease2,
                      target_p11 = cons$p11, adjustment = adj,
                      stringsAsFactors = FALSE)
  pairs <- pairs[order(-abs(pairs$adjustment)), , drop = FALSE]
  rownames(pairs) <- NULL
  list(feasible = sol$objval <= tol, objective = sol$objval, status = sol$status,
       pairs = pairs)
}

#' @export
print.cm_feasibility <- function(x, ...) {
  verdict <- if (isTRUE(x$feasible)) "jointly feasible" else if (isFALSE(x$feasible))
    "NOT jointly feasible" else "no conflict found by the triple screen (necessary condition only)"
  cat(sprintf("<cm_feasibility> %s [%s]\n", verdict, x$method))
  if (nrow(x$triples)) {
    cat("\nConflicting triples:\n")
    print(utils::head(x$triples, 10), row.names = FALSE, digits = 3)
  }
  if (!is.null(x$lp) && isFALSE(x$lp$feasible) && !is.null(x$lp$pairs)) {
    cat(sprintf("\nMinimum total deviation from the pairwise tables: %.3g\n", x$lp$objective))
    cat("Pairs that must move (largest first):\n")
    p <- x$lp$pairs[abs(x$lp$pairs$adjustment) > 1e-9, , drop = FALSE]
    print(utils::head(p, 10), row.names = FALSE, digits = 3)
  }
  invisible(x)
}
