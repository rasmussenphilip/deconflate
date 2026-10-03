#' @export
print.cm_model <- function(x, ...) {
  pt <- pair_tables(x)
  cat("<cm_model>\n")
  cat(sprintf("  Diseases: %d (%s)\n", nrow(x$diseases),
              paste(x$diseases$id, collapse = ", ")))
  if (nrow(pt)) {
    tab <- table(pt$status)
    cat(sprintf("  Disease pairs: %d [%s]\n", nrow(pt),
                paste(sprintf("%s: %d", names(tab), as.integer(tab)), collapse = "; ")))
  }
  if (!is.null(x$associations) && any(x$associations$adjusted)) {
    cat(sprintf("  Covariate-adjusted association measures used as marginal: %d\n",
                sum(x$associations$adjusted)))
  }
  if (!is.null(x$impacts)) {
    cat(sprintf("  Outcomes: %s\n", paste(unique(x$impacts$outcome), collapse = ", ")))
  }
  if (!is.null(x$interactions)) {
    cat(sprintf("  Interactions: %d\n", nrow(x$interactions)))
  }
  invisible(x)
}

#' @export
print.cm_result <- function(x, ...) {
  cat(sprintf("<cm_result> method: %s\n\n", x$method))
  print(x$adjusted[, c("outcome", "disease", "raw", "adjusted", "change")],
        row.names = FALSE, digits = 4)
  cat("\nDiagnostics:\n")
  print(x$diagnostics[, c("outcome", "max_reconstruction_residual",
                          "n_sign_changes", "condition_number")],
        row.names = FALSE, digits = 3)
  invisible(x)
}

#' @export
print.cm_joint <- function(x, ...) {
  cat(sprintf("<cm_joint> %d diseases, %d combinations\n", length(x$diseases),
              length(x$prob)))
  cat(sprintf("  Converged: %s after %d sweeps (max residual %.2e)\n",
              x$converged, x$iterations, x$max_residual))
  cat(sprintf("  Constrained pairs: %d\n", nrow(x$constrained_pairs)))
  invisible(x)
}

#' @export
print.cm_mc <- function(x, ...) {
  cat(sprintf("<cm_mc> method: %s\n", x$method))
  cat(sprintf("  Draws: %d, rejected as infeasible: %d (%.1f%%)\n", x$n_draws,
              x$n_rejected, 100 * x$n_rejected / x$n_draws))
  if (!is.null(x$ess)) cat(sprintf("  Effective sample size: %.1f\n", x$ess))
  invisible(x)
}
