#' @export
print.cm_population <- function(x, ...) {
  cat("<cm_population>\n")
  print_population_lines(x)
  invisible(x)
}

print_population_lines <- function(x) {
  cat(sprintf("  Diseases: %d (%s)\n", nrow(x$diseases),
              paste(x$diseases$id, collapse = ", ")))
  pt <- tryCatch(pair_tables(x), error = function(e) NULL)
  if (!is.null(pt) && nrow(pt)) {
    k <- sum(pt$status == "unknown")
    cat(sprintf("  Disease pairs: %d (%d with an association, %d unknown)\n", nrow(pt),
                nrow(pt) - k, k))
  }
  if (!is.null(x$associations) && any(x$associations$adjusted %in% TRUE)) {
    cat(sprintf("  Covariate-adjusted association measures used as marginal: %d\n",
                sum(x$associations$adjusted %in% TRUE)))
  }
  if (!is.null(x$three_way) && nrow(x$three_way)) {
    cat(sprintf("  Three-way terms: %d\n", nrow(x$three_way)))
  }
}

#' @export
print.cm_model <- function(x, ...) {
  cat("<cm_model>\n")
  print_population_lines(x)
  im <- x$impacts
  if (!is.null(im)) {
    lab <- attr(im, "label")
    un <- attr(im, "units")
    kind <- impact_kind(im)
    cat(sprintf("  Impacts: %s%s (%s)\n", lab %||% "(unlabelled)",
                if (is.null(un) || is.na(un)) "" else sprintf(" [%s]", un),
                if (kind == "event") "event impacts: use event_model = TRUE" else "additive"))
    est <- table(im$estimand)
    cat(sprintf("  Estimands: %s\n",
                paste(sprintf("%s: %d", names(est), as.integer(est)), collapse = "; ")))
    if (kind == "event") {
      ms <- table(im$measure)
      cat(sprintf("  Measures: %s\n", paste(sprintf("%s: %d", names(ms), as.integer(ms)), collapse = "; ")))
    }
  }
  if (!is.null(x$interactions)) {
    cat(sprintf("  Interactions: %d\n", nrow(x$interactions)))
  }
  nd <- length(x$distributions)
  if (nd) cat(sprintf("  Uncertain inputs (with a distribution): %d\n", nd))
  invisible(x)
}

#' @export
print.cm_result <- function(x, ...) {
  cat(sprintf("<cm_result> method: %s", x$method))
  if (!is.null(x$label)) cat(sprintf("; %s", x$label))
  if (!is.null(x$units) && !is.na(x$units)) cat(sprintf(" [%s]", x$units))
  cat("\n\n")
  a <- x$adjusted
  cols <- c("disease", "raw", "adjusted")
  if (!is.null(a$lower)) cols <- c(cols, "lower", "upper")
  cols <- c(cols, "change")
  if (any(a$estimand != "crude")) cols <- c(cols, "estimand")
  print(a[, cols], row.names = FALSE, digits = 4)
  tt <- x$totals
  cat(sprintf("\nRaw sum: %.4g; adjusted total: %.4g", tt$raw_sum, tt$adjusted_total))
  if (!is.null(tt$adjusted_total_lower)) {
    cat(sprintf(" (95%% interval %.4g to %.4g)", tt$adjusted_total_lower, tt$adjusted_total_upper))
  }
  if (isTRUE(tt$interaction_total != 0)) cat(sprintf("; interactions: %.4g", tt$interaction_total))
  cat("\n")
  d <- x$diagnostics
  cat(sprintf("Diagnostics: residual %.2e, condition number %.3g, sign changes %d%s\n",
              d$max_reconstruction_residual, d$condition_number, d$n_sign_changes,
              if (d$n_sign_changes) sprintf(" (%s)", d$sign_changes) else ""))
  print_unknown_pairs(x$unknown_pairs)
  print_draws_and_notes(x)
  invisible(x)
}

print_unknown_pairs <- function(u) {
  if (is.null(u) || !nrow(u)) return(invisible(NULL))
  fo <- u$fitted_or[is.finite(u$fitted_or)]
  fit <- if (length(fo) == 1L) {
    sprintf("; the global fit gave it an odds ratio of %.3g", fo)
  } else if (length(fo)) {
    sprintf("; the global fit gave them odds ratios from %.3g to %.3g", min(fo), max(fo))
  } else ""
  cat(sprintf("Unknown pairs: %d without an association%s (see $unknown_pairs).\n", nrow(u), fit))
  invisible(NULL)
}

#' @export
print.cm_joint <- function(x, ...) {
  cat(sprintf("<cm_joint> %d diseases, %d combinations (%s backend)\n", length(x$diseases),
              length(x$prob), x$backend %||% "exact"))
  cat(sprintf("  Converged: %s after %d %s (max residual %.2e)\n",
              x$converged, x$iterations,
              if (identical(x$backend, "sampled")) "calibration iterations" else "sweeps",
              x$max_residual))
  cat(sprintf("  Constrained pairs: %d\n", nrow(x$constrained_pairs)))
  if (!is.null(x$targets$three_way) && nrow(x$targets$three_way)) {
    cat(sprintf("  Three-way terms: %d\n", nrow(x$targets$three_way)))
  }
  if (identical(x$backend, "sampled")) {
    d <- x$diagnostics$summary
    cat(sprintf("  Sample: %d draws from %d chains, %d distinct combinations\n",
                d$n_samples, d$n_chains, d$n_unique))
    cat(sprintf("  Constraint residual: %.2e in the sample%s; max R-hat %.3f; min ESS %.0f\n",
                d$residual_sample,
                if (is.na(d$residual_calibrated)) "" else sprintf(", %.2e after raking", d$residual_calibrated),
                d$max_rhat, d$min_ess))
  }
  invisible(x)
}
