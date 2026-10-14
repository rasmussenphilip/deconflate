#' @export
print.cm_population <- function(x, ...) {
  cat("<cm_population>\n")
  print_population_lines(x)
  invisible(x)
}

print_population_lines <- function(x) {
  cat(sprintf("  Diseases: %d (%s)\n", nrow(x$diseases),
              paste(x$diseases$disease, collapse = ", ")))
  a <- x$associations
  if (!is.null(a) && nrow(a)) {
    cat(sprintf("  Association rows: %d (%s)\n", nrow(a),
                paste(sprintf("%s: %d", names(table(a$measure)), as.integer(table(a$measure))), collapse = "; ")))
  }
  pt <- tryCatch(pair_tables(x), error = function(e) NULL)
  if (!is.null(pt) && nrow(pt)) {
    k <- sum(pt$status == "unknown")
    o <- sum(pt$status == "other rows")
    cat(sprintf("  Disease pairs: %d (%d with a crude pair row, %s%d unknown)\n", nrow(pt),
                sum(pt$status == "specified"), if (o) sprintf("%d with other rows only, ", o) else "", k))
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
    adj <- sum(!is.na(im$adjusted_for))
    sub <- sum(!is.na(im$among))
    combo <- unique(im$disease[grepl(";", im$disease, fixed = TRUE)])
    cat(sprintf("  Impact rows: %d (%d adjusted, %d within subgroups)%s\n", nrow(im), adj, sub,
                if (length(combo)) sprintf("; interactions: %s", paste(combo, collapse = ", ")) else ""))
    if (kind == "event") {
      ms <- table(im$measure)
      cat(sprintf("  Measures: %s\n", paste(sprintf("%s: %d", names(ms), as.integer(ms)), collapse = "; ")))
    }
  }
  nd <- length(x$distributions)
  if (nd) cat(sprintf("  Uncertain inputs (with a distribution): %d\n", nd))
  invisible(x)
}

#' @export
print.cm_result <- function(x, ...) {
  cat(sprintf("<cm_result> %s", x$method))
  if (!is.null(x$label)) cat(sprintf("; %s", x$label))
  if (!is.null(x$units) && !is.na(x$units)) cat(sprintf(" [%s]", x$units))
  cat("\n\n")
  a <- x$adjusted
  cols <- c("disease", "raw", "adjusted")
  if (!is.null(a$lower)) cols <- c(cols, "lower", "upper")
  cols <- c(cols, "change")
  if (any(a$rows != 1L)) cols <- c(cols, "rows")
  tab <- a[, intersect(cols, names(a)), drop = FALSE]
  print(tab, row.names = FALSE, digits = 4)
  tt <- x$totals
  cat(sprintf("\nRaw sum: %s; adjusted total: %.4g", if (is.finite(tt$raw_sum)) sprintf("%.4g", tt$raw_sum) else "NA (several rows per disease)",
              tt$adjusted_total))
  if (!is.null(tt$adjusted_total_lower)) {
    cat(sprintf(" (95%% interval %.4g to %.4g)", tt$adjusted_total_lower, tt$adjusted_total_upper))
  }
  if (isTRUE(tt$interaction_total != 0)) cat(sprintf("; interactions: %.4g", tt$interaction_total))
  cat("\n")
  print_diagnostics(x)
  print_unknown_pairs(x$unknown_pairs)
  if (!is.null(x$problems)) {
    cat("\n")
    print(x$problems)
  }
  print_draws_and_notes(x)
  invisible(x)
}

print_diagnostics <- function(x) {
  d <- x$diagnostics
  cat(sprintf("Diagnostics: residual %.2e, condition number %.3g, sign changes %d%s\n",
              d$max_reconstruction_residual, d$condition_number, d$n_sign_changes,
              if (d$n_sign_changes) sprintf(" (%s)", d$sign_changes) else ""))
  h <- x$heterogeneity
  if (!is.null(h) && nrow(h)) {
    cat(sprintf("Pooled or over-determined items: %d (see $heterogeneity and $estimates)\n", nrow(h)))
  }
}

print_unknown_pairs <- function(u) {
  if (is.null(u) || !nrow(u)) return(invisible(NULL))
  fo <- u$fitted_or[is.finite(u$fitted_or)]
  fit <- if (length(fo) == 1L) {
    sprintf("; the fit gave it an odds ratio of %.3g", fo)
  } else if (length(fo)) {
    sprintf("; the fit gave them odds ratios from %.3g to %.3g", min(fo), max(fo))
  } else ""
  cat(sprintf("Unknown pairs: %d without an association%s (see $unknown_pairs).\n", nrow(u), fit))
  invisible(NULL)
}

#' @export
print.cm_event_result <- function(x, ...) {
  cat(sprintf("<cm_event_result> %s; %s", x$label %||% "event impacts", x$method))
  cat("\n\n")
  a <- x$adjusted
  tab <- data.frame(disease = a$disease, measure = a$measure, raw = a$raw,
                    adjusted_hr = a$adjusted, stringsAsFactors = FALSE)
  if (!is.null(a$lower)) {
    tab$lower <- a$lower
    tab$upper <- a$upper
  }
  if (any(a$rows != 1L)) tab$rows <- a$rows
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
  cat("\n")
  print_diagnostics(x)
  print_unknown_pairs(x$unknown_pairs)
  if (!is.null(x$problems)) {
    cat("\n")
    print(x$problems)
  }
  print_draws_and_notes(x)
  invisible(x)
}

#' @export
print.cm_joint <- function(x, ...) {
  cat(sprintf("<cm_joint> %d diseases, %d combinations (%s backend)\n", length(x$diseases),
              length(x$prob), x$backend %||% "exact"))
  cat(sprintf("  Converged: %s after %d %s (max residual %.2e)\n",
              x$converged, x$iterations,
              if (identical(x$backend, "sampled")) "calibration iterations" else "sweeps",
              x$max_residual))
  cat(sprintf("  Association rows: %d\n", x$n_rows %||% 0L))
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
