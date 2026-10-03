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
    tab <- table(pt$status)
    cat(sprintf("  Disease pairs: %d [%s]\n", nrow(pt),
                paste(sprintf("%s: %d", names(tab), as.integer(tab)), collapse = "; ")))
  }
  if (!is.null(x$associations) && any(x$associations$adjusted %in% TRUE)) {
    cat(sprintf("  Covariate-adjusted association measures used as marginal: %d\n",
                sum(x$associations$adjusted %in% TRUE)))
  }
  if (!is.null(x$associations) && any(x$associations$corrected %in% TRUE)) {
    cat(sprintf("  Contingency tables with a zero-cell correction: %d\n",
                sum(x$associations$corrected %in% TRUE)))
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
    cat(sprintf("  Impacts: %s%s\n", lab %||% "(unlabelled)",
                if (is.null(un)) "" else sprintf(" [%s]", un)))
    est <- table(im$estimand)
    cat(sprintf("  Estimands: %s\n",
                paste(sprintf("%s: %d", names(est), as.integer(est)), collapse = "; ")))
  }
  if (!is.null(x$interactions)) {
    cat(sprintf("  Interactions: %d\n", nrow(x$interactions)))
  }
  invisible(x)
}

#' @export
print.cm_analyses <- function(x, ...) {
  cat(sprintf("<cm_analyses> %d analyses on one population\n", length(x$models)))
  print_population_lines(x$population)
  for (nm in names(x$models)) {
    im <- x$models[[nm]]$impacts
    un <- attr(im, "units")
    ni <- if (is.null(x$models[[nm]]$interactions)) 0L else nrow(x$models[[nm]]$interactions)
    cat(sprintf("  - %s%s%s\n", nm, if (is.null(un)) "" else sprintf(" [%s]", un),
                if (ni) sprintf(", %d interactions", ni) else ""))
  }
  invisible(x)
}

#' @export
print.cm_result <- function(x, ...) {
  cat(sprintf("<cm_result> method: %s", x$method))
  if (!is.null(x$label)) cat(sprintf("; %s", x$label))
  if (!is.null(x$units)) cat(sprintf(" [%s]", x$units))
  cat("\n\n")
  a <- x$adjusted
  cols <- c("disease", "raw", "adjusted", "change")
  if (any(a$estimand != "crude")) cols <- c(cols, "estimand")
  print(a[, cols], row.names = FALSE, digits = 4)
  tt <- x$totals
  cat(sprintf("\nRaw sum: %.4g; adjusted total: %.4g", tt$raw_sum, tt$adjusted_total))
  if (isTRUE(tt$interaction_total != 0)) cat(sprintf(" (interactions: %.4g)", tt$interaction_total))
  cat("\n")
  d <- x$diagnostics
  cat(sprintf("Diagnostics: residual %.2e, condition number %.3g, sign changes %d%s\n",
              d$max_reconstruction_residual, d$condition_number, d$n_sign_changes,
              if (d$n_sign_changes) sprintf(" (%s)", d$sign_changes) else ""))
  invisible(x)
}

#' @export
print.cm_joint <- function(x, ...) {
  cat(sprintf("<cm_joint> %d diseases, %d combinations\n", length(x$diseases),
              length(x$prob)))
  cat(sprintf("  Converged: %s after %d sweeps (max residual %.2e)\n",
              x$converged, x$iterations, x$max_residual))
  cat(sprintf("  Constrained pairs: %d\n", nrow(x$constrained_pairs)))
  if (!is.null(x$targets$three_way)) {
    cat(sprintf("  Three-way terms: %d\n", nrow(x$targets$three_way)))
  }
  invisible(x)
}

#' @export
print.cm_mc <- function(x, ...) {
  cat(sprintf("<cm_mc> method: %s\n", paste(x$method, collapse = ", ")))
  if (!is.null(x$label)) cat(sprintf("  Analysis: %s%s\n", x$label,
                                     if (is.null(x$units)) "" else sprintf(" [%s]", x$units)))
  cat(sprintf("  Draws: %d, rejected: %d (%.1f%%)\n", x$n_draws,
              x$n_rejected, 100 * x$n_rejected / x$n_draws))
  smp <- x$sampling %||% "random"
  if (identical(smp, "lhs") && !is.null(x$block)) {
    smp <- sprintf("Latin hypercube, %d replicate blocks", length(unique(x$block)))
  }
  if (!is.null(x$proposal)) smp <- paste0(smp, ", importance sampling of ",
                                         paste(names(x$proposal), collapse = ", "))
  cat(sprintf("  Sampling: %s\n", smp))
  if (!is.null(x$ess)) cat(sprintf("  Effective sample size: %.1f\n", x$ess))
  invisible(x)
}
