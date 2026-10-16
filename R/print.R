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

# Printed numbers: 3 significant figures, each value on its own (so a
# column can show 0.0137 above 0.00158). Very large or very small values in
# scientific notation; missing values blank (`na`: the text for them in
# running text, e.g. "NA").
fmt_num <- function(x, digits = 3L, na = "") {
  out <- rep(na, length(x))
  ok <- is.finite(x)
  out[is.infinite(x)] <- ifelse(x[is.infinite(x)] > 0, "Inf", "-Inf")
  v <- signif(x[ok], digits)
  big <- v != 0 & (abs(v) >= 1e5 | abs(v) < 1e-4)
  f <- ifelse(v == 0, "0", formatC(v, digits = digits, format = "fg", flag = "#"))
  f[big] <- formatC(v[big], digits = digits - 1L, format = "e")
  out[ok] <- sub("\\.$", "", trimws(f))
  out
}

# Print a table with the numbers formatted by fmt_num(): double columns to
# 3 significant figures, integers as they are, missing values blank. The
# character columns are right-aligned like the numbers, except those named
# in `left`; `hide` blanks their headers.
print_table <- function(df, left = NULL, digits = 3L, hide = NULL) {
  out <- df
  for (j in seq_along(out)) {
    v <- out[[j]]
    out[[j]] <- if (is.double(v)) fmt_num(v, digits) else ifelse(is.na(v), "", as.character(v))
  }
  for (nm in left) {
    w <- max(nchar(c(nm, out[[nm]])))
    out[[nm]] <- formatC(out[[nm]], width = -w)
    names(out)[names(out) == nm] <- formatC(if (nm %in% hide) "" else nm, width = -w)
  }
  print(out, row.names = FALSE, right = TRUE)
  invisible(df)
}

# The impacts table: one row per disease and interaction (raw, adjusted,
# interval, change), with each disease's contribution and share when
# `contributions` is TRUE.
impact_table <- function(x, contributions = TRUE, rows = FALSE) {
  a <- x$adjusted
  out <- data.frame(disease = a$disease, stringsAsFactors = FALSE)
  if (!is.null(a$measure)) out$measure <- a$measure
  out$raw <- a$raw
  out$adjusted <- a$adjusted
  if (!is.null(a$lower)) {
    out$lower <- a$lower
    out$upper <- a$upper
  }
  out$change <- a$change
  if (rows && any(a$rows != 1L)) out$rows <- a$rows
  ct <- x$contributions
  if (contributions && !is.null(ct)) {
    k <- match(a$disease, ct$disease)
    out$contribution <- ct$contribution[k]
    out$share <- ct$share[k]
  }
  out
}

# The totals in one line (print) for either model.
totals_text <- function(x) {
  iv <- function(lo, hi) if (!is.null(lo) && isTRUE(is.finite(lo))) sprintf(" (95%% interval %s to %s)", fmt_num(lo), fmt_num(hi)) else ""
  f <- function(v) fmt_num(v %||% NA_real_, na = "NA")
  if (inherits(x, "cm_event_result")) {
    s <- x$attributable$summary
    given <- identical(x$risk_given, "baseline")
    txt <- sprintf("Overall risk %s%s; disease-free risk %s%s; attributable to disease %s%s, %s of the overall risk",
                   f(s$overall_risk), if (given) " (computed)" else "",
                   f(s$disease_free_risk), if (given) " (given)" else "",
                   f(s$attributable), iv(s$attributable_lower, s$attributable_upper),
                   f(s$attributable_fraction))
    if (isTRUE(is.finite(s$raw_attributable))) {
      txt <- sprintf("%s; raw attributable risk (diseases unrelated) %s", txt, f(s$raw_attributable))
    }
    return(paste0(txt, "."))
  }
  tt <- x$totals
  txt <- sprintf("Total: raw (diseases unrelated) %s; adjusted %s%s", f(tt$raw_total), f(tt$adjusted_total),
                 iv(tt$adjusted_total_lower, tt$adjusted_total_upper))
  if (isTRUE(tt$interaction_total != 0)) txt <- sprintf("%s, of which interactions %s", txt, f(tt$interaction_total))
  red <- reduction(tt$adjusted_total, tt$raw_total %||% NA_real_)
  if (is.finite(red)) txt <- sprintf("%s; reduction %s", txt, f(red))
  paste0(txt, ".")
}

# One line when the diagnostics need attention (print() leaves the rest to
# summary()): sign changes, a residual of the exact rows above 1e-6 (relative
# to the largest raw value, when that is above 1) or a condition number
# above 1e6.
attention_line <- function(x) {
  d <- x$diagnostics
  why <- character(0)
  if (isTRUE(d$n_sign_changes > 0)) why <- c(why, sprintf("sign changes for %s", d$sign_changes))
  scale <- max(1, abs(x$adjusted$raw), na.rm = TRUE)
  if (isTRUE(d$max_reconstruction_residual > 1e-6 * scale)) {
    why <- c(why, sprintf("the rows are reproduced only to %s", fmt_num(d$max_reconstruction_residual)))
  }
  if (isTRUE(d$condition_number > 1e6)) {
    why <- c(why, sprintf("condition number %s (the results are very sensitive to the inputs)", fmt_num(d$condition_number)))
  }
  if (!length(why)) return(character(0))
  sprintf("Check: %s; see summary().", paste(why, collapse = "; "))
}

unknown_pairs_line <- function(u) {
  if (is.null(u) || !nrow(u)) return(character(0))
  fo <- u$fitted_or[is.finite(u$fitted_or)]
  fit <- if (length(fo) == 1L) {
    sprintf("; the fit gave it an odds ratio of %s", fmt_num(fo))
  } else if (length(fo)) {
    sprintf("; the fit gave them odds ratios from %s to %s", fmt_num(min(fo)), fmt_num(max(fo)))
  } else ""
  sprintf("Unknown pairs: %d without an association%s (see $unknown_pairs).", nrow(u), fit)
}

# Lines of text wrapped to the console width (continuation lines indented
# by `exdent`).
cat_wrapped <- function(lines, exdent = 2L) {
  for (l in lines) cat(strwrap(l, width = getOption("width", 80L), exdent = exdent), sep = "\n")
}

#' @export
print.cm_result <- function(x, ...) {
  cat(sprintf("<cm_result> %s", x$method))
  if (!is.null(x$label)) cat(sprintf("; %s", x$label))
  if (!is.null(x$units) && !is.na(x$units)) cat(sprintf(" [%s]", x$units))
  cat("\n\n")
  print_result_body(x)
  invisible(x)
}

#' @export
print.cm_event_result <- function(x, ...) {
  cat(sprintf("<cm_event_result> %s; %s", x$label %||% "event impacts", x$method))
  cat("\n\n")
  print_result_body(x)
  invisible(x)
}

print_result_body <- function(x) {
  print_table(impact_table(x))
  cat("\n")
  cat_wrapped(c(totals_text(x), attention_line(x), unknown_pairs_line(x$unknown_pairs), uncertainty_line(x$draws)))
  if (!is.null(x$problems)) {
    cat("\n")
    print(x$problems)
  }
  print_notes(x$notes)
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
