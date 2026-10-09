#' Contribution table
#'
#' One row per disease with the raw estimate, the adjusted value and the
#' disease's contribution to the total, with 95% intervals when the result
#' has draws.
#'
#' * Additive impacts: the contribution to the aggregate (its Shapley share,
#'   including any interaction share), in the units of the impacts.
#' * Event impacts: the adjusted hazard ratio and the disease's share of the
#'   risk attributable to disease (Shapley allocation), as a proportion of
#'   animals.
#'
#' @param result A [deconflate()] result.
#' @return A data frame.
#' @export
#' @examples
#' contribution_table(deconflate(example_supplement()))
contribution_table <- function(result) {
  if (inherits(result, "cm_event_result")) {
    a <- result$adjusted
    by <- result$attributable$by_disease
    out <- data.frame(disease = a$disease, measure = a$measure, raw = a$raw,
                      adjusted_hr = a$adjusted, stringsAsFactors = FALSE)
    if (!is.null(a$lower)) {
      out$adjusted_lower <- a$lower
      out$adjusted_upper <- a$upper
    }
    if (!is.null(by)) {
      out$attributable <- by$attributable
      out$share <- by$share
      if (!is.null(by$lower)) {
        out$attributable_lower <- by$lower
        out$attributable_upper <- by$upper
      }
    }
    return(out)
  }
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  a <- result$adjusted
  ct <- result$contributions
  out <- data.frame(disease = a$disease, raw = a$raw, adjusted = a$adjusted, change = a$change,
                    contribution = ct$total, main = ct$main, interaction = ct$interaction,
                    share = ct$share, stringsAsFactors = FALSE)
  if (!is.null(a$lower)) {
    out$adjusted_lower <- a$lower
    out$adjusted_upper <- a$upper
    out$contribution_lower <- ct$lower
    out$contribution_upper <- ct$upper
  }
  out
}

#' Summarise an adjustment result
#'
#' @param object A [deconflate()] result.
#' @param ... Unused.
#' @return A `summary.cm_result` list with `method`, `label`, `units`,
#'   `totals` (additive impacts: naive and adjusted aggregate and the
#'   reduction; event impacts: overall, disease-free and attributable risk),
#'   `diagnostics`, `contributions` (see [contribution_table()]) and `notes`.
#' @export
summary.cm_result <- function(object, ...) {
  tt <- object$totals
  totals <- data.frame(raw_sum = tt$raw_sum, adjusted_total = tt$adjusted_total,
                       interaction_total = tt$interaction_total,
                       reduction = if (abs(tt$raw_sum) > 0) 1 - tt$adjusted_total / tt$raw_sum else NA_real_)
  if (!is.null(tt$adjusted_total_lower)) {
    totals$adjusted_total_lower <- tt$adjusted_total_lower
    totals$adjusted_total_upper <- tt$adjusted_total_upper
  }
  structure(list(method = object$method, label = object$label, units = object$units,
                 totals = totals, diagnostics = object$diagnostics,
                 contributions = contribution_table(object),
                 notes = object$notes),
            class = "summary.cm_result")
}

#' @rdname summary.cm_result
#' @export
summary.cm_event_result <- function(object, ...) {
  s <- object$attributable$summary
  structure(list(method = object$method, label = object$label, units = "hazard ratio; risks as proportions",
                 totals = s, diagnostics = object$diagnostics,
                 contributions = contribution_table(object),
                 notes = object$notes),
            class = "summary.cm_result")
}

#' @export
print.summary.cm_result <- function(x, ...) {
  lab <- if (!is.null(x$label)) sprintf(" - %s", x$label) else ""
  un <- if (!is.null(x$units) && !is.na(x$units)) sprintf(" [%s]", x$units) else ""
  cat(sprintf("Comorbidity adjustment (method: %s)%s%s\n\n", x$method, lab, un))
  print(x$totals, row.names = FALSE, digits = 4)
  cat("\n")
  print(x$contributions, row.names = FALSE, digits = 4)
  d <- x$diagnostics
  if (d$n_sign_changes > 0) cat(sprintf("\nSign changes: %s\n", d$sign_changes))
  cat(sprintf("Feasibility: %s\n", d$feasibility))
  if (length(x$notes)) cat(paste0("Note: ", x$notes, "\n"), sep = "")
  invisible(x)
}
