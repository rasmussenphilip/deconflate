#' Contribution table
#'
#' One row per disease with raw and adjusted impacts and the disease's
#' contribution to the aggregate (its Shapley share, including any
#' interaction share), in the units of the impacts.
#'
#' @param result A [deconflate()] result.
#' @return A data frame.
#' @export
#' @examples
#' contribution_table(deconflate(example_supplement()))
contribution_table <- function(result) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  a <- result$adjusted
  ct <- result$contributions
  data.frame(disease = a$disease, raw = a$raw, adjusted = a$adjusted, change = a$change,
             contribution = ct$total, main = ct$main, interaction = ct$interaction,
             share = ct$share, stringsAsFactors = FALSE)
}

#' Summarise an adjustment result
#'
#' @param object A [deconflate()] result.
#' @param ... Unused.
#' @return A `summary.cm_result` list with `method`, `label`, `units`,
#'   `totals` (naive and adjusted aggregate and the reduction),
#'   `diagnostics` and `contributions`.
#' @export
summary.cm_result <- function(object, ...) {
  tt <- object$totals
  totals <- data.frame(raw_sum = tt$raw_sum, adjusted_total = tt$adjusted_total,
                       interaction_total = tt$interaction_total,
                       reduction = if (abs(tt$raw_sum) > 0) 1 - tt$adjusted_total / tt$raw_sum else NA_real_)
  structure(list(method = object$method, label = object$label, units = object$units,
                 totals = totals, diagnostics = object$diagnostics,
                 contributions = contribution_table(object)),
            class = "summary.cm_result")
}

#' @export
print.summary.cm_result <- function(x, ...) {
  lab <- if (!is.null(x$label)) sprintf(" - %s", x$label) else ""
  un <- if (!is.null(x$units)) sprintf(" [%s]", x$units) else ""
  cat(sprintf("Comorbidity adjustment (method: %s)%s%s\n\n", x$method, lab, un))
  print(x$totals, row.names = FALSE, digits = 4)
  cat("\n")
  print(x$contributions, row.names = FALSE, digits = 4)
  d <- x$diagnostics
  if (d$n_sign_changes > 0) cat(sprintf("\nSign changes: %s\n", d$sign_changes))
  cat(sprintf("Feasibility: %s\n", d$feasibility))
  invisible(x)
}
