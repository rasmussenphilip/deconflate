#' Contribution table
#'
#' One row per disease with raw and adjusted impacts and the disease's
#' contribution to the aggregate (its Shapley share, including any
#' interaction share). With `valuation`, the attributed gap and its value are
#' added (see [productivity_gap()] and [value_losses()]).
#'
#' @param result A [deconflate()] result.
#' @param valuation Optional list with `observed`, `direction`, `effect` and
#'   (optionally) `unit_value`, passed to [productivity_gap()] and
#'   [value_losses()].
#' @return A data frame.
#' @export
#' @examples
#' contribution_table(deconflate(example_supplement()),
#'                    list(observed = 10000, direction = "decrease", effect = "percent",
#'                         unit_value = 0.3))
contribution_table <- function(result, valuation = NULL) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  a <- result$adjusted
  ct <- result$contributions
  out <- data.frame(disease = a$disease, raw = a$raw, adjusted = a$adjusted, change = a$change,
                    contribution = ct$total, main = ct$main, interaction = ct$interaction,
                    share = ct$share, stringsAsFactors = FALSE)
  if (!is.null(valuation)) {
    ev <- evaluate_valuation(result, valuation)
    out$gap <- ev$gap$attribution$gap
    if (!is.null(ev$value)) out$value <- ev$value$by_disease$value
  }
  out
}

# A valuation must say what the impacts mean: there are no defaults, because
# a wrong `effect` (e.g. percent read as a proportion) changes the gap by a
# factor of 100.
check_valuation <- function(valuation) {
  if (!is.list(valuation) || is.null(valuation$observed) || is.null(valuation$direction) ||
      is.null(valuation$effect)) {
    cm_abort("`valuation` must be a list with `observed`, `direction` and `effect` (and optionally `unit_value`).")
  }
  check_scalar_choice(valuation$direction, c("decrease", "increase"), "valuation$direction")
  check_scalar_choice(valuation$effect, c("proportion", "percent", "absolute"), "valuation$effect")
  invisible(valuation)
}

# Gap (and value) for one result from a valuation list.
evaluate_valuation <- function(result, valuation) {
  check_valuation(valuation)
  g <- productivity_gap(result, valuation$observed, direction = valuation$direction,
                        effect = valuation$effect)
  v <- if (!is.null(valuation$unit_value)) value_losses(g, valuation$unit_value) else NULL
  list(gap = g, value = v)
}

#' Summarise an adjustment result
#'
#' @param object A [deconflate()] result.
#' @param valuation Optional valuation list (see [contribution_table()]).
#' @param ... Unused.
#' @return A `summary.cm_result` list with `method`, `label`, `units`,
#'   `totals` (naive and adjusted aggregate, their difference, and the gap
#'   and value with `valuation`), `diagnostics` and `contributions`.
#' @export
summary.cm_result <- function(object, valuation = NULL, ...) {
  tt <- object$totals
  totals <- data.frame(raw_sum = tt$raw_sum, adjusted_total = tt$adjusted_total,
                       interaction_total = tt$interaction_total,
                       reduction = if (abs(tt$raw_sum) > 0) 1 - tt$adjusted_total / tt$raw_sum else NA_real_)
  if (!is.null(valuation)) {
    ev <- evaluate_valuation(object, valuation)
    totals$observed <- ev$gap$summary$observed
    totals$disease_free <- ev$gap$summary$disease_free
    totals$gap <- ev$gap$summary$gap
    if (!is.null(ev$value)) totals$value <- ev$value$value
  }
  structure(list(method = object$method, label = object$label, units = object$units,
                 totals = totals, diagnostics = object$diagnostics,
                 contributions = contribution_table(object, valuation)),
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
