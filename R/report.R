#' Contribution table
#'
#' One row per outcome and disease with raw and adjusted impacts, the
#' disease's burden (Shapley share of the expected loss, including any
#' interaction share) and, if `economics` is given, the attributed gap and
#' its monetary value.
#'
#' @param result A [deconflate()] result.
#' @param economics Optional list with `observed` and `unit_value` (see
#'   [value_losses()]).
#' @return A data frame.
#' @export
#' @examples
#' contribution_table(deconflate(example_uk_dairy_2022(), method = "published"),
#'                    uk_dairy_2022_economics())
contribution_table <- function(result, economics = NULL) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  a <- result$adjusted
  b <- attribute_burden(result)
  out <- merge(a[, c("outcome", "disease", "raw", "adjusted", "change")],
               b[, c("outcome", "disease", "main", "interaction", "total", "share")],
               by = c("outcome", "disease"), sort = FALSE)
  names(out)[names(out) == "total"] <- "burden"
  if (!is.null(economics)) {
    pg <- productivity_gap(result, economics$observed)
    vl <- value_losses(pg, economics$unit_value, economics$additional %||% 0)
    out <- merge(out, vl$by_disease, by = c("outcome", "disease"), all.x = TRUE, sort = FALSE)
  }
  ord <- order(match(out$outcome, unique(a$outcome)), match(out$disease, result$model$diseases$id))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Summarise an adjustment result
#'
#' @param object A [deconflate()] result.
#' @param economics Optional economics list (see [value_losses()]).
#' @param ... Unused.
#' @return A `summary.cm_result` list with `method`, `diagnostics`,
#'   `totals` (expected loss per outcome, raw and adjusted, plus gaps and
#'   values if `economics` is given) and `contributions`
#'   ([contribution_table()]).
#' @export
summary.cm_result <- function(object, economics = NULL, ...) {
  a <- object$adjusted
  P <- stats::setNames(object$model$diseases$prob, object$model$diseases$id)
  b <- attribute_burden(object)
  totals <- do.call(rbind, lapply(split(a, a$outcome), function(x) {
    data.frame(outcome = x$outcome[1],
               raw_loss = sum(x$raw * P[x$disease]),
               adjusted_loss = sum(b$total[b$outcome == x$outcome[1]]),
               stringsAsFactors = FALSE)
  }))
  totals$reduction <- 1 - totals$adjusted_loss / totals$raw_loss
  if (!is.null(economics)) {
    pg <- productivity_gap(object, economics$observed)
    vl <- value_losses(pg, economics$unit_value, economics$additional %||% 0)
    totals <- merge(totals, pg$summary[, c("outcome", "observed", "disease_free", "gap")],
                    by = "outcome", all.x = TRUE, sort = FALSE)
    totals <- merge(totals, vl$by_outcome[, c("outcome", "value")], by = "outcome",
                    all.x = TRUE, sort = FALSE)
    attr(totals, "total_value") <- vl$total
  }
  rownames(totals) <- NULL
  structure(list(method = object$method, diagnostics = object$diagnostics,
                 totals = totals, contributions = contribution_table(object, economics)),
            class = "summary.cm_result")
}

#' @export
print.summary.cm_result <- function(x, ...) {
  cat(sprintf("Comorbidity adjustment (method: %s)\n\n", x$method))
  print(x$totals, row.names = FALSE, digits = 4)
  tv <- attr(x$totals, "total_value")
  if (!is.null(tv)) cat(sprintf("\nTotal value (including additional costs): %.2f\n", tv))
  d <- x$diagnostics
  if (any(d$n_sign_changes > 0)) {
    cat("\nSign changes (raw impacts smaller than associated diseases imply):\n")
    for (r in which(d$n_sign_changes > 0)) cat(sprintf("  %s: %s\n", d$outcome[r], d$sign_changes[r]))
  }
  invisible(x)
}
