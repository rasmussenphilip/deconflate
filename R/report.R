#' Contribution table
#'
#' One row per outcome and disease with raw and adjusted impacts, the
#' disease's burden (Shapley share of the expected loss, including any
#' interaction share) and, if `economics` is given, the attributed gap and
#' its monetary value.
#'
#' @param result A [deconflate()] result.
#' @param economics Optional list with `observed` and `unit_value` (see
#'   [value_losses()]). Hazard-ratio outcomes are valued with
#'   [attributable_risk()] (`observed` = overall risk, `unit_value` = value
#'   per animal removed).
#' @return A data frame. Burden columns are `NA` for hazard-ratio outcomes.
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
               by = c("outcome", "disease"), all.x = TRUE, sort = FALSE)
  names(out)[names(out) == "total"] <- "burden"
  if (!is.null(economics)) {
    ev <- evaluate_economics(result, economics, allocate = TRUE)
    if (!is.null(ev$by_disease)) {
      out <- merge(out, ev$by_disease, by = c("outcome", "disease"), all.x = TRUE, sort = FALSE)
    }
  }
  ord <- order(match(out$outcome, unique(a$outcome)), match(out$disease, result$model$diseases$id))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Gaps and values for every outcome named in economics$observed: productivity
# gaps for additive outcomes and attributable risk for hazard-ratio outcomes
# (observed = overall risk, unit_value = value per animal removed).
evaluate_economics <- function(result, economics, allocate = TRUE) {
  obs <- economics$observed
  if (is.null(obs) || is.null(names(obs))) cm_abort("`economics$observed` must be a named vector.")
  uv <- economics$unit_value
  a <- result$adjusted
  hr_out <- unique(a$outcome[a$scale == "hazard_ratio"])
  add_obs <- obs[!names(obs) %in% hr_out]
  hr_obs <- obs[names(obs) %in% hr_out]
  summ <- NULL
  by_disease <- NULL
  by_outcome <- NULL
  if (length(add_obs)) {
    pg <- productivity_gap(result, add_obs)
    vl <- value_losses(pg, uv, 0)
    summ <- pg$summary[, c("outcome", "observed", "disease_free", "gap")]
    by_outcome <- vl$by_outcome[, c("outcome", "value")]
    by_disease <- vl$by_disease
  }
  for (o in names(hr_obs)) {
    if (is.null(uv) || !(o %in% names(uv))) cm_abort(sprintf("No unit value for outcome '%s'.", o))
    ar <- attributable_risk(result, overall_risk = hr_obs[[o]], outcome = o,
                            unit_value = uv[[o]], allocate = allocate)
    summ <- rbind(summ, data.frame(outcome = o, observed = ar$summary$overall_risk,
                                   disease_free = ar$summary$disease_free_risk,
                                   gap = ar$summary$attributable, stringsAsFactors = FALSE))
    by_outcome <- rbind(by_outcome, data.frame(outcome = o, value = ar$summary$value,
                                               stringsAsFactors = FALSE))
    if (allocate) {
      by_disease <- rbind(by_disease, data.frame(outcome = o, disease = ar$by_disease$disease,
                                                 gap = ar$by_disease$attributable,
                                                 value = ar$by_disease$value,
                                                 stringsAsFactors = FALSE))
    }
  }
  list(summary = summ, by_outcome = by_outcome, by_disease = by_disease,
       total = sum(by_outcome$value) + sum(economics$additional %||% 0))
}

#' Summarise an adjustment result
#'
#' @param object A [deconflate()] result.
#' @param economics Optional economics list (see [value_losses()]). For
#'   outcomes on the hazard-ratio scale, `observed` is the overall risk of
#'   the event (e.g. the annual culling rate as a proportion) and
#'   `unit_value` the value per animal removed; see [attributable_risk()].
#' @param ... Unused.
#' @return A `summary.cm_result` list with `method`, `diagnostics`,
#'   `totals` (expected loss per outcome, raw and adjusted, plus gaps and
#'   values if `economics` is given; hazard-ratio outcomes appear only when
#'   valued), `total_value` and `contributions` ([contribution_table()]).
#' @export
summary.cm_result <- function(object, economics = NULL, ...) {
  a <- object$adjusted
  P <- stats::setNames(object$model$diseases$prob, object$model$diseases$id)
  b <- attribute_burden(object)
  add <- a[a$scale != "hazard_ratio", , drop = FALSE]
  totals <- if (nrow(add)) {
    do.call(rbind, lapply(split(add, add$outcome), function(x) {
      data.frame(outcome = x$outcome[1],
                 raw_loss = sum(x$raw * P[x$disease]),
                 adjusted_loss = sum(b$total[b$outcome == x$outcome[1]]),
                 stringsAsFactors = FALSE)
    }))
  } else {
    data.frame(outcome = character(0), raw_loss = numeric(0), adjusted_loss = numeric(0),
               stringsAsFactors = FALSE)
  }
  totals$reduction <- 1 - totals$adjusted_loss / totals$raw_loss
  total_value <- NULL
  if (!is.null(economics)) {
    ev <- evaluate_economics(object, economics, allocate = FALSE)
    totals <- merge(totals, ev$summary, by = "outcome", all = TRUE, sort = FALSE)
    totals <- merge(totals, ev$by_outcome, by = "outcome", all.x = TRUE, sort = FALSE)
    total_value <- ev$total
    attr(totals, "total_value") <- total_value
  }
  rownames(totals) <- NULL
  structure(list(method = object$method, diagnostics = object$diagnostics,
                 totals = totals, total_value = total_value,
                 contributions = contribution_table(object, economics)),
            class = "summary.cm_result")
}

#' @export
print.summary.cm_result <- function(x, ...) {
  cat(sprintf("Comorbidity adjustment (method: %s)\n\n", x$method))
  print(x$totals, row.names = FALSE, digits = 4)
  tv <- x$total_value %||% attr(x$totals, "total_value")
  if (!is.null(tv)) cat(sprintf("\nTotal value (including additional costs): %.2f\n", tv))
  d <- x$diagnostics
  if (any(d$n_sign_changes > 0)) {
    cat("\nSign changes (raw impacts smaller than associated diseases imply):\n")
    for (r in which(d$n_sign_changes > 0)) cat(sprintf("  %s: %s\n", d$outcome[r], d$sign_changes[r]))
  }
  invisible(x)
}
