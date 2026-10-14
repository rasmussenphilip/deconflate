#' Contribution table
#'
#' One row per disease with the raw estimate (of its row, when it has
#' exactly one), the adjusted value and the disease's contribution to the
#' total, with 95% intervals when the result has draws. Interactions are in
#' the result's `adjusted` table; their contributions are shared among their
#' diseases here.
#'
#' * Additive impacts: the contribution to the aggregate (its Shapley share,
#'   including any interaction share), in the units of the impacts.
#' * Event impacts: the adjusted hazard ratio, the risk of an animal with
#'   only this disease and its risk ratio and risk difference against a
#'   disease-free animal, the relative change from the raw estimate (on the
#'   raw estimate's measure), and the disease's share of the risk
#'   attributable to disease (Shapley allocation), as a proportion of
#'   animals.
#'
#' In both cases the contribution is split into the disease's own part
#' (`main`: its Shapley value with the interactions set to zero) and its
#' share of the interactions (`interaction`: the rest; for event impacts it
#' also holds the small changes the interaction terms make, through
#' non-additive risks, to diseases outside them).
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
    ids <- if (!is.null(by)) by$disease else a$disease[a$kind == "disease"]
    k <- match(ids, a$disease)
    out <- data.frame(disease = ids, measure = a$measure[k], raw = a$raw[k],
                      adjusted_hr = a$adjusted[k], change = a$change[k], stringsAsFactors = FALSE)
    if (!is.null(a$lower)) {
      out$adjusted_lower <- a$lower[k]
      out$adjusted_upper <- a$upper[k]
    }
    if (!is.null(a$risk_alone)) {
      out$risk_alone <- a$risk_alone[k]
      out$rr_alone <- a$rr_alone[k]
      out$rd_alone <- a$rd_alone[k]
    }
    if (!is.null(by)) {
      out$attributable <- by$attributable
      out$main <- by$main
      out$interaction <- by$interaction
      out$share <- by$share
      if (!is.null(by$lower)) {
        out$attributable_lower <- by$lower
        out$attributable_upper <- by$upper
      }
    }
    return(out)
  }
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().", class = "deconflate_invalid_input")
  a <- result$adjusted
  ct <- result$contributions
  k <- match(ct$disease, a$disease)
  out <- data.frame(disease = ct$disease, raw = a$raw[k], adjusted = a$adjusted[k], change = a$change[k],
                    contribution = ct$total, main = ct$main, interaction = ct$interaction,
                    share = ct$share, stringsAsFactors = FALSE)
  if (!is.null(a$lower)) {
    out$adjusted_lower <- a$lower[k]
    out$adjusted_upper <- a$upper[k]
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
#'   `totals` (additive impacts: raw aggregate (the diseases unrelated),
#'   adjusted aggregate and the reduction; event impacts: overall,
#'   disease-free, raw and adjusted attributable risk and the reduction),
#'   `diagnostics`, `contributions` (see [contribution_table()]), `ranks`
#'   (with draws: how often each disease ranks 1st, 2nd, ... by
#'   contribution, and its mean rank; see `$ranks` of [deconflate()]) and
#'   `notes`.
#' @export
summary.cm_result <- function(object, ...) {
  tt <- object$totals
  totals <- data.frame(raw_total = tt$raw_total, adjusted_total = tt$adjusted_total,
                       interaction_total = tt$interaction_total,
                       reduction = if (isTRUE(abs(tt$raw_total) > 0)) 1 - tt$adjusted_total / tt$raw_total else NA_real_)
  if (!is.null(tt$adjusted_total_lower)) {
    totals$adjusted_total_lower <- tt$adjusted_total_lower
    totals$adjusted_total_upper <- tt$adjusted_total_upper
  }
  structure(list(method = object$method, label = object$label, units = object$units,
                 totals = totals, diagnostics = object$diagnostics,
                 contributions = contribution_table(object), ranks = object$ranks,
                 notes = object$notes),
            class = "summary.cm_result")
}

#' @rdname summary.cm_result
#' @export
summary.cm_event_result <- function(object, ...) {
  s <- object$attributable$summary
  if (!is.null(s$raw_attributable)) {
    s$reduction <- if (isTRUE(abs(s$raw_attributable) > 0)) 1 - s$attributable / s$raw_attributable else NA_real_
  }
  structure(list(method = object$method, label = object$label, units = "hazard ratio; risks as proportions",
                 totals = s, diagnostics = object$diagnostics,
                 contributions = contribution_table(object), ranks = object$ranks,
                 notes = c(if (identical(object$risk_given, "baseline"))
                   "The disease-free risk was given; the model computed the overall risk.", object$notes)),
            class = "summary.cm_result")
}

#' @export
print.summary.cm_result <- function(x, ...) {
  lab <- if (!is.null(x$label)) sprintf(" - %s", x$label) else ""
  un <- if (!is.null(x$units) && !is.na(x$units)) sprintf(" [%s]", x$units) else ""
  cat(sprintf("Comorbidity adjustment (%s)%s%s\n\n", x$method, lab, un))
  print(x$totals, row.names = FALSE, digits = 4)
  cat("\n")
  print(x$contributions, row.names = FALSE, digits = 4)
  if (!is.null(x$ranks)) {
    cat("\nRank by contribution over the draws (share of draws in each place):\n")
    r <- x$ranks[order(x$ranks$mean_rank), , drop = FALSE]
    p <- as.matrix(r[, grep("^rank_", names(r)), drop = FALSE])
    tab <- data.frame(disease = r$disease, `mean rank` = sprintf("%.2f", r$mean_rank),
                      matrix(sprintf("%.2f", p), nrow = nrow(p), dimnames = list(NULL, sub("^rank_", "", colnames(p)))),
                      check.names = FALSE, stringsAsFactors = FALSE)
    print(tab, row.names = FALSE, right = TRUE)
  }
  d <- x$diagnostics
  if (d$n_sign_changes > 0) cat(sprintf("\nSign changes: %s\n", d$sign_changes))
  cat(sprintf("Feasibility: %s\n", d$feasibility))
  if (length(x$notes)) cat(paste0("Note: ", x$notes, "\n"), sep = "")
  invisible(x)
}
