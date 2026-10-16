#' Contribution table
#'
#' One row per disease with the raw estimate (of its row, when it has
#' exactly one), the adjusted value and the disease's contribution to the
#' total, at full precision, for export. The columns are the same for every
#' result of a model (intervals are `NA` without draws), so a script that
#' writes the table out does not depend on the run. Interactions are in the
#' result's `adjusted` table; their contributions are shared among their
#' diseases here.
#'
#' * Additive impacts: the contribution to the aggregate (its Shapley
#'   value, including any interaction share), in the units of the impacts.
#' * Event impacts: the adjusted hazard ratio (`adjusted`), the risk of an
#'   animal with only this disease and its risk ratio and risk difference
#'   against a disease-free animal, the relative change from the raw
#'   estimate (on the raw estimate's measure), and the disease's share of the
#'   risk attributable to disease (`contribution`: Shapley allocation, as a
#'   proportion of animals).
#'
#' In both cases the contribution is split into the disease's own part
#' (`main`: its Shapley value with the interactions set to zero) and its
#' share of the interactions (`interaction`: the rest; for event impacts it
#' also holds the small changes the interaction terms make, through
#' non-additive risks, to diseases outside them). `share` is the
#' contribution as a fraction of the total.
#'
#' @param result A [deconflate()] result.
#' @return A data frame with the columns `disease`, (event impacts:
#'   `measure`), `raw`, `adjusted`, `adjusted_lower`, `adjusted_upper`,
#'   `change`, (event impacts: `risk_alone`, `rr_alone`, `rd_alone`),
#'   `contribution`, `contribution_lower`, `contribution_upper`, `main`,
#'   `interaction` and `share`.
#' @export
#' @examples
#' contribution_table(deconflate(example_supplement()))
contribution_table <- function(result) {
  if (!inherits(result, c("cm_result", "cm_event_result"))) {
    cm_abort("`result` must come from deconflate().", class = "deconflate_invalid_input")
  }
  event <- inherits(result, "cm_event_result")
  a <- result$adjusted
  ct <- result$contributions
  ids <- if (!is.null(ct)) ct$disease else a$disease[a$kind == "disease"]
  k <- match(ids, a$disease)
  col <- function(tab, nm, rows) if (is.null(tab) || is.null(tab[[nm]])) rep(NA_real_, length(ids)) else tab[[nm]][rows]
  out <- data.frame(disease = ids, stringsAsFactors = FALSE)
  if (event) out$measure <- a$measure[k]
  out$raw <- a$raw[k]
  out$adjusted <- a$adjusted[k]
  out$adjusted_lower <- col(a, "lower", k)
  out$adjusted_upper <- col(a, "upper", k)
  out$change <- a$change[k]
  if (event) {
    out$risk_alone <- a$risk_alone[k]
    out$rr_alone <- a$rr_alone[k]
    out$rd_alone <- a$rd_alone[k]
  }
  j <- if (is.null(ct)) integer(0) else match(ids, ct$disease)
  out$contribution <- col(ct, "contribution", j)
  out$contribution_lower <- col(ct, "lower", j)
  out$contribution_upper <- col(ct, "upper", j)
  out$main <- col(ct, "main", j)
  out$interaction <- col(ct, "interaction", j)
  out$share <- col(ct, "share", j)
  rownames(out) <- NULL
  out
}

#' Summarise an adjustment result
#'
#' The full report of a [deconflate()] result, in sections that each fit
#' the console: the impacts, (event impacts: the risk of an animal with only
#' these diseases), the contributions, the totals, the ranks over the draws,
#' the diagnostics and the notes. `print(result)` gives the short overview;
#' [contribution_table()] the table for export.
#'
#' @param object A [deconflate()] result.
#' @param ... Unused.
#' @return A `summary.cm_result` list with `method`, `label`, `units`,
#'   `event`, `impacts` (one row per disease and interaction: raw, adjusted,
#'   interval, change and, when some term has several rows, `rows`),
#'   `alone` (event impacts: `risk_alone`, `rr_alone`, `rd_alone`),
#'   `contributions` (per disease: `contribution` with `lower` and `upper`,
#'   `main`, `interaction`, `share`), `totals` (one row per quantity:
#'   `quantity`, `value`, `lower`, `upper`; additive impacts: raw total,
#'   adjusted total, interactions and the reduction; event impacts: overall,
#'   disease-free, attributable and raw attributable risk, the attributable
#'   fraction, the reduction and, when the Shapley allocation is incomplete,
#'   the part not allocated), `ranks` (with draws: how often each
#'   disease ranks 1st, 2nd, ... by contribution, and its mean rank; see
#'   `$ranks` of [deconflate()]), `diagnostics` (lines of text) and `notes`.
#' @export
summary.cm_result <- function(object, ...) {
  tt <- object$totals
  totals <- totals_rows(
    c("raw total", "adjusted total", "interactions", "reduction"),
    c(tt$raw_total %||% NA_real_, tt$adjusted_total, tt$interaction_total,
      reduction(tt$adjusted_total, tt$raw_total %||% NA_real_)),
    lower = c(NA, tt$adjusted_total_lower %||% NA, NA, NA),
    upper = c(NA, tt$adjusted_total_upper %||% NA, NA, NA))
  summary_object(object, totals, event = FALSE)
}

#' @rdname summary.cm_result
#' @export
summary.cm_event_result <- function(object, ...) {
  s <- object$attributable$summary
  q <- c("overall risk", "disease-free risk", "attributable risk", "attributable fraction")
  v <- c(s$overall_risk, s$disease_free_risk, s$attributable, s$attributable_fraction)
  lo <- c(NA, NA, s$attributable_lower %||% NA, NA)
  hi <- c(NA, NA, s$attributable_upper %||% NA, NA)
  if (!is.null(s$raw_attributable)) {
    q <- c(q, "raw attributable risk", "reduction")
    v <- c(v, s$raw_attributable, reduction(s$attributable, s$raw_attributable))
    lo <- c(lo, NA, NA)
    hi <- c(hi, NA, NA)
  }
  if (isTRUE(abs(s$unallocated) > max(1e-8 * abs(s$attributable), 1e-12))) {
    q <- c(q, "not allocated")
    v <- c(v, s$unallocated)
    lo <- c(lo, NA)
    hi <- c(hi, NA)
  }
  out <- summary_object(object, totals_rows(q, v, lo, hi), event = TRUE)
  out$units <- "hazard ratio; risks as proportions"
  if (identical(object$risk_given, "baseline")) {
    out$notes <- c("The disease-free risk was given; the model computed the overall risk.", out$notes)
  }
  out
}

reduction <- function(adjusted, raw) if (isTRUE(abs(raw) > 0)) 1 - adjusted / raw else NA_real_

totals_rows <- function(quantity, value, lower, upper) {
  data.frame(quantity = quantity, value = as.numeric(value), lower = as.numeric(lower), upper = as.numeric(upper),
             stringsAsFactors = FALSE)
}

summary_object <- function(object, totals, event) {
  ct <- object$contributions
  contributions <- if (!is.null(ct)) {
    out <- data.frame(disease = ct$disease, contribution = ct$contribution, stringsAsFactors = FALSE)
    if (!is.null(ct$lower)) {
      out$lower <- ct$lower
      out$upper <- ct$upper
    }
    out$main <- ct$main
    out$interaction <- ct$interaction
    out$share <- ct$share
    out
  }
  a <- object$adjusted
  alone <- if (event) data.frame(disease = a$disease, risk_alone = a$risk_alone, rr_alone = a$rr_alone,
                                 rd_alone = a$rd_alone, stringsAsFactors = FALSE)
  d <- object$diagnostics
  h <- object$heterogeneity
  diagnostics <- c(
    sprintf("Residual %s; condition number %s; sign changes %d%s.", fmt_num(d$max_reconstruction_residual),
            fmt_num(d$condition_number), d$n_sign_changes,
            if (d$n_sign_changes) sprintf(" (%s)", d$sign_changes) else ""),
    if (!is.na(d$feasibility)) sprintf("Feasibility: %s.", d$feasibility),
    if (!is.null(h) && nrow(h)) sprintf("Pooled or over-determined items: %d (see $heterogeneity and $estimates).", nrow(h)),
    unknown_pairs_line(object$unknown_pairs),
    uncertainty_line(object$draws))
  structure(list(method = object$method, label = object$label, units = object$units, event = event,
                 impacts = impact_table(object, contributions = FALSE, rows = TRUE), alone = alone,
                 contributions = contributions, totals = totals, ranks = object$ranks,
                 diagnostics = diagnostics, notes = object$notes),
            class = "summary.cm_result")
}

#' @export
print.summary.cm_result <- function(x, ...) {
  lab <- if (!is.null(x$label)) sprintf(" - %s", x$label) else ""
  un <- if (!is.null(x$units) && !is.na(x$units)) sprintf(" [%s]", x$units) else ""
  cat_wrapped(sprintf("Comorbidity adjustment (%s)%s%s", x$method, lab, un))
  cat("\nImpacts\n")
  print_table(x$impacts)
  if (!is.null(x$alone)) {
    cat("\nRisk with only these diseases (rr, rd: against a disease-free animal)\n")
    print_table(x$alone)
  }
  if (!is.null(x$contributions)) {
    cat(if (isTRUE(x$event)) "\nContributions to the attributable risk (Shapley allocation)\n" else
      "\nContributions to the adjusted total (Shapley values)\n")
    print_table(x$contributions)
  }
  cat("\nTotals\n")
  tot <- x$totals
  if (all(is.na(tot$lower))) tot$lower <- tot$upper <- NULL
  print_table(tot, left = "quantity", hide = "quantity")
  if (!is.null(x$ranks)) {
    cat("\nRank by contribution over the draws (share of draws in each place)\n")
    r <- x$ranks[order(x$ranks$mean_rank), , drop = FALSE]
    p <- as.matrix(r[, grep("^rank_", names(r)), drop = FALSE])
    tab <- data.frame(disease = r$disease, `mean rank` = sprintf("%.2f", r$mean_rank),
                      matrix(sprintf("%.2f", p), nrow = nrow(p), dimnames = list(NULL, sub("^rank_", "", colnames(p)))),
                      check.names = FALSE, stringsAsFactors = FALSE)
    print(tab, row.names = FALSE, right = TRUE)
  }
  cat("\nDiagnostics\n")
  cat_wrapped(x$diagnostics)
  print_notes(x$notes, title = "Notes")
  invisible(x)
}
