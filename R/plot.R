#' Plots
#'
#' Base-graphics plots for the main result types:
#' * `plot(<cm_result>)`: raw and adjusted impacts by disease, with 95%
#'   intervals of the adjusted impacts when the result has draws.
#' * `plot(<cm_event_result>)`: adjusted hazard ratios by disease (with
#'   intervals when the result has draws), and the raw hazard and rate
#'   ratios.
#' * `plot_burden()`: each disease's share of the aggregate (including its
#'   share of interaction effects), or of the risk attributable to disease.
#' * `plot(<cm_screen>)`: the most influential scenarios from
#'   [screen_associations()], [screen_interactions()] or
#'   [screen_three_way()].
#' * `plot(<cm_oat>)`: tornado plot from [sensitivity_oat()].
#'
#' @param x,result The object to plot.
#' @param top Number of rows to show.
#' @param ... Passed to the underlying graphics function.
#' @return The input, invisibly.
#' @name plots
#' @examples
#' res <- deconflate(example_supplement())
#' plot(res)
#' plot_burden(res)
NULL

impact_axis_label <- function(units) {
  if (is.null(units) || is.na(units)) "Impact" else sprintf("Impact (%s)", units)
}

#' @rdname plots
#' @export
plot.cm_result <- function(x, ...) {
  a <- x$adjusted
  cols <- grDevices::hcl.colors(2, "Blues 3")
  mat <- rbind(raw = a$raw, adjusted = a$adjusted)
  colnames(mat) <- a$disease
  ylim <- range(c(0, mat, a$lower, a$upper), finite = TRUE)
  mids <- graphics::barplot(mat, beside = TRUE, col = cols, las = 2, ylim = ylim,
                            ylab = impact_axis_label(x$units),
                            main = sprintf("%s (%s)", x$label %||% "Impacts", x$method),
                            legend.text = c("raw", "adjusted"),
                            args.legend = list(bty = "n"), ...)
  if (!is.null(a$lower)) {
    graphics::arrows(mids[2, ], a$lower, mids[2, ], a$upper, angle = 90, code = 3, length = 0.04)
  }
  graphics::abline(h = 0, col = "grey40")
  invisible(x)
}

#' @rdname plots
#' @export
plot.cm_event_result <- function(x, ...) {
  a <- x$adjusted
  k <- nrow(a)
  lo <- a$lower %||% a$adjusted
  hi <- a$upper %||% a$adjusted
  raw_hr <- ifelse(a$measure %in% c("HR", "rate_ratio"), a$raw, NA_real_)
  op <- graphics::par(mar = c(4, 6, 3, 1))
  on.exit(graphics::par(op))
  xl <- range(c(1, lo, hi, raw_hr), finite = TRUE)
  graphics::plot(a$adjusted, seq_len(k), xlim = xl, yaxt = "n", pch = 19, log = "x",
                 xlab = "Hazard ratio", ylab = "",
                 main = sprintf("%s: adjusted hazard ratios", x$label %||% "Event impacts"), ...)
  if (!is.null(a$lower)) graphics::segments(lo, seq_len(k), hi, seq_len(k))
  graphics::points(raw_hr, seq_len(k), pch = 1, col = "grey40")
  graphics::axis(2, at = seq_len(k), labels = a$disease, las = 1)
  graphics::abline(v = 1, lty = 2, col = "grey50")
  graphics::legend("bottomright", legend = c("adjusted", "raw (hazard and rate ratios)"),
                   pch = c(19, 1), col = c("black", "grey40"), bty = "n", cex = 0.8)
  invisible(x)
}

#' @rdname plots
#' @export
plot_burden <- function(result, ...) {
  if (inherits(result, "cm_event_result")) {
    by <- result$attributable$by_disease
    if (is.null(by)) cm_abort("The result has no allocation of the attributable risk.")
    ids <- by$disease
    share <- 100 * by$share
    main <- "Risk attributable to disease, by disease"
    ylab <- "Share of the attributable risk (%)"
  } else if (inherits(result, "cm_result")) {
    ids <- result$adjusted$disease
    share <- 100 * result$contributions$share
    main <- "Burden by disease"
    ylab <- "Share of the aggregate (%)"
  } else {
    cm_abort("`result` must come from deconflate().")
  }
  mat <- matrix(share, nrow = length(ids), dimnames = list(ids, result$label %||% "impacts"))
  cols <- grDevices::hcl.colors(nrow(mat), "Set 2")
  op <- graphics::par(mar = c(5, 4, 3, 8), xpd = TRUE)
  on.exit(graphics::par(op))
  graphics::barplot(mat, col = cols, ylab = ylab, main = main, ...)
  graphics::legend("topright", inset = c(-0.25, 0), legend = rownames(mat), fill = cols,
                   bty = "n", cex = 0.8)
  invisible(result)
}

#' @rdname plots
#' @export
plot.cm_screen <- function(x, top = 15, ...) {
  d <- x[!is.na(x$rel_change), , drop = FALSE]
  if (!nrow(d)) cm_abort("No scenario has a relative change to plot.")
  d <- utils::head(d, top)
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  op <- graphics::par(mar = c(4, 12, 3, 1))
  on.exit(graphics::par(op))
  graphics::barplot(100 * d$rel_change, horiz = TRUE, names.arg = paste(d$pair, d$scenario),
                    las = 1, col = ifelse(d$rel_change < 0, "#4477AA", "#CC6677"),
                    xlab = sprintf("Change in the %s (%%)", attr(x, "metric") %||% "total"),
                    main = "Most influential scenarios", ...)
  graphics::abline(v = 0)
  invisible(x)
}

#' @rdname plots
#' @export
plot.cm_oat <- function(x, top = 15, ...) {
  d <- utils::head(x[!is.na(x$swing), , drop = FALSE], top)
  if (!nrow(d)) cm_abort("No input has a swing to plot.")
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  base <- attr(x, "baseline_total")
  lo <- pmin(d$total_low, d$total_high) - base
  hi <- pmax(d$total_low, d$total_high) - base
  k <- nrow(d)
  op <- graphics::par(mar = c(4, 12, 3, 1))
  on.exit(graphics::par(op))
  graphics::plot(0, 0, type = "n", xlim = range(c(lo, hi, 0)), ylim = c(0.5, k + 0.5),
                 yaxt = "n", xlab = sprintf("Change in the %s", attr(x, "metric") %||% "total"), ylab = "",
                 main = "One-at-a-time sensitivity", ...)
  graphics::rect(lo, seq_len(k) - 0.4, hi, seq_len(k) + 0.4, col = "#88CCEE", border = NA)
  graphics::axis(2, at = seq_len(k), labels = d$input, las = 1, cex.axis = 0.8)
  graphics::abline(v = 0)
  invisible(x)
}

#' @export
print.cm_screen <- function(x, ...) {
  cat(sprintf("<cm_screen> %d scenarios; baseline total %.4g\n", nrow(x),
              attr(x, "baseline_total")))
  d <- x
  class(d) <- "data.frame"
  print(utils::head(d, 15), row.names = FALSE, digits = 3)
  if (nrow(x) > 15) cat(sprintf("... %d more rows\n", nrow(x) - 15))
  invisible(x)
}
