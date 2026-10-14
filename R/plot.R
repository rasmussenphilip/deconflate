#' Plots
#'
#' Base-graphics plots for the main result types:
#' * `plot(<cm_result>)`: raw and adjusted impacts by disease, with 95%
#'   intervals of the adjusted impacts when the result has draws.
#' * `plot(<cm_event_result>)`: adjusted hazard ratios by disease (with
#'   intervals when the result has draws), and the raw hazard and rate
#'   ratios.
#' * `plot(<cm_sensitivity>)`: tornado plot from [cm_sensitivity()]:
#'   the range of the total as each input moves from low to high.
#'
#' @param x The object to plot.
#' @param ... Passed to the underlying graphics function.
#' @return The input, invisibly.
#' @name plots
#' @examples
#' res <- deconflate(example_supplement())
#' plot(res)
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

