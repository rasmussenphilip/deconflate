#' Plots
#'
#' Base-graphics plots for the main result types:
#' * `plot(<cm_result>)`: raw and adjusted impacts by disease;
#'   `plot(<cm_results>)` one panel per analysis.
#' * `plot_burden()`: each disease's share of the aggregate (including its
#'   share of interaction effects); for several analyses, one bar per
#'   analysis.
#' * `plot(<cm_mc>)`: Monte Carlo means and intervals of adjusted impacts.
#' * `plot(<cm_screen>)`: the most influential scenarios from
#'   [screen_associations()], [screen_interactions()] or
#'   [screen_three_way()].
#' * `plot(<cm_oat>)`: tornado plot from [sensitivity_oat()].
#'
#' @param x,result The object to plot.
#' @param top Number of rows to show.
#' @param probs Interval bounds for Monte Carlo plots.
#' @param method For Monte Carlo runs with several methods: the method to
#'   show (default the first).
#' @param ... Passed to the underlying graphics function.
#' @return The input, invisibly.
#' @name plots
#' @examples
#' res <- deconflate(example_uk_dairy_2022())
#' plot(res$yield)
#' plot_burden(res)
NULL

impact_axis_label <- function(units) {
  if (is.null(units)) "Impact" else sprintf("Impact (%s)", units)
}

#' @rdname plots
#' @export
plot.cm_result <- function(x, ...) {
  a <- x$adjusted
  cols <- grDevices::hcl.colors(2, "Blues 3")
  mat <- rbind(raw = a$raw, adjusted = a$adjusted)
  colnames(mat) <- a$disease
  graphics::barplot(mat, beside = TRUE, col = cols, las = 2,
                    ylab = impact_axis_label(x$units),
                    main = sprintf("%s (%s)", x$label %||% "Impacts", x$method),
                    legend.text = c("raw", "adjusted"),
                    args.legend = list(bty = "n"), ...)
  graphics::abline(h = 0, col = "grey40")
  invisible(x)
}

#' @rdname plots
#' @export
plot.cm_results <- function(x, ...) {
  op <- graphics::par(mfrow = c(length(x), 1), mar = c(5, 4, 3, 1))
  on.exit(graphics::par(op))
  for (nm in names(x)) plot(x[[nm]], ...)
  invisible(x)
}

#' @rdname plots
#' @export
plot_burden <- function(result, ...) {
  if (inherits(result, "cm_result")) {
    res <- list(result)
    names(res) <- result$label %||% "impacts"
  } else if (inherits(result, "cm_results")) {
    res <- unclass(result)
  } else {
    cm_abort("`result` must come from deconflate().")
  }
  ids <- res[[1]]$adjusted$disease
  mat <- vapply(res, function(r) 100 * r$contributions$share, numeric(length(ids)))
  ylab <- "Share of the aggregate (%)"
  mat <- matrix(mat, nrow = length(ids), dimnames = list(ids, names(res)))
  cols <- grDevices::hcl.colors(nrow(mat), "Set 2")
  op <- graphics::par(mar = c(5, 4, 3, 8), xpd = TRUE)
  on.exit(graphics::par(op))
  graphics::barplot(mat, col = cols, ylab = ylab, main = "Burden by disease", ...)
  graphics::legend("topright", inset = c(-0.25, 0), legend = rownames(mat), fill = cols,
                   bty = "n", cex = 0.8)
  invisible(result)
}

#' @rdname plots
#' @export
plot.cm_mc <- function(x, probs = c(0.025, 0.975), method = NULL, ...) {
  s <- summary(x, probs = c(probs[1], 0.5, probs[2]), diagnose = FALSE)
  method <- method %||% s$method[1]
  s <- s[s$method == method, , drop = FALSE]
  lo <- s[[paste0("q", probs[1])]]
  hi <- s[[paste0("q", probs[2])]]
  mu <- s$mean
  k <- nrow(s)
  op <- graphics::par(mar = c(4, 6, 3, 1))
  on.exit(graphics::par(op))
  graphics::plot(mu, seq_len(k), xlim = range(c(lo, hi, 0), na.rm = TRUE), yaxt = "n",
                 pch = 19, xlab = sub("^Impact", "Adjusted impact", impact_axis_label(x$units)),
                 ylab = "",
                 main = sprintf("%s (%s): mean and %g%% interval (%d draws, %d rejected)",
                                x$label %||% "Impacts", method, 100 * (probs[2] - probs[1]),
                                x$n_draws, x$n_rejected), ...)
  graphics::segments(lo, seq_len(k), hi, seq_len(k))
  graphics::axis(2, at = seq_len(k), labels = s$disease, las = 1)
  graphics::abline(v = 0, lty = 2, col = "grey50")
  invisible(x)
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
                    xlab = "Change in total burden (%)", main = "Most influential scenarios", ...)
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
                 yaxt = "n", xlab = "Change in total burden", ylab = "",
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
