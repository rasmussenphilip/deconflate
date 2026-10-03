#' Plots
#'
#' Base-graphics plots for the main result types:
#' * `plot(<cm_result>)`: raw and adjusted impacts by disease, one panel per
#'   outcome.
#' * `plot_burden()`: each disease's share of the burden per outcome (main
#'   effects and interaction shares), or monetary losses if `economics` is
#'   given.
#' * `plot(<cm_mc>)`: Monte Carlo means and intervals of adjusted impacts.
#' * `plot(<cm_screen>)`: the most influential pairs from
#'   [screen_associations()] or [screen_interactions()].
#' * `plot(<cm_oat>)`: tornado plot from [sensitivity_oat()].
#'
#' @param x,result The object to plot.
#' @param outcome Outcome(s) to show (default all, or the first for
#'   Monte Carlo results).
#' @param economics Optional economics list (see [value_losses()]).
#' @param top Number of rows to show.
#' @param probs Interval bounds for Monte Carlo plots.
#' @param method For Monte Carlo runs with several methods: the method to
#'   show (default the first).
#' @param ... Passed to the underlying graphics function.
#' @return The input, invisibly.
#' @name plots
#' @examples
#' res <- deconflate(example_uk_dairy_2022(), method = "published")
#' plot(res, outcome = "yield")
#' plot_burden(res)
NULL

#' @rdname plots
#' @export
plot.cm_result <- function(x, outcome = NULL, ...) {
  a <- x$adjusted
  outs <- outcome %||% unique(a$outcome)
  op <- graphics::par(mfrow = c(length(outs), 1), mar = c(5, 4, 3, 1))
  on.exit(graphics::par(op))
  cols <- grDevices::hcl.colors(2, "Blues 3")
  for (o in outs) {
    d <- a[a$outcome == o, , drop = FALSE]
    mult <- if (d$scale[1] == "proportion") 100 else 1
    mat <- rbind(raw = d$raw, adjusted = d$adjusted) * mult
    colnames(mat) <- d$disease
    graphics::barplot(mat, beside = TRUE, col = cols, las = 2,
                      ylab = if (mult == 100) "Impact (%)" else if (d$scale[1] == "hazard_ratio") "Hazard ratio" else "Impact",
                      main = sprintf("%s (%s)", o, x$method),
                      legend.text = c("raw", "adjusted"),
                      args.legend = list(bty = "n"), ...)
    graphics::abline(h = 0, col = "grey40")
  }
  invisible(x)
}

#' @rdname plots
#' @export
plot_burden <- function(result, economics = NULL, ...) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  ids <- result$model$diseases$id
  if (!is.null(economics)) {
    vl <- value_losses(productivity_gap(result, economics$observed), economics$unit_value)
    d <- vl$by_disease
    mat <- tapply(d$value, list(d$disease, d$outcome), sum)
    mat[is.na(mat)] <- 0
    mat <- mat[intersect(ids, rownames(mat)), , drop = FALSE]
    ylab <- "Value"
  } else {
    b <- attribute_burden(result)
    mat <- tapply(b$share, list(b$disease, b$outcome), sum)
    mat <- mat[ids, , drop = FALSE] * 100
    ylab <- "Share of burden (%)"
  }
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
plot.cm_mc <- function(x, outcome = NULL, probs = c(0.025, 0.975), method = NULL, ...) {
  s <- summary(x, probs = c(probs[1], 0.5, probs[2]), diagnose = FALSE)
  outcome <- outcome %||% s$outcome[1]
  method <- method %||% s$method[1]
  s <- s[s$outcome == outcome & s$method == method, , drop = FALSE]
  sc <- x$draws$scale[x$draws$outcome == outcome][1]
  mult <- if (is.null(sc) || is.na(sc) || sc == "proportion") 100 else 1
  ref <- if (!is.null(sc) && !is.na(sc) && sc == "hazard_ratio") 1 else 0
  lo <- s[[paste0("q", probs[1])]] * mult
  hi <- s[[paste0("q", probs[2])]] * mult
  mu <- s$mean * mult
  k <- nrow(s)
  op <- graphics::par(mar = c(4, 6, 3, 1))
  on.exit(graphics::par(op))
  graphics::plot(mu, seq_len(k), xlim = range(c(lo, hi, ref), na.rm = TRUE), yaxt = "n",
                 pch = 19,
                 xlab = if (mult == 100) "Adjusted impact (%)" else if (ref == 1) "Adjusted hazard ratio" else "Adjusted impact",
                 ylab = "",
                 main = sprintf("%s (%s): mean and %g%% interval (%d draws, %d rejected)",
                                outcome, method, 100 * (probs[2] - probs[1]),
                                x$n_draws, x$n_rejected), ...)
  graphics::segments(lo, seq_len(k), hi, seq_len(k))
  graphics::axis(2, at = seq_len(k), labels = s$disease, las = 1)
  graphics::abline(v = ref, lty = 2, col = "grey50")
  invisible(x)
}

#' @rdname plots
#' @export
plot.cm_screen <- function(x, top = 15, ...) {
  d <- x[!is.na(x$rel_change), , drop = FALSE]
  d <- utils::head(d, top)
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  op <- graphics::par(mar = c(4, 12, 3, 1))
  on.exit(graphics::par(op))
  graphics::barplot(100 * d$rel_change, horiz = TRUE, names.arg = paste(d$pair, d$scenario),
                    las = 1, col = ifelse(d$rel_change < 0, "#4477AA", "#CC6677"),
                    xlab = "Change in total burden (%)", main = "Most influential pairs", ...)
  graphics::abline(v = 0)
  invisible(x)
}

#' @rdname plots
#' @export
plot.cm_oat <- function(x, top = 15, ...) {
  d <- utils::head(x[!is.na(x$swing), , drop = FALSE], top)
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
