#' Plots
#'
#' Plots of the results, the sensitivity analyses and the threshold
#' searches, in base graphics. Each plot draws the data frame that
#' [plot_data()] returns for the same object and `type`, and returns it
#' invisibly, so the same plot can be drawn with other tools (e.g. ggplot2).
#'
#' * `plot(<result>, type = "impacts")`: the raw and adjusted impacts of
#'   each disease and interaction, with 95% intervals of the adjusted values
#'   when the result has draws. For event impacts, the adjusted hazard
#'   ratios on a log scale; raw values are shown for hazard and rate ratios
#'   (the other measures are on another scale).
#' * `type = "contributions"`: each disease's contribution to the aggregate
#'   (event impacts: to the attributable risk), with its interval and its
#'   raw contribution (the diseases unrelated).
#' * `type = "ranks"`: how often each disease ranks 1st, 2nd, ... by
#'   contribution over the draws (needs draws).
#' * `type = "rows"`: each impact row as entered, with the central 95% of
#'   its distribution when it has one, and as the model fits it; rows the fit disagrees
#'   with (standardised residual of 2 or more) are marked.
#' * `plot(<cm_sensitivity>)`: the aggregate (or attributable risk) at each
#'   input's low and high value, for the inputs with the largest swings.
#' * `plot(<cm_threshold>)`: the compared quantity over the searched range,
#'   the line it is compared with, and the thresholds.
#' * `plot(<cm_threshold_set>)`: one panel per input, on a common y-axis.
#' * `plot(<cm_threshold_prob>)`: the probability that the conclusion holds
#'   over each input's range, and where it crosses the chosen levels.
#'
#' Titles and axes use the impact table's optional `outcome` and `units`
#' (see [cm_impacts()]); without them they use general wording, and the
#' axes show no unit.
#'
#' @param x The object to plot.
#' @param type What to plot (see above).
#' @param main Optional title, instead of the one the plot chooses.
#' @param ... Not used.
#' @return The plotted data ([plot_data()]), invisibly.
#' @name plots
#' @seealso [plot_data()] for the data of each plot.
#' @examples
#' res <- deconflate(example_supplement(), n_draws = 0)
#' plot(res)
#' plot(res, type = "contributions")
#' plot(cm_sensitivity(res))
#' plot(cm_threshold(res, "pair:d1:d2", c(0.2, 20), conclusion = "rank"))
NULL

#' Data of a plot
#'
#' The data frame that a `plot()` of the same object and `type` draws (see
#' [plots]), in long format (one row per point, bar or line segment), to
#' draw the plot with other tools.
#'
#' Results of [deconflate()] (`type`):
#' * `"impacts"`: `term`, `kind` (`"disease"` or `"interaction"`),
#'   `measure` (event impacts), `raw`, `adjusted`, `lower`, `upper`
#'   (intervals from the draws, else `NA`) and `same_scale` (whether `raw`
#'   is on the scale of `adjusted`: always for additive impacts; for event
#'   impacts, hazard and rate ratios).
#' * `"contributions"`: `disease`, `raw_contribution`, `contribution`,
#'   `lower`, `upper`, `main`, `interaction` and `share`.
#' * `"ranks"`: `disease`, `rank`, `share` (of the accepted draws in which
#'   the disease has that rank) and `mean_rank`.
#' * `"rows"`: one row per impact row: `row`, `term`, `among`,
#'   `adjusted_for`, `measure`, `study`, `value`, `lower` and `upper` (the
#'   2.5% and 97.5% quantiles of its distribution, when it has one),
#'   `fitted`, `std_residual`, `disagrees` (`|std_residual| >= 2`) and
#'   `log_scale` (ratios of event impacts, drawn on a log axis).
#'
#' [cm_sensitivity()] (`"tornado"`): `input`, `label`, `kind`, `low`,
#' `high` (the input's values), `total_low`, `total_high` (the aggregate or
#' attributable risk at them), `central` and `swing`, largest swing first.
#'
#' [cm_threshold()] (one input or a set of inputs): `"curves"`: `input`, `label`,
#' `value` (of the input), `item`, `quantity` (the compared quantity: the
#' loss of the set, the adjusted value or the difference in score),
#' `level` (what it is compared with) and `status` (of the grid point;
#' `quantity` is `NA` where it could not be evaluated); `"thresholds"`:
#' `input`, `label`, `item`, `status`, `threshold`, `lower`, `upper` and
#' `level`.
#'
#' [cm_threshold_prob()]: `"curves"`: `input`, `label`, `value`, `item`,
#' `probability` and `n` (the draws used); `"crossings"`: its `crossings`
#' table.
#'
#' @param x A [deconflate()] result, or a result of [cm_sensitivity()],
#'   [cm_threshold()] (one input or a set of inputs) or
#'   [cm_threshold_prob()].
#' @param type Which plot's data (see Details); the first is the default.
#' @param ... Not used.
#' @return A data frame.
#' @export
#' @examples
#' res <- deconflate(example_supplement(), n_draws = 0)
#' plot_data(res, "contributions")
#' # The same plot with ggplot2:
#' # d <- plot_data(res, "contributions")
#' # ggplot2::ggplot(d, ggplot2::aes(contribution, reorder(disease, contribution))) +
#' #   ggplot2::geom_col() + ggplot2::geom_point(ggplot2::aes(raw_contribution), shape = 1)
plot_data <- function(x, ...) UseMethod("plot_data")

#' @export
plot_data.default <- function(x, ...) {
  cm_abort(sprintf("plot_data() takes a result of deconflate(), cm_sensitivity(), cm_threshold() or cm_threshold_prob(), not a %s.",
                   class(x)[1]), class = "deconflate_invalid_input")
}

# Style -----------------------------------------------------------------------

# Colours: blue for the estimates, grey rings for raw values, orange for
# what is unknown or disagrees, aqua as a third series; inks for text and
# a hairline grid. The ramp orders ranks (1st darkest).
plot_pal <- list(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", grey = "#9a9993", light = "#d9d8d3",
                 ink = "#1f1f1d", ink2 = "#52514e", muted = "#8a8984", grid = "#e8e7e3", surface = "#ffffff",
                 ramp = c("#104281", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef", "#b7d3f6"))

# Graphical parameters of the plots; returns the old values.
cm_par <- function(mar = c(3.6, 4, 4.8, 1.2), ...) {
  graphics::par(bg = plot_pal$surface, fg = plot_pal$ink2, col.axis = plot_pal$ink2, col.lab = plot_pal$ink2,
                col.main = plot_pal$ink, mar = mar, mgp = c(2.1, 0.5, 0), tcl = -0.2, las = 1, cex.axis = 0.8,
                cex.lab = 0.85, xpd = FALSE, ...)
}

# The left margin (in lines) that row labels need.
label_margin <- function(labels, cex = 0.8) {
  w <- max(graphics::strwidth(labels, units = "inches", cex = cex), 0)
  max(3, w / graphics::par("csi") + 1.2)
}

# The header of a figure, in device coordinates: a bold title, a subtitle
# and a row of keys, left-aligned. Text is shrunk to fit the width.
cm_header <- function(main, sub = NULL, keys = list()) {
  op <- graphics::par(xpd = NA)
  on.exit(graphics::par(op))
  din <- graphics::par("din")
  x0 <- 0.025
  gx <- function(v) graphics::grconvertX(v, "ndc", "user")
  gy <- function(v) graphics::grconvertY(v, "ndc", "user")
  fit <- function(txt, cex) {
    w <- graphics::strwidth(txt, units = "inches", cex = cex)
    if (w > 0.95 * din[1]) cex * max(0.6, 0.95 * din[1] / w) else cex
  }
  y <- 1 - 0.3 / din[2]
  if (length(main) && !is.na(main) && nzchar(main)) {
    graphics::text(gx(x0), gy(y), main, adj = c(0, 0.5), font = 2, cex = fit(main, 0.95), col = plot_pal$ink)
    y <- y - 0.24 / din[2]
  }
  if (!is.null(sub)) {
    graphics::text(gx(x0), gy(y), sub, adj = c(0, 0.5), cex = fit(sub, 0.78), col = plot_pal$ink2)
    y <- y - 0.24 / din[2]
  }
  if (length(keys)) {
    x <- x0
    for (k in keys) {
      w <- graphics::strwidth(k$label, units = "inches", cex = 0.75) / din[1]
      if (x > x0 && x + 0.03 + w > 0.98) {
        x <- x0
        y <- y - 0.2 / din[2]
      }
      k$draw(gx(x + 0.012), gy(y))
      graphics::text(gx(x + 0.03), gy(y), k$label, adj = c(0, 0.5), cex = 0.75, col = plot_pal$ink2)
      x <- x + 0.03 + w + 0.03
    }
  }
  invisible(NULL)
}

# Tick labels without padded zeros ("0.2", "1", "20").
tick_labels <- function(at) format(at, trim = TRUE, drop0trailing = TRUE)

# The ticks of an axis within `lim` (log axes: 1-2-5 steps).
axis_ticks <- function(lim, log = FALSE) {
  at <- if (log) grDevices::axisTicks(log10(lim), log = TRUE) else pretty(lim)
  at[at >= min(lim) - 1e-12 * diff(range(lim)) & at <= max(lim) + 1e-12 * diff(range(lim))]
}

# An empty panel with one row per label: hairline grid at the x ticks and
# the row labels on the left.
cm_frame <- function(xlim, labels, xlab, log = FALSE) {
  k <- length(labels)
  if (!all(is.finite(xlim))) xlim <- c(0, 1)
  if (diff(xlim) == 0) xlim <- xlim + if (log) c(0.9, 1.1) * xlim - xlim else c(-0.5, 0.5)
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = c(0.5, k + 0.5), log = if (log) "x" else "")
  at <- axis_ticks(xlim, log)
  graphics::abline(v = at, col = plot_pal$grid, lwd = 1)
  graphics::axis(1, at = at, labels = tick_labels(at), lwd = 0, lwd.ticks = 0, line = -0.4)
  graphics::axis(2, at = seq_len(k), labels = labels, lwd = 0, lwd.ticks = 0, col.axis = plot_pal$ink, line = -0.6)
  graphics::title(xlab = xlab)
}

key_dot <- function(bg) {
  function(x, y) graphics::points(x, y, pch = 21, bg = bg, col = plot_pal$surface, cex = 1.3, lwd = 1.5)
}
key_ring <- function() {
  function(x, y) graphics::points(x, y, pch = 21, bg = plot_pal$surface, col = plot_pal$grey, cex = 1.15, lwd = 1.8)
}
key_box <- function(fill) function(x, y) graphics::points(x, y, pch = 22, bg = fill, col = fill, cex = 1.5)
key_line <- function(col) {
  function(x, y) {
    # Measured on the device, so that log axes give the same length.
    xn <- graphics::grconvertX(x, "user", "ndc")
    graphics::segments(graphics::grconvertX(xn - 0.012, "ndc", "user"), y,
                       graphics::grconvertX(xn + 0.012, "ndc", "user"), y, col = col, lwd = 2)
  }
}

# "Title: outcome", or the title alone without an outcome.
with_outcome <- function(title, outcome) if (is.null(outcome)) title else sprintf("%s: %s", title, outcome)

# "Label (units)", or the label alone without units.
with_units <- function(label, units) if (is.null(units) || is.na(units)) label else sprintf("%s (%s)", label, units)

# The colours of k ranks: the ramp, stretched when there are more ranks.
rank_colours <- function(k) {
  r <- plot_pal$ramp
  if (k <= length(r)) r[unique(round(seq(1, length(r), length.out = k)))] else grDevices::colorRampPalette(r)(k)
}

# Dark fills take white text.
text_on <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  ifelse(0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ] < 0.5, "white", plot_pal$ink)
}

ordinal <- function(i) {
  suf <- ifelse(i %% 100 %in% 11:13, "th", c("th", "st", "nd", "rd", rep("th", 6))[i %% 10 + 1])
  paste0(i, suf)
}

# Results ---------------------------------------------------------------------


#' @rdname plot_data
#' @export
plot_data.cm_result <- function(x, type = c("impacts", "contributions", "ranks", "rows"), ...) {
  type <- cm_match_arg(type)
  switch(type, impacts = impacts_data(x), contributions = contributions_data(x), ranks = ranks_data(x),
         rows = rows_data(x))
}

#' @rdname plot_data
#' @export
plot_data.cm_event_result <- plot_data.cm_result

impacts_data <- function(x) {
  a <- x$adjusted
  event <- inherits(x, "cm_event_result")
  n <- nrow(a)
  measure <- if (event) a$measure else rep(NA_character_, n)
  data.frame(term = a$disease, kind = a$kind %||% rep("disease", n), measure = measure, raw = a$raw,
             adjusted = a$adjusted, lower = a$lower %||% rep(NA_real_, n), upper = a$upper %||% rep(NA_real_, n),
             same_scale = if (event) measure %in% c("HR", "rate_ratio") else rep(TRUE, n),
             stringsAsFactors = FALSE)
}

contributions_data <- function(x) {
  ct <- x$contributions
  if (is.null(ct)) cm_abort("The result has no contributions to plot.", class = "deconflate_invalid_input")
  n <- nrow(ct)
  data.frame(disease = ct$disease, raw_contribution = ct$raw_contribution %||% rep(NA_real_, n),
             contribution = ct$contribution, lower = ct$lower %||% rep(NA_real_, n),
             upper = ct$upper %||% rep(NA_real_, n), main = ct$main, interaction = ct$interaction,
             share = ct$share, stringsAsFactors = FALSE)
}

ranks_data <- function(x) {
  r <- x$ranks
  if (is.null(r)) {
    cm_abort("The result has no draws, so no ranks: rerun deconflate() with n_draws (e.g. 500) and inputs with distributions.",
             class = "deconflate_invalid_input")
  }
  cols <- grep("^rank_", names(r), value = TRUE)
  k <- as.integer(sub("^rank_", "", cols))
  out <- data.frame(disease = rep(r$disease, each = length(cols)), rank = rep(k, nrow(r)),
                    share = as.vector(t(as.matrix(r[, cols, drop = FALSE]))),
                    mean_rank = rep(r$mean_rank, each = length(cols)), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

rows_data <- function(x) {
  e <- x$estimates
  e <- e[e$table == "impacts", , drop = FALSE]
  n <- nrow(e)
  log_scale <- inherits(x, "cm_event_result") & !is.na(e$measure) & e$measure != "RD"
  # The 2.5% and 97.5% quantiles of each row's distribution.
  q <- vapply(e$row, function(r) {
    d <- row_dist(x$model, "impacts", r)
    if (is.null(d) || identical(d$type, "fixed")) return(c(NA_real_, NA_real_))
    tryCatch(as.numeric(d$q(c(0.025, 0.975))), error = function(err) c(NA_real_, NA_real_))
  }, numeric(2))
  lo <- if (n) q[1, ] else numeric(0)
  hi <- if (n) q[2, ] else numeric(0)
  out <- data.frame(row = as.integer(e$row), term = e$item, among = e$among, adjusted_for = e$adjusted_for,
                    measure = e$measure, study = e$study, value = e$value, lower = lo, upper = hi,
                    fitted = e$fitted, std_residual = e$std_residual,
                    disagrees = !is.na(e$std_residual) & abs(e$std_residual) >= 2, log_scale = log_scale,
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' @rdname plots
#' @export
plot.cm_result <- function(x, type = c("impacts", "contributions", "ranks", "rows"), main = NULL, ...) {
  type <- cm_match_arg(type)
  d <- plot_data(x, type)
  switch(type, impacts = plot_impacts(x, d, main), contributions = plot_contributions(x, d, main),
         ranks = plot_ranks(x, d, main), rows = plot_rows(x, d, main))
  invisible(d)
}

#' @rdname plots
#' @export
plot.cm_event_result <- plot.cm_result

plot_impacts <- function(x, d, main) {
  event <- inherits(x, "cm_event_result")
  d <- d[order(d$adjusted), , drop = FALSE]
  raw <- ifelse(d$same_scale, d$raw, NA_real_)
  if (event) raw[!is.na(raw) & raw <= 0] <- NA_real_
  xr <- range(c(d$adjusted, d$lower, d$upper, raw, if (event) 1 else 0), finite = TRUE)
  op <- cm_par(mar = c(3.6, label_margin(d$term), 4.8, 1.2))
  on.exit(graphics::par(op))
  cm_frame(xr, d$term, xlab = if (event) "Hazard ratio (log scale)" else with_units("Impact", x$units), log = event)
  graphics::abline(v = if (event) 1 else 0, col = plot_pal$muted, lwd = 1)
  y <- seq_len(nrow(d))
  ok <- is.finite(d$lower) & is.finite(d$upper)
  if (any(ok)) graphics::segments(d$lower[ok], y[ok], d$upper[ok], y[ok], col = plot_pal$blue, lwd = 2, lend = 1)
  graphics::points(raw, y, pch = 21, bg = plot_pal$surface, col = plot_pal$grey, cex = 1.3, lwd = 1.8)
  graphics::points(d$adjusted, y, pch = 21, bg = plot_pal$blue, col = plot_pal$surface, cex = 1.5, lwd = 1.5)
  keys <- list(list(label = if (event) "raw (hazard and rate ratios)" else "raw", draw = key_ring()),
               list(label = if (any(ok)) "adjusted, with 95% interval" else "adjusted", draw = key_dot(plot_pal$blue)))
  cm_header(main %||% with_outcome(if (event) "Adjusted hazard ratios" else "Adjusted impacts", x$outcome),
            "Each disease and interaction, raw and adjusted for the associations", keys)
}

plot_contributions <- function(x, d, main) {
  event <- inherits(x, "cm_event_result")
  d <- d[order(d$contribution), , drop = FALSE]
  xr <- range(c(0, d$contribution, d$lower, d$upper, d$raw_contribution), finite = TRUE)
  xr[2] <- xr[2] + 0.12 * diff(xr)
  op <- cm_par(mar = c(3.6, label_margin(d$disease), 4.8, 1.2))
  on.exit(graphics::par(op))
  cm_frame(xr, d$disease, xlab = if (event) "Attributable risk" else with_units("Contribution", x$units))
  y <- seq_len(nrow(d))
  graphics::rect(0, y - 0.22, d$contribution, y + 0.22, col = plot_pal$blue, border = NA)
  graphics::abline(v = 0, col = plot_pal$muted)
  ok <- is.finite(d$lower) & is.finite(d$upper)
  if (any(ok)) graphics::segments(d$lower[ok], y[ok], d$upper[ok], y[ok], col = plot_pal$ink, lwd = 1.2)
  graphics::points(d$raw_contribution, y, pch = 21, bg = plot_pal$surface, col = plot_pal$grey, cex = 1.2, lwd = 1.8)
  right <- pmax(d$contribution, d$upper, d$raw_contribution, 0, na.rm = TRUE)
  graphics::text(right, y, fmt_num(d$contribution), pos = 4, cex = 0.72, col = plot_pal$ink2, offset = 0.5)
  tot <- if (event) x$attributable$summary$attributable else x$totals$adjusted_total
  raw <- if (event) x$attributable$summary$raw_attributable else x$totals$raw_total
  sub <- sprintf("%s %s%s", if (event) "Attributable risk" else "Adjusted total", fmt_num(tot),
                 if (isTRUE(is.finite(raw))) sprintf("; raw (diseases unrelated) %s", fmt_num(raw)) else "")
  keys <- list(list(label = if (any(ok)) "adjusted, with 95% interval" else "adjusted", draw = key_box(plot_pal$blue)))
  if (any(is.finite(d$raw_contribution))) keys <- c(keys, list(list(label = "raw (diseases unrelated)", draw = key_ring())))
  cm_header(main %||% with_outcome(if (event) "Contributions to the attributable risk" else "Contributions to the aggregate",
                                   x$outcome), sub, keys)
}

plot_ranks <- function(x, d, main) {
  ids <- unique(d$disease)
  mr <- d$mean_rank[match(ids, d$disease)]
  ids <- ids[order(-mr)]
  k <- max(d$rank)
  cols <- rank_colours(k)
  op <- cm_par(mar = c(3.6, label_margin(ids), 4.8, 1.2))
  on.exit(graphics::par(op))
  cm_frame(c(0, 1), ids, xlab = "Share of draws")
  for (i in seq_along(ids)) {
    di <- d[d$disease == ids[i], , drop = FALSE]
    di <- di[order(di$rank), , drop = FALSE]
    left <- 0
    for (j in seq_len(nrow(di))) {
      p <- di$share[j]
      fill <- cols[di$rank[j]]
      if (p > 0) graphics::rect(left, i - 0.3, left + p, i + 0.3, col = fill, border = plot_pal$surface, lwd = 2)
      if (p >= 0.08) graphics::text(left + p / 2, i, sprintf("%.2f", p), cex = 0.68, col = text_on(fill))
      left <- left + p
    }
  }
  show <- if (k <= 8L) seq_len(k) else c(1:3, k)
  keys <- lapply(show, function(j) list(label = ordinal(j), draw = key_box(cols[j])))
  dr <- x$draws
  cm_header(main %||% with_outcome("Rank by contribution", x$outcome),
            sprintf("Share of the %d accepted draws in which each disease ranks 1st, 2nd, ...",
                    as.integer(dr$n_draws - dr$n_rejected)), keys)
}

plot_rows <- function(x, d, main) {
  if (!nrow(d)) cm_abort("The result has no impact rows to plot.", class = "deconflate_invalid_input")
  lab <- sprintf("row %d: %s%s%s%s", d$row, d$term, ifelse(is.na(d$among), "", paste0(" among ", d$among)),
                 ifelse(is.na(d$adjusted_for), "", paste0(", adj. for ", d$adjusted_for)),
                 ifelse(is.na(d$measure), "", paste0(" (", d$measure, ")")))
  groups <- if (any(d$log_scale) && any(!d$log_scale)) list(which(d$log_scale), which(!d$log_scale)) else list(seq_len(nrow(d)))
  op <- cm_par(mar = c(3.6, label_margin(lab), 4.8, 1.2))
  on.exit(graphics::par(op))
  if (length(groups) > 1L) {
    graphics::layout(matrix(1:2, 2), heights = c(length(groups[[1]]) + 3.2, length(groups[[2]]) + 2))
    on.exit(graphics::layout(1), add = TRUE)
  }
  for (g in seq_along(groups)) {
    i <- rev(groups[[g]])
    lg <- all(d$log_scale[i])
    lo <- d$lower[i]
    if (lg) lo[!is.na(lo) & lo <= 0] <- NA_real_
    if (g == 2L) graphics::par(mar = c(3.6, label_margin(lab), 1, 1.2))
    xr <- range(c(lo, d$upper[i], d$value[i], d$fitted[i]), finite = TRUE)
    xlab <- if (inherits(x, "cm_event_result")) {
      if (lg) "Ratio (log scale)" else "Risk difference"
    } else with_units("Value", x$units)
    cm_frame(xr, lab[i], xlab = xlab, log = lg)
    y <- seq_along(i)
    ok <- is.finite(lo) & is.finite(d$upper[i])
    if (any(ok)) graphics::segments(lo[ok], y[ok], d$upper[i][ok], y[ok], col = plot_pal$grey, lwd = 2)
    graphics::points(d$value[i], y, pch = 21, bg = plot_pal$surface, col = plot_pal$grey, cex = 1.3, lwd = 1.8)
    graphics::points(d$fitted[i], y, pch = 21, bg = ifelse(d$disagrees[i], plot_pal$orange, plot_pal$blue),
                     col = plot_pal$surface, cex = 1.5, lwd = 1.5)
    if (g == 1L) {
      keys <- list(list(label = "entered (95% range)", draw = key_ring()),
                   list(label = "fitted", draw = key_dot(plot_pal$blue)))
      if (any(d$disagrees)) {
        keys <- c(keys, list(list(label = "fitted, disagrees (standardised residual 2 or more)",
                                  draw = key_dot(plot_pal$orange))))
      }
      cm_header(main %||% with_outcome("Impact rows", x$outcome),
                "Each row as entered and as the model fits it", keys)
    }
  }
}

# Sensitivity -------------------------------------------------------------------

#' @rdname plot_data
#' @param n For `cm_sensitivity`: the number of inputs (largest swings
#'   first).
#' @export
plot_data.cm_sensitivity <- function(x, type = "tornado", n = Inf, ...) {
  type <- cm_match_arg(type, "tornado")
  d <- x$inputs[is.finite(x$inputs$swing), , drop = FALSE]
  d <- utils::head(d, n)
  out <- data.frame(input = d$input, label = d$label, kind = d$kind, low = d$low, high = d$high,
                    total_low = d$total_low, total_high = d$total_high, central = rep(x$central$total, nrow(d)),
                    swing = d$swing, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' @rdname plots
#' @param n For `cm_sensitivity` plots: the number of inputs to show.
#' @export
plot.cm_sensitivity <- function(x, n = 10, main = NULL, ...) {
  d <- plot_data(x, n = n)
  if (!nrow(d)) cm_abort("No input has a swing to plot.", class = "deconflate_invalid_input")
  r <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  lo <- pmin(r$total_low, r$total_high)
  hi <- pmax(r$total_low, r$total_high)
  base <- x$central$total
  what <- if (isTRUE(x$event)) "attributable risk" else "aggregate"
  op <- cm_par(mar = c(3.6, label_margin(r$label), 4.8, 1.2))
  on.exit(graphics::par(op))
  cm_frame(range(lo, hi, base, finite = TRUE), r$label,
           xlab = if (isTRUE(x$event)) "Attributable risk" else with_units("Aggregate", x$units))
  y <- seq_len(nrow(r))
  graphics::rect(lo, y - 0.28, hi, y + 0.28, col = ifelse(r$kind == "unknown", plot_pal$orange, plot_pal$blue),
                 border = NA)
  graphics::abline(v = base, col = plot_pal$ink, lwd = 1)
  kinds <- intersect(c("entered", "unknown"), r$kind)
  keys <- lapply(kinds, function(k) {
    list(label = if (k == "entered") "entered inputs" else "unknown inputs",
         draw = key_box(if (k == "entered") plot_pal$blue else plot_pal$orange))
  })
  n_all <- sum(is.finite(x$inputs$swing))
  cm_header(main %||% with_outcome(sprintf("Sensitivity of the %s", what), x$outcome),
            sprintf("The %s at each input's low and high value; central value %s (black line)%s", what,
                    fmt_num(base), if (n_all > nrow(d)) sprintf("; the %d largest of %d swings", nrow(d), n_all) else ""),
            keys)
  invisible(d)
}

# Thresholds --------------------------------------------------------------------

# The name of the compared quantity, for the y-axis.
threshold_ylab <- function(x) {
  event <- isTRUE(x$event_model)
  un <- if (event) NULL else x$units
  switch(x$conclusion,
    total = with_units(tools_title(sub("^the ", "", total_words(list(set = if (length(x$diseases)) x$diseases), event))), un),
    sign = if (event) "Adjusted hazard ratio" else with_units("Adjusted impact", un),
    rank = if (x$by == "contribution") {
      if (event) "Difference in attributable risk" else with_units("Difference in contribution", un)
    } else if (event) "Difference in adjusted hazard ratio" else with_units("Difference in adjusted impact", un))
}

# Curves and thresholds of one search, on the scale of the plot: for
# "total" the loss itself (not its distance to the level), for "sign" of
# event impacts the hazard ratio.
threshold_curves <- function(s, label = s$label %||% s$input) {
  v <- s$values
  items <- colnames(v)
  off <- if (s$conclusion == "total") s$level else 0
  q <- v + off
  hr <- s$conclusion == "sign" && isTRUE(s$event_model)
  if (hr) q <- exp(v)
  level <- if (s$conclusion == "total") s$level else if (hr) 1 else 0
  out <- data.frame(input = s$input, label = label, value = rep(s$scan$value, length(items)),
                    item = rep(items, each = nrow(v)), quantity = as.vector(q), level = level,
                    status = rep(s$scan$status, length(items)), stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

threshold_points <- function(s, label = s$label %||% s$input) {
  th <- s$thresholds
  hr <- s$conclusion == "sign" && isTRUE(s$event_model)
  level <- if (s$conclusion == "total") s$level else if (hr) 1 else 0
  data.frame(input = rep(s$input, nrow(th)), label = rep(label, nrow(th)), item = th$item, status = th$status,
             threshold = th$threshold, lower = th$lower, upper = th$upper, level = rep(level, nrow(th)),
             stringsAsFactors = FALSE)
}

#' @rdname plot_data
#' @export
plot_data.cm_threshold <- function(x, type = c("curves", "thresholds"), ...) {
  type <- cm_match_arg(type)
  if (type == "curves") threshold_curves(x) else threshold_points(x)
}

#' @rdname plot_data
#' @export
plot_data.cm_threshold_set <- function(x, type = c("curves", "thresholds"), ...) {
  type <- cm_match_arg(type)
  f <- if (type == "curves") threshold_curves else threshold_points
  ok <- vapply(x$searches, inherits, logical(1), "cm_threshold")
  parts <- lapply(which(ok), function(i) f(x$searches[[i]], label = x$inputs$label[i]))
  out <- do.call(rbind, parts)
  if (is.null(out)) {
    out <- if (type == "curves") {
      data.frame(input = character(0), label = character(0), value = numeric(0), item = character(0),
                 quantity = numeric(0), level = numeric(0), status = character(0), stringsAsFactors = FALSE)
    } else {
      data.frame(input = character(0), label = character(0), item = character(0), status = character(0),
                 threshold = numeric(0), lower = numeric(0), upper = numeric(0), level = numeric(0),
                 stringsAsFactors = FALSE)
    }
  }
  rownames(out) <- NULL
  out
}

# Colours of the items of a threshold plot: one item in blue; with several,
# the items that cross in blue, orange and aqua (then darker inks), the
# others light grey.
item_colours <- function(items, crossing) {
  if (length(items) == 1L) return(stats::setNames(plot_pal$blue, items))
  hc <- c(plot_pal$blue, plot_pal$orange, plot_pal$aqua, plot_pal$ink, plot_pal$ink2, plot_pal$muted)
  cols <- stats::setNames(rep(plot_pal$light, length(items)), items)
  hit <- items[items %in% crossing]
  if (!length(hit)) hit <- utils::head(items, length(hc))
  cols[hit] <- rep(hc, length.out = length(hit))
  cols
}

# One threshold panel: the curves, the level, "now", the thresholds and
# the grid points that could not be evaluated.
threshold_panel <- function(cv, pts, cols, xlim, ylim, log_x, log_y, baseline, ylab = NULL, yaxis = TRUE,
                            xlab = NULL) {
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim, log = paste0(if (log_x) "x" else "", if (log_y) "y" else ""))
  at <- axis_ticks(xlim, log_x)
  aty <- axis_ticks(ylim, log_y)
  graphics::abline(h = aty, col = plot_pal$grid)
  graphics::axis(1, at = at, labels = tick_labels(at), lwd = 0, line = -0.4)
  if (yaxis) graphics::axis(2, at = aty, labels = tick_labels(aty), lwd = 0, line = -0.4)
  usr <- graphics::par("usr")
  top <- if (log_y) 10^usr[4] else usr[4]
  bottom <- if (log_y) 10^usr[3] else usr[3]
  if (isTRUE(is.finite(baseline)) && (!log_x || baseline > 0)) {
    graphics::abline(v = baseline, col = plot_pal$muted, lwd = 1)
    graphics::text(baseline, top, "now", adj = c(-0.15, 1.2), cex = 0.68, col = plot_pal$ink2)
  }
  level <- cv$level[1]
  if (is.finite(level) && (!log_y || level > 0)) graphics::abline(h = level, col = plot_pal$ink2, lwd = 1)
  items <- names(cols)
  grey_first <- items[order(cols != plot_pal$light, decreasing = FALSE)]
  for (it in grey_first) {
    ci <- cv[cv$item == it, , drop = FALSE]
    graphics::lines(ci$value, ci$quantity, col = cols[[it]], lwd = if (cols[[it]] == plot_pal$light) 1.5 else 2)
  }
  bad <- unique(cv$value[cv$status != "ok"])
  if (length(bad)) graphics::points(bad, rep(bottom, length(bad)), pch = 124, col = plot_pal$orange, cex = 0.8)
  th <- pts[pts$status == "threshold" & pts$item %in% items, , drop = FALSE]
  if (nrow(th)) {
    graphics::segments(th$threshold, bottom, th$threshold, th$level, col = plot_pal$ink2, lwd = 1)
    graphics::points(th$threshold, th$level, pch = 21, bg = plot_pal$ink, col = plot_pal$surface, cex = 1.3, lwd = 1.5)
    # Labels that would overlap go one line up.
    th <- th[order(th$threshold), , drop = FALSE]
    lab <- fmt_num(th$threshold)
    xin <- graphics::grconvertX(th$threshold, "user", "inches")
    w <- graphics::strwidth(lab, units = "inches", cex = 0.7) + 0.08
    lift <- integer(nrow(th))
    for (k in seq_len(nrow(th))[-1]) {
      used <- lift[seq_len(k - 1L)][xin[seq_len(k - 1L)] + w[seq_len(k - 1L)] > xin[k]]
      lift[k] <- if (length(used)) min(setdiff(0:length(used), used)) else 0L
    }
    ly <- graphics::grconvertY(graphics::grconvertY(bottom, "user", "inches") + lift * 0.16, "inches", "user")
    graphics::text(th$threshold, ly, lab, adj = c(-0.1, -0.4), cex = 0.7, col = plot_pal$ink)
  }
  if (!is.null(xlab)) graphics::title(xlab = xlab)
  if (!is.null(ylab)) graphics::title(ylab = ylab)
}

# The y-range of threshold curves, with room under them for the threshold
# labels.
threshold_ylim <- function(cv, log_y) {
  v <- c(cv$quantity, cv$level[1])
  v <- v[is.finite(v) & (!log_y | v > 0)]
  if (!length(v)) return(if (log_y) c(0.5, 2) else c(-1, 1))
  r <- range(v)
  if (log_y) {
    r <- log(r)
    if (diff(r) == 0) r <- r + c(-0.1, 0.1)
    return(exp(r + c(-0.12, 0.04) * diff(r)))
  }
  if (diff(r) == 0) r <- r + c(-0.5, 0.5) * max(abs(r), 1)
  r + c(-0.12, 0.04) * diff(r)
}

threshold_keys <- function(cols, bad, dots = TRUE) {
  keys <- list()
  if (length(cols) > 1L) {
    hit <- names(cols)[cols != plot_pal$light]
    keys <- lapply(hit, function(it) list(label = it, draw = key_line(cols[[it]])))
    if (any(cols == plot_pal$light)) {
      keys <- c(keys, list(list(label = "no threshold", draw = key_line(plot_pal$light))))
    }
  }
  if (dots) keys <- c(keys, list(list(label = "threshold", draw = key_dot(plot_pal$ink))))
  if (bad) {
    keys <- c(keys, list(list(label = "could not be evaluated", draw = function(x, y) {
      graphics::points(x, y, pch = 124, col = plot_pal$orange, cex = 0.8)
    })))
  }
  keys
}

#' @rdname plots
#' @param items For threshold plots: the items to show (default all).
#' @export
plot.cm_threshold <- function(x, items = NULL, main = NULL, ...) {
  cv <- plot_data(x, "curves")
  pts <- plot_data(x, "thresholds")
  all_items <- colnames(x$values)
  items <- check_plot_items(items, all_items)
  cv <- cv[cv$item %in% items, , drop = FALSE]
  log_y <- x$conclusion == "sign" && isTRUE(x$event_model)
  cols <- item_colours(items, pts$item[pts$status == "threshold"])
  op <- cm_par(mar = c(3.6, 4.2, 4.8, 1.2))
  on.exit(graphics::par(op))
  label <- x$label %||% x$input
  threshold_panel(cv, pts, cols, xlim = x$range, ylim = threshold_ylim(cv, log_y), log_x = isTRUE(x$log_scale),
                  log_y = log_y, baseline = x$baseline, ylab = threshold_ylab(x), xlab = label)
  cm_header(main %||% tools_title(threshold_question(x)),
            sprintf("%s from %s to %s, everything else at the model's values", label,
                    format(signif(x$range[1], 3)), format(signif(x$range[2], 3))),
            threshold_keys(cols, any(cv$status != "ok"), any(pts$status == "threshold" & pts$item %in% items)))
  out <- plot_data(x, "curves")
  invisible(out[out$item %in% items, , drop = FALSE])
}

check_plot_items <- function(items, available) {
  items <- items %||% available
  unk <- setdiff(items, available)
  if (length(unk)) {
    cm_abort(sprintf("Unknown items: %s (available: %s).", paste(unk, collapse = ", "),
                     paste(available, collapse = ", ")), class = "deconflate_invalid_input")
  }
  items
}

#' @rdname plots
#' @export
plot.cm_threshold_set <- function(x, items = NULL, main = NULL, ...) {
  cv <- plot_data(x, "curves")
  pts <- plot_data(x, "thresholds")
  n <- nrow(x$inputs)
  if (!n) cm_abort("The threshold set has no inputs to plot.", class = "deconflate_invalid_input")
  first <- Find(function(s) inherits(s, "cm_threshold"), x$searches)
  all_items <- if (!is.null(first)) colnames(first$values) else character(0)
  items <- check_plot_items(items, all_items)
  cv <- cv[cv$item %in% items, , drop = FALSE]
  log_y <- x$conclusion == "sign" && isTRUE(x$event_model)
  cols <- item_colours(items, pts$item[pts$status == "threshold"])
  ylim <- threshold_ylim(cv, log_y)
  nc <- min(n, 3L)
  nr <- ceiling(n / nc)
  op <- cm_par(mar = c(3.6, 1, 2, 0.8), mfrow = c(nr, nc), oma = c(0, 3.4, 4.4, 0))
  on.exit(graphics::par(op))
  # A grid of panels shrinks the text; keep it at full size.
  graphics::par(cex = 1)
  for (i in seq_len(n)) {
    s <- x$searches[[i]]
    key <- x$inputs$input[i]
    if (inherits(s, "cm_threshold")) {
      threshold_panel(cv[cv$input == key, , drop = FALSE], pts[pts$input == key, , drop = FALSE], cols,
                      xlim = s$range, ylim = ylim, log_x = isTRUE(s$log_scale), log_y = log_y,
                      baseline = s$baseline, yaxis = (i - 1L) %% nc == 0L, xlab = "")
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, x$inputs$result[i], cex = 0.75, col = plot_pal$ink2)
    }
    graphics::mtext(x$inputs$label[i], side = 3, line = 0.4, adj = 0, cex = 0.72, col = plot_pal$ink)
  }
  graphics::mtext(threshold_ylab(x), side = 2, line = 2, outer = TRUE, las = 0, cex = 0.72, col = plot_pal$ink2)
  cm_header(main %||% tools_title(threshold_question(x)),
            if (isTRUE(x$auto)) {
              sprintf("The %d input%s that moved it most in cm_sensitivity(), each over its range", n,
                      if (n == 1L) "" else "s")
            } else "Each input over its range, everything else at the model's values",
            threshold_keys(cols, any(cv$status != "ok"), any(pts$status == "threshold" & pts$item %in% items)))
  invisible(cv)
}

# Probabilistic thresholds --------------------------------------------------------

#' @rdname plot_data
#' @export
plot_data.cm_threshold_prob <- function(x, type = c("curves", "crossings"), ...) {
  type <- cm_match_arg(type)
  if (type == "crossings") return(x$crossings)
  cv <- x$curves
  out <- data.frame(input = cv$input, label = x$inputs$label[match(cv$input, x$inputs$input)], value = cv$value,
                    item = cv$item, probability = cv$probability, n = cv$n, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' @rdname plots
#' @export
plot.cm_threshold_prob <- function(x, main = NULL, ...) {
  cv <- plot_data(x, "curves")
  cr <- x$crossings
  ins <- unique(x$inputs$input[x$inputs$input %in% cv$input])
  if (!length(ins)) cm_abort("No input has probabilities to plot.", class = "deconflate_invalid_input")
  items <- unique(cv$item)
  cols <- item_colours(items, items)
  n <- length(ins)
  nc <- min(n, 3L)
  nr <- ceiling(n / nc)
  op <- cm_par(mar = c(3.6, 1, 2, 2), mfrow = c(nr, nc), oma = c(0, 3.4, 4.4, 0))
  on.exit(graphics::par(op))
  # A grid of panels shrinks the text; keep it at full size.
  graphics::par(cex = 1)
  for (i in seq_len(n)) {
    ci <- cv[cv$input == ins[i], , drop = FALSE]
    inp <- x$inputs[x$inputs$input == ins[i], , drop = FALSE]
    log_x <- isTRUE(inp$log_scale[1]) && all(ci$value > 0)
    xlim <- range(ci$value)
    graphics::plot.new()
    graphics::plot.window(xlim = xlim, ylim = c(0, 1), log = if (log_x) "x" else "")
    graphics::abline(h = c(0, 0.25, 0.5, 0.75, 1), col = plot_pal$grid)
    graphics::abline(h = x$levels, col = plot_pal$muted, lwd = 1)
    if (i %% nc == 0L || i == n) {
      op2 <- graphics::par(xpd = NA)
      usr <- graphics::par("usr")
      graphics::text(if (log_x) 10^usr[2] else usr[2], x$levels, format(x$levels), pos = 4, cex = 0.65,
                     col = plot_pal$ink2, offset = 0.3)
      graphics::par(op2)
    }
    at <- axis_ticks(xlim, log_x)
    graphics::axis(1, at = at, labels = tick_labels(at), lwd = 0, line = -0.4)
    if ((i - 1L) %% nc == 0L) graphics::axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), lwd = 0, line = -0.4)
    b <- inp$baseline[1]
    if (isTRUE(is.finite(b)) && (!log_x || b > 0)) {
      graphics::abline(v = b, col = plot_pal$muted)
      graphics::text(b, 1, "now", adj = c(-0.15, -0.4), cex = 0.68, col = plot_pal$ink2)
    }
    for (it in items) {
      cj <- ci[ci$item == it, , drop = FALSE]
      cj <- cj[order(cj$value), , drop = FALSE]
      graphics::lines(cj$value, cj$probability, col = cols[[it]], lwd = 2)
      graphics::points(cj$value, cj$probability, pch = 21, bg = cols[[it]], col = plot_pal$surface, cex = 0.9, lwd = 1)
    }
    ri <- cr[cr$input == ins[i] & is.finite(cr$threshold), , drop = FALSE]
    if (nrow(ri)) {
      graphics::points(ri$threshold, ri$probability, pch = 21, bg = plot_pal$ink, col = plot_pal$surface, cex = 1.3,
                       lwd = 1.5)
    }
    graphics::mtext(inp$label[1], side = 3, line = 0.4, adj = 0, cex = 0.72, col = plot_pal$ink)
  }
  graphics::mtext("Probability", side = 2, line = 2, outer = TRUE, las = 0, cex = 0.72, col = plot_pal$ink2)
  keys <- if (length(items) > 1L) lapply(items, function(it) list(label = prob_statement(x, it), draw = key_line(cols[[it]])))
  keys <- c(keys, list(list(label = sprintf("crosses %s", paste(format(x$levels), collapse = ", ")),
                            draw = key_dot(plot_pal$ink))))
  cm_header(main %||% if (length(items) == 1L) prob_statement(x, items) else "Probability that each conclusion holds",
            sprintf("As each input varies, with the other uncertain inputs drawn (%d draws)", x$n_draws), keys)
  invisible(cv)
}
