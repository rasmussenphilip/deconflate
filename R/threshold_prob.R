#' Probabilistic thresholds: how likely is a conclusion as an input varies?
#'
#' The probabilistic version of [cm_threshold()]: one input is set to each
#' value of a grid while all other uncertain inputs are drawn from their
#' distributions, as in the draws of [deconflate()]. At each value it gives
#' the probability that the conclusion holds, and it finds the input values
#' at which that probability crosses chosen levels (by default 0.5, 0.8 and
#' 0.95). Where [cm_threshold()] says "SCK overtakes MET when the LAM-SCK
#' odds ratio passes 2.2", this says "SCK ranks above MET in at least 80%
#' of the draws once the odds ratio passes 3.1".
#'
#' The conclusions are those of [cm_threshold()], each as a probability over
#' the draws:
#' * `"total"`: the loss of the set of diseases (the sum of their
#'   contributions; all diseases when `diseases` is left out) is below
#'   `target` (a value, or with `relative_to` a fraction of the set's
#'   central adjusted or raw value, which is then fixed for every draw);
#' * `"rank"`: the first disease of each pair ranks above the second (by
#'   contribution, or by adjusted impact with `by = "adjusted"`);
#' * `"sign"`: a disease's adjusted impact is above no effect (above 0; for
#'   event impacts an adjusted hazard ratio above 1).
#'
#' Inputs are named, or chosen automatically, and searched over their
#' default ranges, as in [cm_threshold()] (see its Search ranges).
#'
#' @section Method:
#' The draws are those of [deconflate()]: every distribution of the model,
#' the distribution of the overall or baseline risk, and the between-study
#' variation of pooled association rows, with the pooling frozen at the
#' central estimate. The varied input itself is fixed at the grid value
#' (its own distribution, if any, is not drawn). The same draws are used at
#' every grid value (common random numbers from `seed`), so the
#' probability changes only through the input, and a curve is smooth up to
#' the steps of single draws. A draw that cannot be adjusted (e.g. an
#' infeasible set of inputs) is left out at that grid value; the
#' probability is over the accepted draws, whose number is kept with the
#' curves. A crossing is found by linear interpolation between the grid
#' values on either side of the level (on the log scale for a log-spaced
#' grid). The Monte Carlo standard error of a probability p is
#' sqrt(p (1 - p) / n), at most 0.5 / sqrt(n) with n accepted draws.
#'
#' The run time is about `n_grid * n_draws` adjustments per input (4,200 by
#' default), so the automatic choice searches only the `top` inputs.
#'
#' @inheritParams cm_threshold
#' @param levels Probabilities whose crossings are reported (default 0.5,
#'   0.8 and 0.95).
#' @param n_draws Number of draws of the other uncertain inputs (default
#'   200).
#' @param n_grid Number of values of the input (default 21).
#' @param seed Seed of the draws (as for [set.seed()]); `NULL` draws one.
#' @param sampling `"random"` (default) or `"lhs"` (Latin hypercube), as in
#'   [deconflate()].
#' @param ... `pool`, or arguments of [fit_joint()] (e.g.
#'   `backend = "sampled"`), with a model.
#' @return A `cm_threshold_prob` object: `curves` (one row per input, grid
#'   value and item: `input`, `value`, `item`, `probability`, `n` (accepted
#'   draws)), `crossings` (one row per input, item and level crossed:
#'   `input`, `label`, `item`, `probability` (the level), `threshold` (the
#'   input value where the probability crosses it), `direction` (`"rises"`
#'   or `"falls"`; a level crossed next to a grid value where no draw could
#'   be adjusted is not reported), and the grid values `lower` and `upper`
#'   around it), `inputs` (as for a `cm_threshold_set`, with the probability
#'   at both ends of the range, `p_low` and `p_high`, when one item is
#'   compared, and whether the grid was log-spaced, `log_scale`),
#'   and the settings (`conclusion`, `target`, `level` (the value the loss
#'   is compared with), `reference`, `diseases`, `relative_to`, `by`,
#'   `levels`, `n_draws`, `seed`, `sampling`).
#' @seealso [cm_threshold()] for the deterministic thresholds,
#'   [cm_sensitivity()] for which inputs matter.
#' @export
#' @examples
#' m <- cm_model(example_supplement(), example_supplement()$impacts,
#'               distributions = list("impacts:1" = dist_normal(2.5, 0.5),
#'                                    "impacts:2" = dist_normal(5, 1)))
#' res <- deconflate(m, n_draws = 0)
#' # How likely is the aggregate below 2.1, as the d2-d3 odds ratio varies?
#' cm_threshold_prob(res, "pair:d2:d3", c(1, 10), conclusion = "total", target = 2.1,
#'                   n_draws = 50, n_grid = 7, seed = 1)
cm_threshold_prob <- function(model, input = NULL, range = NULL, conclusion = c("rank", "sign", "total"),
                              target = NULL, diseases = NULL, relative_to = NULL,
                              by = c("contribution", "adjusted"), levels = c(0.5, 0.8, 0.95),
                              n_draws = 200L, n_grid = 21L, seed = NULL, sampling = c("random", "lhs"),
                              top = 3, inputs = c("entered", "unknown"), ranges = NULL, sensitivity = NULL,
                              three_way = NULL, log_scale = NULL, event_model = FALSE, overall_risk = NULL,
                              baseline_risk = NULL, ...) {
  dots <- list(...)
  retired_args(dots)
  if (identical(conclusion, "change")) {
    cm_abort("conclusion = 'change' has been replaced: use conclusion = 'total', relative_to = 'adjusted' and target = 1 + the change (e.g. target = 0.9 for 10% below the current total).",
             class = "deconflate_unsupported")
  }
  sampling <- cm_match_arg(sampling)
  if (!is.numeric(levels) || !length(levels) || any(!is.finite(levels)) || any(levels <= 0 | levels >= 1)) {
    cm_abort("`levels` must be probabilities strictly between 0 and 1.", class = "deconflate_invalid_input")
  }
  levels <- sort(unique(levels))
  if (!is.numeric(n_draws) || length(n_draws) != 1L || !is.finite(n_draws) || n_draws < 10 || n_draws != round(n_draws)) {
    cm_abort("`n_draws` must be a whole number, at least 10.", class = "deconflate_invalid_input")
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 3 || n_grid != round(n_grid)) {
    cm_abort("`n_grid` must be a whole number, at least 3.", class = "deconflate_invalid_input")
  }
  check_draw_args(n_draws, seed, sampling)
  if (!is.null(log_scale) && (!is.logical(log_scale) || length(log_scale) != 1L || is.na(log_scale))) {
    cm_abort("`log_scale` must be NULL, TRUE or FALSE.", class = "deconflate_invalid_input")
  }
  if (!is.null(input) && (!is.character(input) || !length(input) || anyNA(input))) {
    cm_abort("`input` must be an input key such as \"pair:d1:d2\", several keys, or left out.",
             class = "deconflate_invalid_input")
  }
  th <- threshold_setup(model, event_model, overall_risk, baseline_risk, dots)
  spec <- threshold_spec(conclusion, target, diseases, relative_to, by, th$model)
  res <- threshold_central(th)
  if (!result_is_finite(res)) {
    cm_abort("The central result is not finite, so the draws have nothing to start from.", class = "deconflate_nonfinite")
  }
  th$res <- res
  ctx <- result_context(res)
  th$dots <- utils::modifyList(ctx$dots, th$dots[intersect(names(th$dots), joint_arg_names)])
  th$st$pool <- ctx$plan$pool %||% "random"
  sel <- threshold_inputs(th, spec, input, range, top, inputs, ranges, sensitivity, three_way)
  lv <- threshold_level(spec, ctx$model, th$st, ctx$dots, res)
  items <- threshold_items(spec, lv$level)

  # The draws: every uncertain input of deconflate(), drawn once.
  specs <- draw_specs(ctx$model, ctx$risk$dist, res$pooling$frozen)
  if (!length(specs)) {
    cm_abort("No input has a distribution, so there is nothing to draw: use cm_threshold().",
             class = "deconflate_invalid_input")
  }
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  set.seed(seed)
  U <- draw_uniforms(as.integer(n_draws), names(specs), sampling, 10L)$U
  params <- vapply(names(specs), function(k) vapply(U[, k], specs[[k]]$q, numeric(1)), numeric(n_draws))
  params <- matrix(params, nrow = n_draws, dimnames = list(NULL, names(specs)))

  curves <- list()
  crossings <- list()
  ends <- matrix(NA_real_, nrow(sel$inputs), 2L)
  scales <- rep(NA, nrow(sel$inputs))
  for (i in seq_len(nrow(sel$inputs))) {
    key <- sel$inputs$input[i]
    if (!is.finite(sel$inputs$low[i])) next
    pc <- prob_curve(ctx, th, spec, items, key, c(sel$inputs$low[i], sel$inputs$high[i]), params,
                     as.integer(n_grid), log_scale)
    if (inherits(pc, "condition")) {
      if (!is.null(input) && length(input) == 1L) stop(pc)
      sel$inputs$how[i] <- sprintf("could not be searched: %s", conditionMessage(pc))
      next
    }
    sel$inputs$value[i] <- pc$baseline %||% sel$inputs$value[i]
    scales[i] <- isTRUE(pc$log_scale)
    curves[[length(curves) + 1L]] <- data.frame(input = key, value = rep(pc$xs, nrow(items)),
                                                item = rep(items$item, each = length(pc$xs)),
                                                probability = as.vector(pc$p), n = rep(pc$n, nrow(items)),
                                                stringsAsFactors = FALSE)
    if (nrow(items) == 1L) ends[i, ] <- pc$p[c(1L, length(pc$xs)), 1L]
    for (k in seq_len(nrow(items))) {
      cr <- prob_crossings(pc$xs, pc$p[, k], levels, pc$log_scale)
      if (nrow(cr)) {
        crossings[[length(crossings) + 1L]] <- cbind(data.frame(input = key, label = sel$inputs$label[i],
                                                                item = items$item[k], stringsAsFactors = FALSE), cr)
      }
    }
  }
  curves <- if (length(curves)) do.call(rbind, curves) else
    data.frame(input = character(0), value = numeric(0), item = character(0), probability = numeric(0),
               n = integer(0), stringsAsFactors = FALSE)
  crossings <- if (length(crossings)) do.call(rbind, crossings) else
    data.frame(input = character(0), label = character(0), item = character(0), probability = numeric(0),
               threshold = numeric(0), direction = character(0), lower = numeric(0), upper = numeric(0),
               stringsAsFactors = FALSE)
  rownames(curves) <- NULL
  rownames(crossings) <- NULL
  inp <- sel$inputs
  names(inp)[names(inp) == "value"] <- "baseline"
  inp$p_low <- ends[, 1]
  inp$p_high <- ends[, 2]
  inp$log_scale <- scales
  structure(list(curves = curves, crossings = crossings, inputs = inp, items = items$item, auto = sel$auto,
                 conclusion = spec$conclusion, target = spec$target, level = lv$level, reference = lv$reference,
                 diseases = spec$diseases, relative_to = spec$relative_to, by = spec$by, levels = levels,
                 n_draws = as.integer(n_draws), seed = seed, sampling = sampling,
                 event_model = isTRUE(th$st$event_model), outcome = impact_meta(th$model$impacts)$outcome,
                 units = impact_meta(th$model$impacts)$units),
            class = "cm_threshold_prob")
}

# The probability curve of one input: at each grid value, the share of the
# accepted draws (the varied input fixed) in which each item's conclusion
# holds. Returns list(xs, p (grid x items), n (accepted draws), baseline,
# log_scale), or the condition if the input cannot be set up.
prob_curve <- function(ctx, th, spec, items, key, range, params, n_grid, log_scale) {
  model <- ctx$model
  st <- th$st
  res <- ctx$res
  tryCatch({
    check_numeric(range, "range")
    if (length(range) != 2L || range[1] >= range[2]) cm_abort("The range must be two increasing numbers.", class = "deconflate_invalid_input")
    setter <- threshold_setter(model, key, st)
    if (setter$type == "inter" && isTRUE(st$event_model) && range[1] <= 0) {
      cm_abort("An interaction of event impacts is a hazard ratio: the range must be positive.",
               class = "deconflate_invalid_input")
    }
    ls <- log_scale %||% (setter$ratio_scale && range[1] > 0)
    if (ls && range[1] <= 0) cm_abort("A log-spaced grid needs a positive range.", class = "deconflate_invalid_input")
    xs <- if (ls) exp(seq(log(range[1]), log(range[2]), length.out = n_grid)) else
      seq(range[1], range[2], length.out = n_grid)
    xs[c(1L, n_grid)] <- range
    # The varied input is not drawn.
    ek <- entered_key(key, model)
    ck <- if (is.null(ek)) NULL else if (ek == "risk") "risk" else canonical_key(ek, model)
    keep <- setdiff(colnames(params), ck)
    if (!length(keep)) {
      cm_abort(sprintf("No input other than %s has a distribution: use cm_threshold().", key),
               class = "deconflate_invalid_input")
    }
    pv <- params[, keep, drop = FALSE]
    # Inputs that only change a value keep the central pooling and row SDs;
    # those that add or replace a row (an unknown pair, an interaction, a
    # three-way association) are refitted at each grid value.
    value_setter <- grepl("^(prob|associations|impacts):", key) || key == "risk" ||
      (setter$type %in% c("assoc", "impact") && !is.null(setter$baseline_set))
    fj <- ctx$dots[intersect(names(ctx$dots), joint_arg_names)]
    event <- isTRUE(st$event_model)
    allocate <- spec$conclusion == "total" || (spec$conclusion == "rank" && spec$by == "contribution")
    joint_start <- if (!is.null(res$joint) && identical(res$joint$backend, "exact")) res$joint
    beta_start <- if (event) log(res$adjusted$adjusted)
    holds <- function(r) {
      v <- threshold_metric(r, spec, items)
      if (spec$conclusion == "total") v < 0 else v > 0
    }
    p <- matrix(NA_real_, n_grid, nrow(items), dimnames = list(NULL, items$item))
    n_ok <- integer(n_grid)
    for (g in seq_len(n_grid)) {
      m2 <- tryCatch(setter$set(xs[g]), deconflate_error = function(e) e)
      if (inherits(m2, "condition")) next
      plan <- plan_model(m2$model, st$event_model, dots = fj, pool = st$pool, risk_kind = st$risk_kind)
      if (!is.null(ctx$plan$fit_args$seed)) plan$fit_args$seed <- ctx$plan$fit_args$seed
      frozen <- list(associations = res$pooling$frozen, impacts = res$row_se)
      if (!value_setter) {
        # The input adds or replaces a row: the pooling and the impact-row
        # SDs are those of the changed model at its point values, as
        # deconflate() would fix them for that model.
        c2 <- tryCatch(withCallingHandlers(
          run_point(m2$model, plan, risk = m2$st$risk, joint = NULL, warn = FALSE, dots = fj,
                    joint_start = joint_start, start = beta_start, allocate = FALSE, raw = FALSE),
          deconflate_nonconvergence = function(w) {
            if (inherits(w, "warning")) invokeRestart("muffleWarning")
          }), deconflate_error = function(e) NULL)
        if (is.null(c2)) next
        frozen <- list(associations = c2$pooling$frozen, impacts = c2$row_se)
      }
      hit <- matrix(NA, nrow(pv), nrow(items))
      for (d in seq_len(nrow(pv))) {
        x <- pv[d, ]
        r <- tryCatch(withCallingHandlers({
          m3 <- set_inputs(m2$model, x[names(x) != "risk"])
          rk <- if ("risk" %in% names(x)) x[["risk"]] else m2$st$risk
          if (!is.null(rk) && !(is.finite(rk) && rk > 0 && rk < 1)) {
            cm_abort("The drawn risk is not a proportion.", class = "deconflate_infeasible")
          }
          run_point(m3, plan, risk = rk, joint = NULL, warn = FALSE, dots = fj, joint_start = joint_start,
                    start = beta_start, frozen = frozen, allocate = allocate, raw = FALSE)
        }, deconflate_nonconvergence = function(w) {
          if (inherits(w, "warning")) invokeRestart("muffleWarning")
        }), deconflate_error = function(e) NULL)
        if (is.null(r) || !result_is_finite(r)) next
        hit[d, ] <- holds(r)
      }
      ok <- stats::complete.cases(hit)
      n_ok[g] <- sum(ok)
      if (n_ok[g]) p[g, ] <- colMeans(hit[ok, , drop = FALSE])
    }
    list(xs = xs, p = p, n = n_ok, baseline = setter$baseline, log_scale = ls)
  }, deconflate_error = function(e) e)
}

# Where a probability curve crosses each level: linear interpolation
# between neighbouring grid values (on the log scale for a log grid).
prob_crossings <- function(xs, p, levels, log_scale) {
  out <- list()
  tx <- if (log_scale) log(xs) else xs
  for (lv in levels) {
    for (g in seq_len(length(xs) - 1L)) {
      a <- p[g]
      b <- p[g + 1L]
      if (!is.finite(a) || !is.finite(b)) next
      reached_a <- a >= lv
      reached_b <- b >= lv
      if (reached_a == reached_b) next
      t <- if (b == a) tx[g] else tx[g] + (lv - a) / (b - a) * (tx[g + 1L] - tx[g])
      out[[length(out) + 1L]] <- data.frame(probability = lv, threshold = if (log_scale) exp(t) else t,
                                            direction = if (reached_b) "rises" else "falls",
                                            lower = xs[g], upper = xs[g + 1L], stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(probability = numeric(0), threshold = numeric(0), direction = character(0),
                      lower = numeric(0), upper = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, out)
}

# The probability statement of an item: "P(the aggregate < 2.00)", "P(LAM
# ranks above SCK)", ...
prob_statement <- function(x, item) {
  spec <- list(conclusion = x$conclusion, set = if (x$conclusion == "total" && length(x$diseases)) x$diseases,
               relative_to = x$relative_to, target = x$target)
  switch(x$conclusion,
    total = sprintf("P(%s < %s)", total_words(spec, isTRUE(x$event_model)),
                    target_words(spec, x$level, x$reference)),
    rank = {
      ab <- strsplit(item, " vs ", fixed = TRUE)[[1]]
      sprintf("P(%s ranks above %s by %s)", ab[1], ab[2], x$by)
    },
    sign = sprintf("P(the adjusted value of %s is above no effect)", item))
}

#' @export
print.cm_threshold_prob <- function(x, digits = 3, ...) {
  f <- function(v) fmt_num(v, digits, na = "NA")
  inp <- x$inputs
  one <- length(x$items) == 1L
  head <- if (one) prob_statement(x, x$items) else switch(x$conclusion,
    rank = sprintf("P(one disease ranks above another by %s), for each pair", x$by),
    sign = "P(a disease's adjusted value is above no effect), for each disease")
  cat_wrapped(sprintf("<cm_threshold_prob> %s, as each input varies; the other uncertain inputs drawn (%d draws, %s sampling, seed %d)%s",
                      head, x$n_draws, if (x$sampling == "lhs") "Latin hypercube" else "random", x$seed,
                      if (isTRUE(x$auto)) "; inputs chosen by their effect in cm_sensitivity()" else ""))
  cat("\n")
  tab <- data.frame(input = inp$label, now = f(inp$baseline), from = f(inp$low), to = f(inp$high),
                    stringsAsFactors = FALSE)
  tab$now[!is.finite(inp$baseline)] <- "-"
  if (one) {
    tab$`P at from` <- f(inp$p_low)
    tab$`P at to` <- f(inp$p_high)
  }
  rest <- sum(vapply(names(tab)[-1], function(nm) max(nchar(c(nm, tab[[nm]]))) + 1L, numeric(1)))
  room <- max(16L, getOption("width", 80L) - rest - 2L)
  long <- nchar(tab$input) > room
  tab$input[long] <- paste0(substr(tab$input[long], 1L, room - 3L), "...")
  print_table(tab, left = "input")
  cat("\nLevels crossed:\n")
  for (i in seq_len(nrow(inp))) {
    k <- inp$input[i]
    if (!any(x$curves$input == k)) {
      cat(strwrap(sprintf("%s: %s", inp$label[i], inp$how[i]), width = getOption("width", 80L), indent = 2L,
                  exdent = 4L), sep = "\n")
      next
    }
    for (it in x$items) {
      cr <- x$crossings[x$crossings$input == k & x$crossings$item == it, , drop = FALSE]
      cv <- x$curves[x$curves$input == k & x$curves$item == it, , drop = FALSE]
      lab <- if (one) inp$label[i] else sprintf("%s, %s", inp$label[i], prob_statement(x, it))
      if (!nrow(cv)) {
        txt <- inp$how[i]
      } else if (!nrow(cr)) {
        rng <- range(cv$probability, na.rm = TRUE)
        txt <- sprintf("no level crossed (P from %s to %s)", f(rng[1]), f(rng[2]))
      } else {
        txt <- paste(sprintf("%s %s at %s", ifelse(cr$direction == "rises", "rises to", "falls below"),
                             as.character(cr$probability), f(cr$threshold)), collapse = "; ")
      }
      cat(strwrap(sprintf("%s: %s", lab, txt), width = getOption("width", 80L), indent = 2L, exdent = 4L),
          sep = "\n")
    }
  }
  n <- x$curves$n[x$curves$n > 0]
  if (length(n)) {
    cat_wrapped(sprintf("Each probability is over the accepted draws (%s); its Monte Carlo standard error is at most %s. Curves in $curves.",
                        if (min(n) == max(n)) min(n) else sprintf("%d to %d", min(n), max(n)), f(0.5 / sqrt(min(n)))),
                exdent = 0L)
  }
  invisible(x)
}
