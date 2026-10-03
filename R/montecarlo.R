#' Monte Carlo propagation of input uncertainty
#'
#' Draws input sets with a sampler (usually from [cm_sampler()]) and adjusts
#' each with [deconflate()]. Every draw is a complete [cm_model()], so the
#' uncertainty of probabilities, associations and impacts is propagated
#' jointly.
#'
#' Draws that fail are rejected and counted by type (see `rejections`):
#' infeasible inputs (an impossible probability or association, or jointly
#' infeasible pairs), singular or non-identifiable systems, numerical
#' non-convergence, unsupported combinations, and non-finite results.
#' Conditioning on acceptance changes the effective input distribution, so
#' the rejection rate is part of the result and should be reported.
#'
#' @section Several methods:
#' With more than one `method`, every draw is adjusted with each method, so
#' the methods are compared on identical inputs. A draw is rejected if any
#' method fails on it. Use [compare_methods()] on the result.
#'
#' @section Batch runs:
#' With a [cm_batch_sampler()], each draw of the shared disease and
#' association inputs is used by every analysis, with the same draw
#' identifiers; each analysis draws its own impacts. A draw that fails in one
#' analysis is rejected for that analysis only. Batch runs use simple random
#' sampling.
#'
#' @section Variance reduction:
#' * `sampling = "lhs"`: Latin hypercube sampling, in `lhs_replicates`
#'   independent blocks. LHS draws are not independent, so the Monte Carlo
#'   standard error is estimated from the spread of the block means.
#' * `proposal`: importance sampling. Named inputs are drawn from proposal
#'   distributions, and each draw is weighted by the ratio of the input's own
#'   density to the proposal density (self-normalised). Each proposal must
#'   cover the support of the input's own distribution, which is checked from
#'   the distributions' supports; point masses (fixed values) cannot be
#'   importance-sampled. Self-normalised estimates have a small finite-sample
#'   bias. A defensive mixture (see [cm_suggest_proposal()]) guarantees
#'   support but not a finite variance or better precision: compare the
#'   standard errors.
#'
#' @param sampler A function of the draw index returning a [cm_model()]
#'   (typically from [cm_sampler()]), or a [cm_batch_sampler()].
#' @param n_draws Number of draws.
#' @param method Adjustment method(s) passed to [deconflate()].
#' @param seed Optional random seed, for reproducibility.
#' @param progress Logical: print progress every 10% of draws?
#' @param sampling `"random"` (default) or `"lhs"` (needs a [cm_sampler()]).
#' @param lhs_replicates Number of independent Latin hypercube blocks.
#' @param proposal Optional named list of `cm_dist` objects (importance
#'   sampling proposals), keyed as in `params` (e.g. `"impact:SCK"`). Needs a
#'   [cm_sampler()].
#' @param ... Passed to [deconflate()].
#'
#' @return A `cm_mc` object (or `cm_mc_batch` for a batch sampler) with:
#'   * `draws`: raw and adjusted impacts and contributions by draw, method
#'     and disease;
#'   * `totals`: naive and adjusted aggregate by draw and method;
#'   * `params`: one row per accepted draw with the sampled inputs;
#'   * `weights` (normalised), `log_weights`, `ess` (Kish effective sample
#'     size), `block` (LHS block of each accepted draw) and `n_blocks` (the
#'     number of LHS blocks sampled, including blocks with no accepted
#'     draw);
#'   * `n_draws`, `n_rejected`, `rejections` (draw, type, reason),
#'     `sign_changes`, `method`, `sampling`, `proposal`, `specs`, `label`
#'     and `units`.
#' @export
#' @examples
#' s <- cm_sampler(example_supplement(), impacts = list(d1 = dist_normal(2.5, 0.5)))
#' mc <- cm_monte_carlo(s, 100, method = c("published", "simultaneous"), seed = 1)
#' compare_methods(mc)
cm_monte_carlo <- function(sampler, n_draws, method = "simultaneous", seed = NULL,
                           progress = FALSE, sampling = c("random", "lhs"),
                           lhs_replicates = 10L, proposal = NULL, ...) {
  method <- unique(match.arg(method, c("simultaneous", "published", "global"), several.ok = TRUE))
  sampling <- match.arg(sampling)
  if (inherits(sampler, "cm_batch_sampler")) {
    if (sampling != "random" || !is.null(proposal)) {
      cm_abort("Batch runs use simple random sampling (no `sampling = 'lhs'` or `proposal`).",
               class = "deconflate_unsupported")
    }
    return(mc_batch(sampler, n_draws, method, seed, progress, ...))
  }
  if (!is.function(sampler)) cm_abort("`sampler` must be a function returning a cm_model().")
  specs <- attr(sampler, "specs")
  if ((sampling == "lhs" || !is.null(proposal)) && is.null(specs)) {
    cm_abort("Latin hypercube sampling and importance sampling need a sampler from cm_sampler().")
  }
  if (!is.null(proposal)) check_proposal(proposal, specs)
  if (!is.null(seed)) set.seed(seed)

  block <- NULL
  U <- NULL
  if (sampling == "lhs") {
    R <- min(as.integer(lhs_replicates), n_draws %/% 2L)
    if (R < 2L) cm_abort("Latin hypercube sampling needs at least 4 draws (2 blocks of 2).")
    block <- rep(seq_len(R), length.out = n_draws)
    block <- sort(block)
    keys <- names(specs)
    U <- matrix(NA_real_, n_draws, length(keys), dimnames = list(NULL, keys))
    for (b in seq_len(R)) {
      idx <- which(block == b)
      nb <- length(idx)
      for (k in keys) U[idx, k] <- (sample.int(nb) - stats::runif(nb)) / nb
    }
  }

  acc <- new_accumulator(n_draws)
  logw <- rep(0, n_draws)
  step <- max(1L, n_draws %/% 10L)
  for (d in seq_len(n_draws)) {
    if (progress && d %% step == 0L) message(sprintf("Draw %d of %d", d, n_draws))
    vals <- NULL
    if (!is.null(proposal)) {
      vals <- vapply(names(proposal), function(k) {
        if (!is.null(U)) proposal[[k]]$q(U[d, k]) else proposal[[k]]$r(1)
      }, numeric(1))
      logw[d] <- sum(vapply(names(proposal), function(k) {
        specs[[k]]$logd(vals[[k]]) - proposal[[k]]$logd(vals[[k]])
      }, numeric(1)))
    }
    u <- NULL
    if (!is.null(U)) {
      cols <- setdiff(colnames(U), names(proposal))
      if (length(cols)) u <- stats::setNames(as.numeric(U[d, cols]), cols)
    }
    res <- run_draw(function() {
      if (is.null(u) && is.null(vals)) sampler(d) else sampler(d, u = u, values = vals)
    }, method, ...)
    acc <- record_draw(acc, d, res, method)
  }
  finish_mc(acc, n_draws, logw, block, method, sampling, proposal, specs)
}

# Check that each proposal covers the support of the input's own distribution.
check_proposal <- function(proposal, specs) {
  if (is.null(names(proposal)) || !all(vapply(proposal, inherits, logical(1), "cm_dist"))) {
    cm_abort("`proposal` must be a named list of cm_dist objects.")
  }
  unk <- setdiff(names(proposal), names(specs))
  if (length(unk)) cm_abort(sprintf("No sampled input named: %s.", paste(unk, collapse = ", ")))
  for (k in names(proposal)) {
    tg <- specs[[k]]
    pr <- proposal[[k]]
    if (isTRUE(tg$discrete) || isTRUE(pr$discrete)) {
      cm_abort(sprintf("Importance sampling of '%s' is not supported: point masses (fixed values) cannot be reweighted against densities.", k),
               class = "deconflate_unsupported")
    }
    # The whole support must be covered, not just its end points: a mixture
    # can span the range and still leave a gap that is never sampled.
    if (!support_covers(pr, tg)) {
      cm_abort(sprintf("The proposal for '%s' has support %s, which does not cover the input's support %s; use a defensive mixture that includes the input's own distribution (see cm_suggest_proposal()).",
                       k, format_support(pr), format_support(tg)), class = "deconflate_unsupported")
    }
  }
  invisible(TRUE)
}

# Run one draw: sample a model and adjust it with every method. Returns the
# results or a classed failure.
run_draw <- function(make_model, method, ...) {
  tryCatch(
    withCallingHandlers({
      model <- make_model()
      rs <- lapply(method, function(m) deconflate(model, method = m, warn = FALSE, ...))
      bad <- !vapply(rs, result_is_finite, logical(1))
      if (any(bad)) {
        cm_abort(sprintf("Non-finite adjusted impacts (method %s).", paste(method[bad], collapse = ", ")),
                 class = "deconflate_nonfinite")
      }
      list(model = model, results = rs)
    }, deconflate_nonconvergence = function(w) {
      if (inherits(w, "warning")) invokeRestart("muffleWarning")
    }),
    deconflate_error = function(e) e
  )
}

new_accumulator <- function(n) {
  list(draws = vector("list", n), totals = vector("list", n), params = vector("list", n),
       sign = integer(n), rej_draw = integer(0), rej_type = character(0),
       rej_reason = character(0), label = NULL, units = NULL)
}

record_draw <- function(acc, d, res, method) {
  if (inherits(res, "condition")) {
    acc$rej_draw <- c(acc$rej_draw, d)
    acc$rej_type <- c(acc$rej_type, condition_type(res))
    acc$rej_reason <- c(acc$rej_reason, conditionMessage(res))
    return(acc)
  }
  acc$draws[[d]] <- do.call(rbind, lapply(seq_along(method), function(j) {
    r <- res$results[[j]]
    data.frame(draw = d, method = method[j], disease = r$adjusted$disease,
               raw = r$adjusted$raw, adjusted = r$adjusted$adjusted,
               contribution = r$contributions$total, stringsAsFactors = FALSE)
  }))
  acc$totals[[d]] <- do.call(rbind, lapply(seq_along(method), function(j) {
    r <- res$results[[j]]
    data.frame(draw = d, method = method[j], raw_sum = r$totals$raw_sum,
               adjusted_total = r$totals$adjusted_total, stringsAsFactors = FALSE)
  }))
  acc$params[[d]] <- cbind(draw = d, flatten_model(res$model))
  acc$sign[d] <- sum(vapply(res$results, function(r) r$diagnostics$n_sign_changes, numeric(1)))
  acc$label <- acc$label %||% res$results[[1]]$label
  acc$units <- acc$units %||% res$results[[1]]$units
  acc
}

finish_mc <- function(acc, n_draws, logw, block, method, sampling, proposal, specs) {
  keep <- !vapply(acc$draws, is.null, logical(1))
  lw <- logw[keep]
  w <- numeric(0)
  if (any(keep)) {
    if (any(is.nan(lw)) || any(lw == Inf)) {
      cm_abort("Importance weights are not finite; check that each proposal covers the input's distribution.",
               class = "deconflate_nonfinite")
    }
    if (all(lw == -Inf)) cm_abort("Every accepted draw has zero importance weight.",
                                  class = "deconflate_nonfinite")
    w <- exp(lw - max(lw))
    w <- w / sum(w)
  }
  structure(list(
    draws = do.call(rbind, acc$draws[keep]),
    totals = do.call(rbind, acc$totals[keep]),
    params = do.call(rbind, acc$params[keep]),
    weights = w,
    log_weights = lw,
    ess = if (length(w)) 1 / sum(w^2) else 0,
    block = if (is.null(block)) NULL else block[keep],
    # The number of replicate blocks as sampled, including blocks whose
    # draws were all rejected (they are replicates with zero weight).
    n_blocks = if (is.null(block)) NULL else length(unique(block)),
    n_draws = n_draws, n_rejected = length(acc$rej_draw),
    rejections = data.frame(draw = acc$rej_draw, type = acc$rej_type, reason = acc$rej_reason,
                            stringsAsFactors = FALSE),
    sign_changes = data.frame(draw = which(keep), n = acc$sign[keep]),
    method = method, sampling = sampling, proposal = proposal, specs = specs,
    label = acc$label, units = acc$units
  ), class = "cm_mc")
}

# Batch run: shared population draws, one accumulator per analysis.
mc_batch <- function(bs, n_draws, method, seed, progress, ...) {
  if (!is.null(seed)) set.seed(seed)
  nms <- names(bs$samplers)
  accs <- stats::setNames(lapply(nms, function(nm) new_accumulator(n_draws)), nms)
  step <- max(1L, n_draws %/% 10L)
  for (d in seq_len(n_draws)) {
    if (progress && d %% step == 0L) message(sprintf("Draw %d of %d", d, n_draws))
    pop_vals <- vapply(bs$population_keys, function(k) bs$population_specs[[k]]$r(1), numeric(1))
    for (nm in nms) {
      s <- bs$samplers[[nm]]
      res <- run_draw(function() s(d, values = pop_vals), method, ...)
      accs[[nm]] <- record_draw(accs[[nm]], d, res, method)
    }
  }
  out <- lapply(nms, function(nm) {
    finish_mc(accs[[nm]], n_draws, rep(0, n_draws), NULL, method, "random", NULL,
              attr(bs$samplers[[nm]], "specs"))
  })
  names(out) <- nms
  structure(list(analyses = out, n_draws = n_draws, population_keys = bs$population_keys),
            class = "cm_mc_batch")
}

# One-row data frame of a model's numeric inputs, keyed as in cm_sampler().
flatten_model <- function(model) {
  v <- stats::setNames(model$diseases$value, paste0("prob:", model$diseases$id))
  a <- model$associations
  if (!is.null(a) && nrow(a)) {
    v <- c(v, stats::setNames(a$value, paste0("assoc:", a$disease1, ":", a$disease2)))
  }
  i <- model$impacts
  if (!is.null(i) && nrow(i)) v <- c(v, stats::setNames(i$value, paste0("impact:", i$disease)))
  x <- model$interactions
  if (!is.null(x) && nrow(x)) {
    v <- c(v, stats::setNames(x$value, paste0("inter:", x$disease1, ":", x$disease2)))
  }
  tw <- model$three_way
  if (!is.null(tw) && nrow(tw)) {
    v <- c(v, stats::setNames(tw$ratio, paste0("three:", tw$disease1, ":", tw$disease2, ":", tw$disease3)))
  }
  as.data.frame(as.list(v), check.names = FALSE)
}

#' Reweight Monte Carlo draws (importance sampling)
#'
#' Reweights accepted draws from [cm_monte_carlo()] so that summaries reflect
#' a different input distribution, without re-running the adjustment.
#' Weights are self-normalised, so densities need only be known up to a
#' constant. [cm_scenario()] builds the log-ratio from distributions and
#' checks supports; use `cm_reweight()` directly for custom targets.
#'
#' Reweighting cannot recover rejected draws: if the target puts weight
#' where draws were rejected, results are conditional on acceptance. Always
#' check the effective sample size.
#'
#' @param mc A `cm_mc` object.
#' @param log_ratio A function taking `mc$params` (one row per accepted draw)
#'   and returning `log(target density / density of the inputs' own
#'   distributions)` for each row. Importance weights from the original run
#'   are kept.
#' @return `mc` with updated `weights` and `ess` (Kish effective sample size,
#'   `1 / sum(w^2)`).
#' @export
cm_reweight <- function(mc, log_ratio) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  lw <- log_ratio(mc$params)
  if (length(lw) != nrow(mc$params)) cm_abort("`log_ratio` must return one value per accepted draw.")
  if (length(mc$log_weights) == length(lw)) lw <- lw + mc$log_weights
  set_weights(mc, lw)
}

set_weights <- function(mc, lw) {
  if (anyNA(lw) || any(lw == Inf)) {
    cm_abort("Log weights must be finite or -Inf (zero weight); check the supports.",
             class = "deconflate_nonfinite")
  }
  if (all(lw == -Inf)) cm_abort("The target gives zero weight to every draw.",
                                class = "deconflate_nonfinite")
  w <- exp(lw - max(lw))
  w <- w / sum(w)
  mc$weights <- w
  mc$ess <- 1 / sum(w^2)
  if (mc$ess < 0.1 * length(w)) {
    cm_warn(sprintf("Effective sample size is low (%.1f of %d draws).", mc$ess, length(w)),
            class = "deconflate_low_ess")
  }
  mc
}

#' Scenario analysis by reweighting Monte Carlo draws
#'
#' Replaces the distributions of some inputs with scenario distributions and
#' reweights the existing draws accordingly (see [cm_reweight()]). Keys are
#' those of `mc$params`, e.g. `"assoc:LAM:SCK"` or `"impact:LAM"`. Each
#' scenario distribution must lie within the support of the distribution the
#' input was sampled from.
#'
#' @param mc A `cm_mc` object from a [cm_sampler()]-based run.
#' @param changes Named list of `cm_dist` objects (the scenario
#'   distributions).
#' @return The reweighted `mc` (with `ess`).
#' @export
cm_scenario <- function(mc, changes) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  specs <- mc$specs
  if (is.null(specs)) cm_abort("`mc` has no recorded distributions; use cm_sampler() or cm_reweight().")
  if (is.null(names(changes)) || !all(vapply(changes, inherits, logical(1), "cm_dist"))) {
    cm_abort("`changes` must be a named list of cm_dist objects.")
  }
  unk <- setdiff(names(changes), names(specs))
  if (length(unk)) cm_abort(sprintf("No sampled input named: %s.", paste(unk, collapse = ", ")))
  sampled <- function(k) mc$proposal[[k]] %||% specs[[k]]
  for (k in names(changes)) {
    s <- sampled(k)
    g <- changes[[k]]
    if (isTRUE(s$discrete) || isTRUE(g$discrete)) {
      cm_abort(sprintf("Cannot reweight '%s': point masses (fixed values) cannot be reweighted against densities.", k),
               class = "deconflate_unsupported")
    }
    if (!support_covers(s, g)) {
      cm_abort(sprintf("The scenario for '%s' (support %s) is not covered by the distribution it was sampled from (support %s); widen the sampling distribution (e.g. with dist_mixture()).",
                       k, format_support(g), format_support(s)),
               class = "deconflate_unsupported")
    }
  }
  keys <- union(names(mc$proposal), names(changes))
  lw <- rep(0, nrow(mc$params))
  for (k in keys) {
    x <- mc$params[[k]]
    target <- changes[[k]] %||% specs[[k]]
    lw <- lw + target$logd(x) - sampled(k)$logd(x)
  }
  set_weights(mc, lw)
}

#' Summarise Monte Carlo results
#'
#' Summarises adjusted impacts, contributions or aggregates over the accepted
#' draws, with Monte Carlo standard errors and stability checks. When some
#' estimates look unstable, a message explains why and suggests what to do
#' (see [cm_diagnose()]).
#'
#' @param object A `cm_mc` object.
#' @param what `"adjusted"` (adjusted impacts by disease), `"contribution"`
#'   (contributions to the aggregate by disease), `"total"` (naive and
#'   adjusted aggregate) or `"rejections"` (counts by type).
#' @param probs Quantiles to report.
#' @param trim Fraction trimmed from each tail for `trimmed_mean` (a
#'   different estimand from the mean).
#' @param diagnose Logical: print suggestions when estimates look unstable?
#' @param ... Unused.
#' @return A data frame with (weighted) means, SDs, the Monte Carlo standard
#'   error of the mean (`mcse`), quantiles, a trimmed mean and stability
#'   diagnostics:
#'   * `mcse`: for independent draws, the self-normalised importance-sampling
#'     estimate `sqrt(sum(w^2 (x - mean)^2))` (with equal weights,
#'     `sd / sqrt(n)`); for Latin hypercube runs, the standard error of the
#'     pooled mean as a ratio estimator over the R replicate blocks, with
#'     `S_b` the weighted sum and `W_b` the weight of block b:
#'     `sqrt(sum((S_b - mean * W_b)^2) / (R (R - 1))) / mean(W_b)`. Blocks whose
#'     draws were all rejected, or have zero weight, count in R with
#'     `S_b = W_b = 0`. With equal weights and every block present, this is
#'     the standard deviation of the block means over `sqrt(R)`;
#'   * `rel_mcse`: `mcse` relative to the absolute mean;
#'   * `tail_share`: the share of the variance contributed by the most
#'     extreme 1% of draws;
#'   * `stability`: `"ok"`, `"imprecise"` (`rel_mcse` above 5%),
#'     `"insufficient_info"` (precision cannot be assessed: fewer than two
#'     Latin hypercube blocks with positive weight, or a mean of zero with a
#'     positive standard error; `mcse` gives the absolute precision where
#'     available), `"heavy_tail"` (`tail_share` above 60%) or
#'     `"possible_pole"` (the published approximation's denominator changes
#'     sign within the sampled inputs in a way the raw impact does not
#'     explain, so the estimate has a pole inside the input distribution and
#'     its mean may not exist). A missing precision is never reported as
#'     `"ok"`.
#' @export
summary.cm_mc <- function(object, what = c("adjusted", "contribution", "total", "rejections"),
                          probs = c(0.025, 0.5, 0.975), trim = 0.05, diagnose = TRUE, ...) {
  what <- match.arg(what)
  if (what == "rejections") {
    rj <- object$rejections
    if (!nrow(rj)) return(data.frame(type = character(0), reason = character(0), n = integer(0)))
    tab <- stats::aggregate(list(n = rep(1L, nrow(rj))), rj[, c("type", "reason")], sum)
    return(tab[order(-tab$n), , drop = FALSE])
  }
  if (is.null(object$draws)) cm_abort("All draws were rejected.")
  w_by_draw <- stats::setNames(object$weights, object$params$draw)
  blk <- if (!is.null(object$block)) stats::setNames(object$block, object$params$draw) else NULL
  if (what == "total") {
    tt <- object$totals
    d <- rbind(data.frame(draw = tt$draw, method = tt$method, disease = "adjusted_total",
                          x = tt$adjusted_total, raw = NA_real_, adjusted = NA_real_,
                          stringsAsFactors = FALSE),
               data.frame(draw = tt$draw[tt$method == object$method[1]], method = "raw",
                          disease = "raw_sum", x = tt$raw_sum[tt$method == object$method[1]],
                          raw = NA_real_, adjusted = NA_real_, stringsAsFactors = FALSE))
  } else {
    d <- object$draws
    d$x <- if (what == "adjusted") d$adjusted else d$contribution
  }
  d$w <- w_by_draw[as.character(d$draw)]
  d$block <- if (is.null(blk)) NA_integer_ else blk[as.character(d$draw)]
  grp <- list(factor(d$disease, levels = unique(d$disease)),
              factor(d$method, levels = unique(d$method)))
  groups <- split(d, grp, drop = TRUE)
  rows <- lapply(groups, function(g) mc_stats(g, probs, trim, check_pole = what != "total",
                                               n_blocks = object$n_blocks))
  out <- do.call(rbind, lapply(rows, `[[`, "row"))
  info <- do.call(rbind, lapply(rows, `[[`, "info"))
  rownames(out) <- NULL
  if (what == "total") names(out)[names(out) == "disease"] <- "quantity"
  attr(out, "n_draws") <- object$n_draws
  attr(out, "n_rejected") <- object$n_rejected
  attr(out, "ess") <- object$ess
  diag <- diagnosis_table(out, info, object)
  attr(out, "diagnosis") <- diag
  if (diagnose && nrow(diag)) message(format_diagnosis(diag))
  out
}

# Weighted summary and stability statistics for one group of draws.
mc_stats <- function(g, probs, trim, check_pole = TRUE, n_blocks = NULL) {
  w <- g$w / sum(g$w)
  x <- g$x
  mu <- sum(w * x)
  s <- sqrt(sum(w * (x - mu)^2))
  q <- weighted_quantile(x, w, probs)
  lim <- weighted_quantile(x, w, c(trim, 1 - trim))
  inside <- x >= lim[1] & x <= lim[2]
  tmean <- if (any(inside)) sum(w[inside] * x[inside]) / sum(w[inside]) else NA_real_
  if (all(!is.na(g$block))) {
    # The pooled mean is a ratio estimator, mu = sum_b S_b / sum_b W_b, over
    # the R replicate blocks (S_b: weighted sum, W_b: weight of block b).
    # Its standard error is sqrt(sum_b (S_b - mu W_b)^2 / (R (R - 1))) / mean(W_b).
    # Blocks whose draws were all rejected or have zero weight are replicates
    # with S_b = W_b = 0: they count in R. With equal block weights this is
    # sd(block means) / sqrt(R). Fewer than two blocks with positive weight
    # give no precision information.
    ix_b <- split(seq_along(x), g$block)
    Wb <- vapply(ix_b, function(ix) sum(w[ix]), numeric(1))
    Sb <- vapply(ix_b, function(ix) sum(w[ix] * x[ix]), numeric(1))
    R <- max(n_blocks %||% 0L, length(ix_b))
    mcse <- if (sum(Wb > 0) >= 2L && R >= 2L) {
      sqrt(sum((Sb - mu * Wb)^2) / (R * (R - 1))) / (sum(Wb) / R)
    } else NA_real_
  } else {
    mcse <- sqrt(sum(w^2 * (x - mu)^2))
  }
  rel <- if (isTRUE(abs(mu) > 1e-12)) mcse / abs(mu) else NA_real_
  cdev <- w * (x - mu)^2
  k <- max(1L, ceiling(0.01 * length(x)))
  tail_share <- if (isTRUE(sum(cdev) > 0)) sum(sort(cdev, decreasing = TRUE)[seq_len(k)]) / sum(cdev) else 0
  # The published approximation is m^2 / (m + c). With c = 0 (no associated
  # impacts) it reduces to m and the singularity is removable. A pole lies
  # inside the sampled inputs when m + c (recovered as m^2 / adjusted) takes
  # both signs and differs in sign from m in some draws.
  pole <- NA_real_
  if (check_pole && identical(g$method[1], "published") && all(!is.na(g$raw))) {
    m <- g$raw
    a <- g$adjusted
    ok <- is.finite(m) & is.finite(a) & abs(m) > 1e-12 & abs(a) > 1e-12
    if (sum(ok) > 1) {
      den <- m[ok]^2 / a[ok]
      both <- any(den > 0) && any(den < 0)
      flip <- mean(sign(den) != sign(m[ok]))
      pole <- if (both && flip > 0) flip else 0
    }
  }
  stability <- if (!is.na(pole) && pole > 0) {
    "possible_pole"
  } else if (length(x) >= 50 && tail_share > 0.6) {
    "heavy_tail"
  } else if (is.na(mcse) || (is.na(rel) && mcse > 0)) {
    # Precision cannot be judged: too few usable replicate blocks, or a mean
    # of (nearly) zero with a positive standard error.
    "insufficient_info"
  } else if (!is.na(rel) && rel > 0.05) {
    "imprecise"
  } else {
    "ok"
  }
  row <- data.frame(disease = g$disease[1], method = g$method[1],
                    mean = mu, sd = s, mcse = mcse,
                    t(stats::setNames(q, paste0("q", probs))),
                    trimmed_mean = tmean, rel_mcse = rel, tail_share = tail_share,
                    stability = stability, check.names = FALSE, stringsAsFactors = FALSE)
  list(row = row, info = data.frame(pole = pole, n = length(x), stringsAsFactors = FALSE))
}

diagnosis_table <- function(out, info, object) {
  key <- if (!is.null(out$quantity)) out$quantity else out$disease
  bad <- out$stability != "ok"
  rows <- list()
  for (r in which(bad)) {
    st <- out$stability[r]
    if (st == "possible_pole") {
      detail <- sprintf("the published approximation divides by m + c, which changes sign within the sampled inputs (in %.1f%% of draws its sign differs from m), so the estimate has a pole inside the input distribution and its mean may not exist",
                        100 * info$pole[r])
      sugg <- "report quantiles (e.g. the median), or use method = \"simultaneous\" (exact, no division). A trimmed mean is a different estimand. If the mean does not exist, more draws or importance sampling will not make it converge."
    } else if (st == "heavy_tail") {
      detail <- sprintf("the most extreme 1%% of draws contribute %.0f%% of the variance", 100 * out$tail_share[r])
      sugg <- sprintf("check which inputs produce the extreme draws; importance sampling may reduce the error if they come from one input region: prop <- cm_suggest_proposal(mc, \"%s\", method = \"%s\"), then cm_monte_carlo(sampler, n_draws, method = \"%s\", proposal = prop), and compare the standard errors. Report quantiles as well.",
                      key[r], out$method[r], out$method[r])
    } else if (st == "insufficient_info") {
      detail <- if (is.na(out$mcse[r])) {
        "precision cannot be assessed: fewer than two Latin hypercube blocks have draws with positive weight"
      } else {
        sprintf("the mean is (nearly) zero, so relative precision is undefined; the absolute Monte Carlo standard error is %.3g",
                out$mcse[r])
      }
      sugg <- "judge the absolute standard error against the size of effect that matters, increase n_draws or lhs_replicates, or check why draws are rejected or have zero weight (summary(mc, what = \"rejections\"))."
    } else {
      need <- ceiling(info$n[r] * (out$rel_mcse[r] / 0.02)^2)
      detail <- sprintf("the Monte Carlo standard error is %.1f%% of the mean", 100 * out$rel_mcse[r])
      sugg <- sprintf("increase n_draws to about %s for a 2%% standard error, or use sampling = \"lhs\" (Latin hypercube).",
                      format(need, big.mark = ",", scientific = FALSE))
    }
    rows[[length(rows) + 1L]] <- data.frame(item = key[r], method = out$method[r], stability = st,
                                            detail = detail, suggestion = sugg,
                                            stringsAsFactors = FALSE)
  }
  nf <- object$rejections[object$rejections$type == "nonfinite", , drop = FALSE]
  if (nrow(nf)) {
    rows[[length(rows) + 1L]] <- data.frame(
      item = "all", method = paste(object$method, collapse = ", "), stability = "non_finite",
      detail = sprintf("%d draw(s) gave non-finite results and were rejected", nrow(nf)),
      suggestion = "check the published approximation for division by zero, or near-singular inputs; the exact method avoids the division.",
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) {
    return(data.frame(item = character(0), method = character(0), stability = character(0),
                      detail = character(0), suggestion = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

format_diagnosis <- function(diag) {
  lines <- sprintf("* %s (%s): %s.\n    Suggestion: %s", diag$item, diag$method, diag$detail,
                   diag$suggestion)
  paste0(sprintf("%d Monte Carlo estimate(s) may be unstable:\n", nrow(diag)),
         paste(lines, collapse = "\n"),
         "\n(See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)")
}

#' Diagnose unstable Monte Carlo estimates
#'
#' Lists the estimates that [summary.cm_mc()] flags, with the reason and a
#' suggested remedy:
#' * `"possible_pole"`: with the published approximation, `m^2 / (m + c)`
#'   has a pole where `m + c = 0` and `m != 0`. When the sampled inputs reach
#'   both sides of it, the estimate's distribution can have tails so heavy
#'   that its mean does not exist; no sampling scheme then makes the mean
#'   converge. Report quantiles, or use the exact (`"simultaneous"`) method.
#'   (With `c = 0`, e.g. independent diseases, the formula reduces to `m` and
#'   is not flagged.)
#' * `"heavy_tail"`: a few extreme draws dominate the variance. Importance
#'   sampling may help if they come from one region of one input. The exact
#'   method can also be heavy-tailed when sampled inputs make the conflation
#'   matrix nearly singular.
#' * `"imprecise"`: the Monte Carlo standard error is large relative to the
#'   mean. Use more draws, or Latin hypercube sampling.
#' * `"insufficient_info"`: precision cannot be assessed (too few usable
#'   Latin hypercube blocks, or a mean of zero); judge the absolute standard
#'   error instead.
#' * `"non_finite"`: draws gave non-finite results and were rejected.
#'
#' @param mc A `cm_mc` object.
#' @param what Passed to [summary.cm_mc()].
#' @return A data frame (class `cm_diagnosis`) with one row per flagged
#'   estimate: item, method, stability, detail and suggestion.
#' @export
cm_diagnose <- function(mc, what = c("adjusted", "contribution", "total")) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  what <- match.arg(what)
  s <- summary(mc, what = what, diagnose = FALSE)
  out <- attr(s, "diagnosis")
  class(out) <- c("cm_diagnosis", "data.frame")
  out
}

#' @export
print.cm_diagnosis <- function(x, ...) {
  if (!nrow(x)) {
    cat("No unstable Monte Carlo estimates.\n")
  } else {
    d <- x
    class(d) <- "data.frame"
    cat(format_diagnosis(d), "\n")
  }
  invisible(x)
}

#' Suggest an importance-sampling proposal for an unstable estimate
#'
#' Finds the sampled input most associated with the most extreme draws of
#' one adjusted impact, and builds a defensive mixture proposal for it: part
#' the input's own distribution and part a uniform distribution over the
#' range of that input in the extreme draws (within the input's support).
#' The mixture always covers the input's support and bounds the weights by
#' `1 / (1 - weight)`. It may reduce the Monte Carlo error when the extreme
#' values come from that region and the mean exists; it does not guarantee a
#' finite variance or better precision, so compare standard errors.
#'
#' @param mc A `cm_mc` object from a [cm_sampler()]-based run.
#' @param disease The disease whose adjusted impact is unstable.
#' @param method Adjustment method (default: the first in `mc`).
#' @param top Fraction of draws treated as extreme.
#' @param weight Weight of the uniform component in the mixture.
#' @return A named list with one `cm_dist`, for the `proposal` argument of
#'   [cm_monte_carlo()]. The attribute `"explanation"` describes it.
#' @export
cm_suggest_proposal <- function(mc, disease, method = NULL, top = 0.02, weight = 0.5) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  if (is.null(mc$specs)) cm_abort("`mc` has no recorded distributions; use cm_sampler().")
  if (!(weight > 0 && weight < 1)) cm_abort("`weight` must be between 0 and 1.")
  method <- method %||% mc$method[1]
  d <- mc$draws[mc$draws$disease == disease & mc$draws$method == method, , drop = FALSE]
  if (nrow(d) < 20L) cm_abort("Too few accepted draws for this estimate.")
  keys <- names(mc$specs)
  keys <- keys[vapply(keys, function(k) !isTRUE(mc$specs[[k]]$discrete), logical(1))]
  if (!length(keys)) cm_abort("No sampled inputs to build a proposal for.")
  dev <- abs(d$adjusted - stats::median(d$adjusted))
  k <- max(5L, ceiling(top * nrow(d)))
  top_draws <- d$draw[order(-dev)][seq_len(min(k, nrow(d)))]
  p <- mc$params
  ks <- vapply(keys, function(key) {
    all_v <- p[[key]]
    top_v <- p[[key]][p$draw %in% top_draws]
    grid <- sort(unique(all_v))
    max(abs(stats::ecdf(top_v)(grid) - stats::ecdf(all_v)(grid)))
  }, numeric(1))
  key <- keys[which.max(ks)]
  spec <- mc$specs[[key]]
  top_v <- p[[key]][p$draw %in% top_draws]
  all_v <- p[[key]]
  span <- diff(range(all_v))
  if (span <= 0) span <- abs(mean(all_v)) + 1
  lo <- max(min(top_v) - 0.1 * span, spec$lower %||% -Inf)
  hi <- min(max(top_v) + 0.1 * span, spec$upper %||% Inf)
  if (!(hi > lo)) cm_abort("Could not find a range for the proposal.")
  prop <- dist_mixture(spec, dist_uniform(lo, hi), weights = c(1 - weight, weight))
  out <- stats::setNames(list(prop), key)
  attr(out, "explanation") <- sprintf(
    "Proposal for %s: %.0f%% its own distribution, %.0f%% uniform on [%.4g, %.4g], the range of %s in the %d most extreme draws of %s (%s).",
    key, 100 * (1 - weight), 100 * weight, lo, hi, key, length(top_draws), disease, method)
  message(attr(out, "explanation"))
  out
}

#' Productivity gaps over Monte Carlo draws (optional helper)
#'
#' Applies [productivity_gap()] (and a unit value) to every accepted draw of
#' a Monte Carlo run, and summarises the gap and its value.
#'
#' @param mc A `cm_mc` object.
#' @param observed Observed mean of the outcome (a number or a `cm_dist`,
#'   drawn per draw).
#' @param direction,effect As in [productivity_gap()].
#' @param unit_value Optional value per unit of gap (a number or a `cm_dist`,
#'   drawn per draw).
#' @param seed Optional seed for drawing `observed` and `unit_value`.
#' @return A list with `draws` (gap and value by draw, method and disease,
#'   and the totals) and `summary` (weighted mean and quantiles of the total
#'   gap and value per method).
#' @export
cm_mc_gap <- function(mc, observed, direction = c("decrease", "increase"),
                      effect = c("proportion", "percent", "absolute"), unit_value = NULL,
                      seed = NULL) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  direction <- match.arg(direction)
  effect <- match.arg(effect)
  given <- list(observed = observed, unit_value = unit_value)
  for (nm in names(given)) {
    v <- given[[nm]]
    if (is.null(v) || inherits(v, "cm_dist")) next
    check_numeric(v, nm)
    if (length(v) != 1L) cm_abort(sprintf("`%s` must be a single number or a cm_dist.", nm))
  }
  if (!is.null(seed)) set.seed(seed)
  draws <- unique(mc$params$draw)
  pick <- function(v) if (inherits(v, "cm_dist")) v$r(length(draws)) else rep(as.numeric(v), length(draws))
  obs <- stats::setNames(pick(observed), draws)
  uv <- if (is.null(unit_value)) NULL else stats::setNames(pick(unit_value), draws)
  scale <- if (effect == "percent") 100 else 1
  d <- mc$draws
  tt <- mc$totals
  L <- tt$adjusted_total / scale
  x <- obs[as.character(tt$draw)]
  factor <- switch(effect, absolute = rep(1, length(L)),
                   if (direction == "decrease") x / (1 - L) else x / (1 + L))
  if (effect != "absolute" && any(if (direction == "decrease") L >= 1 else L <= -1)) {
    cm_abort("Some draws give an aggregate proportional change outside the valid range.")
  }
  tt$gap <- unname(if (effect == "absolute") L else if (direction == "decrease") x / (1 - L) - x else x - x / (1 + L))
  fkey <- paste(tt$draw, tt$method)
  d$gap <- unname(factor[match(paste(d$draw, d$method), fkey)] * d$contribution / scale)
  if (!is.null(uv)) {
    tt$value <- unname(tt$gap * uv[as.character(tt$draw)])
    d$value <- unname(d$gap * uv[as.character(d$draw)])
  }
  w <- stats::setNames(mc$weights, mc$params$draw)
  summ <- do.call(rbind, lapply(split(tt, tt$method), function(g) {
    ww <- w[as.character(g$draw)]
    ww <- ww / sum(ww)
    out <- data.frame(method = g$method[1], gap_mean = sum(ww * g$gap),
                      gap_q0.025 = weighted_quantile(g$gap, ww, 0.025),
                      gap_q0.975 = weighted_quantile(g$gap, ww, 0.975),
                      stringsAsFactors = FALSE)
    if (!is.null(g$value)) {
      out$value_mean <- sum(ww * g$value)
      out$value_q0.025 <- weighted_quantile(g$value, ww, 0.025)
      out$value_q0.975 <- weighted_quantile(g$value, ww, 0.975)
    }
    out
  }))
  rownames(summ) <- NULL
  list(draws = d, totals = tt, summary = summ)
}

#' @export
summary.cm_mc_batch <- function(object, ...) {
  parts <- lapply(names(object$analyses), function(nm) {
    a <- object$analyses[[nm]]
    s <- tryCatch(summary(a, ...), deconflate_error = function(e) NULL)
    if (is.null(s) || !nrow(s)) return(NULL)
    cbind(analysis = nm, s, stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(data.frame(analysis = character(0), stringsAsFactors = FALSE))
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' @export
print.cm_mc_batch <- function(x, ...) {
  cat(sprintf("<cm_mc_batch> %d draws, %d analyses\n", x$n_draws, length(x$analyses)))
  for (nm in names(x$analyses)) {
    a <- x$analyses[[nm]]
    cat(sprintf("  %s: %d accepted, %d rejected\n", nm, x$n_draws - a$n_rejected, a$n_rejected))
  }
  invisible(x)
}
