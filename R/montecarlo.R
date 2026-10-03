#' Monte Carlo propagation of input uncertainty
#'
#' Draws input sets with a sampler (usually from [cm_sampler()]) and adjusts
#' each with [deconflate()]. Because every draw is a complete [cm_model()],
#' disease probabilities and associations are shared across all outcomes
#' within a draw, so the outcomes' uncertainty stays dependent.
#'
#' Draws whose inputs are infeasible (an impossible probability or odds ratio,
#' an association incompatible with the sampled marginals, or a jointly
#' infeasible set for the global method) are rejected and counted.
#' Conditioning on feasibility changes the effective input distribution, so
#' the rejection rate is part of the result and should be reported.
#'
#' @param sampler A function of one argument (the draw index) returning a
#'   [cm_model()], typically from [cm_sampler()].
#' @param n_draws Number of draws.
#' @param method Adjustment method passed to [deconflate()].
#' @param economics Optional list with named elements `observed` and
#'   `unit_value` (and optionally `additional`), as in [value_losses()].
#'   Elements may be numbers or `cm_dist` objects (drawn per draw). When
#'   given, productivity gaps and monetary losses are computed for every
#'   draw.
#' @param seed Optional random seed, for reproducibility.
#' @param progress Logical: print progress every 10% of draws?
#' @param ... Passed to [deconflate()].
#'
#' @return A `cm_mc` object with elements:
#'   * `draws`: long data frame of raw and adjusted impacts by draw;
#'   * `losses`: long data frame of gaps and values by draw, outcome and
#'     disease (if `economics` was given);
#'   * `params`: one row per accepted draw with the sampled inputs, keyed
#'     `prob:<id>`, `assoc:<d1>:<d2>`, `impact:<outcome>:<disease>`,
#'     `inter:<outcome>:<d1>:<d2>`, `observed:<outcome>`, `unit_value:<outcome>`;
#'   * `weights` (initially equal), `n_draws`, `n_rejected`, `rejections`
#'     (draw and reason), `sign_changes` and `method`.
#' @export
cm_monte_carlo <- function(sampler, n_draws, method = "simultaneous",
                           economics = NULL, seed = NULL, progress = FALSE, ...) {
  if (!is.function(sampler)) cm_abort("`sampler` must be a function returning a cm_model().")
  if (!is.null(seed)) set.seed(seed)
  draws <- vector("list", n_draws)
  losses <- vector("list", n_draws)
  params <- vector("list", n_draws)
  sign_changes <- integer(n_draws)
  rej_draw <- integer(0)
  rej_reason <- character(0)
  step <- max(1L, n_draws %/% 10L)

  for (d in seq_len(n_draws)) {
    if (progress && d %% step == 0L) message(sprintf("Draw %d of %d", d, n_draws))
    res <- tryCatch(
      withCallingHandlers({
        model <- sampler(d)
        r <- deconflate(model, method = method, warn = FALSE, ...)
        eco <- if (!is.null(economics)) draw_economics(economics) else NULL
        lv <- if (!is.null(eco)) {
          value_losses(productivity_gap(r, eco$observed), eco$unit_value, eco$additional)
        } else NULL
        list(model = model, result = r, eco = eco, lv = lv)
      }, deconflate_nonconvergence = function(w) invokeRestart("muffleWarning")),
      deconflate_infeasible = function(e) e
    )
    if (inherits(res, "condition")) {
      rej_draw <- c(rej_draw, d)
      rej_reason <- c(rej_reason, conditionMessage(res))
      next
    }
    draws[[d]] <- cbind(draw = d, res$result$adjusted[, c("outcome", "disease", "raw", "adjusted")])
    p <- flatten_model(res$model)
    if (!is.null(res$eco)) {
      p <- cbind(p, as.data.frame(as.list(c(
        stats::setNames(res$eco$observed, paste0("observed:", names(res$eco$observed))),
        stats::setNames(res$eco$unit_value, paste0("unit_value:", names(res$eco$unit_value))))),
        check.names = FALSE))
      losses[[d]] <- cbind(draw = d, res$lv$by_disease)
    }
    params[[d]] <- cbind(draw = d, p)
    sign_changes[d] <- sum(res$result$diagnostics$n_sign_changes)
  }
  keep <- !vapply(draws, is.null, logical(1))
  structure(list(
    draws = do.call(rbind, draws[keep]),
    losses = if (!is.null(economics)) do.call(rbind, losses[keep]) else NULL,
    params = do.call(rbind, params[keep]),
    weights = if (any(keep)) rep(1 / sum(keep), sum(keep)) else numeric(0),
    n_draws = n_draws, n_rejected = length(rej_draw),
    rejections = data.frame(draw = rej_draw, reason = rej_reason, stringsAsFactors = FALSE),
    sign_changes = data.frame(draw = which(keep), n = sign_changes[keep]),
    method = method,
    specs = attr(sampler, "specs"),
    correlated = attr(sampler, "correlated"),
    economics_specs = economics
  ), class = "cm_mc")
}

draw_economics <- function(economics) {
  pick <- function(x) {
    if (is.null(x)) return(NULL)
    vapply(names(x), function(k) {
      v <- x[[k]]
      if (inherits(v, "cm_dist")) v$r(1) else as.numeric(v)
    }, numeric(1))
  }
  obs <- economics$observed
  uv <- economics$unit_value
  if (is.null(obs) || is.null(uv)) cm_abort("`economics` needs `observed` and `unit_value`.")
  list(observed = pick(as.list(obs)), unit_value = pick(as.list(uv)),
       additional = economics$additional %||% 0)
}

# One-row data frame of a model's numeric inputs, keyed as in cm_sampler().
flatten_model <- function(model) {
  v <- stats::setNames(model$diseases$value, paste0("prob:", model$diseases$id))
  a <- model$associations
  if (!is.null(a) && nrow(a)) {
    v <- c(v, stats::setNames(a$value, paste0("assoc:", a$disease1, ":", a$disease2)))
  }
  i <- model$impacts
  if (!is.null(i) && nrow(i)) {
    div <- if (is.null(i$input_scale)) 1 else ifelse(i$input_scale %in% "percent", 100, 1)
    v <- c(v, stats::setNames(i$value * div, paste0("impact:", i$outcome, ":", i$disease)))
  }
  x <- model$interactions
  if (!is.null(x) && nrow(x)) {
    div <- if (is.null(x$input_scale)) 1 else ifelse(x$input_scale %in% "percent", 100, 1)
    v <- c(v, stats::setNames(x$value * div, paste0("inter:", x$outcome, ":", x$disease1, ":", x$disease2)))
  }
  as.data.frame(as.list(v), check.names = FALSE)
}

#' Reweight Monte Carlo draws (importance sampling)
#'
#' Reweights accepted draws from [cm_monte_carlo()] so that summaries reflect
#' a different input distribution, without re-running the adjustment. Weights
#' are self-normalised, so densities need only be known up to a constant.
#' [cm_scenario()] builds the log-ratio from distributions; use
#' `cm_reweight()` directly for custom targets.
#'
#' The sampling (proposal) distribution must cover the scenario: give the
#' sampler wider spread than the base case (e.g. [dist_mixture()]).
#' Reweighting cannot recover rejected draws: if the scenario puts weight
#' where draws were infeasible, results are conditional on feasibility.
#' Always check the effective sample size.
#'
#' @param mc A `cm_mc` object.
#' @param log_ratio A function taking `mc$params` (one row per accepted draw)
#'   and returning `log(target density / proposal density)` for each row.
#' @return `mc` with updated `weights` and `ess` (effective sample size,
#'   `1 / sum(w^2)`).
#' @export
cm_reweight <- function(mc, log_ratio) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  lw <- log_ratio(mc$params)
  if (length(lw) != nrow(mc$params)) cm_abort("`log_ratio` must return one value per accepted draw.")
  if (anyNA(lw)) cm_abort("`log_ratio` returned missing values.")
  if (all(!is.finite(lw) & lw < 0)) cm_abort("The scenario gives zero weight to every draw.")
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
#' those of `mc$params`, e.g. `"assoc:LAM:SCK"` or `"impact:yield:LAM"`.
#'
#' Inputs drawn with an outcome correlation (copula) cannot be reweighted
#' individually, because changing a marginal also changes the copula
#' coordinates; re-run [cm_monte_carlo()] for such scenarios.
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
  bad <- intersect(names(changes), mc$correlated)
  if (length(bad)) {
    cm_abort(sprintf("Cannot reweight inputs drawn with outcome correlation: %s.", paste(bad, collapse = ", ")))
  }
  cm_reweight(mc, function(p) {
    lr <- rep(0, nrow(p))
    for (k in names(changes)) {
      x <- p[[k]]
      lr <- lr + changes[[k]]$logd(x) - specs[[k]]$logd(x)
    }
    lr
  })
}

#' Summarise Monte Carlo results
#'
#' @param object A `cm_mc` object.
#' @param what `"adjusted"` (adjusted impacts by outcome and disease),
#'   `"loss"` (monetary losses by outcome and disease), `"total"` (total
#'   losses per outcome and overall) or `"rejections"`.
#' @param probs Quantiles to report.
#' @param ... Unused.
#' @return A data frame of (weighted) means, SDs, quantiles and the Monte
#'   Carlo standard error of the mean (based on the effective sample size).
#' @export
summary.cm_mc <- function(object, what = c("adjusted", "loss", "total", "rejections"),
                          probs = c(0.025, 0.5, 0.975), ...) {
  what <- match.arg(what)
  if (what == "rejections") {
    if (!nrow(object$rejections)) return(data.frame(reason = character(0), n = integer(0)))
    tab <- table(object$rejections$reason)
    return(data.frame(reason = names(tab), n = as.integer(tab), stringsAsFactors = FALSE))
  }
  if (is.null(object$draws)) cm_abort("All draws were rejected.")
  w_by_draw <- stats::setNames(object$weights, object$params$draw)
  ess <- 1 / sum(object$weights^2)
  if (what == "adjusted") {
    d <- object$draws
    d$x <- d$adjusted
    keys <- c("outcome", "disease")
  } else {
    if (is.null(object$losses)) cm_abort("No losses: run cm_monte_carlo() with `economics`.")
    d <- object$losses
    d$x <- d$value
    keys <- c("outcome", "disease")
    if (what == "total") {
      by_out <- stats::aggregate(x ~ draw + outcome, data = d, FUN = sum)
      all_out <- stats::aggregate(x ~ draw, data = d, FUN = sum)
      all_out$outcome <- "total"
      d <- rbind(by_out, all_out[, c("draw", "outcome", "x")])
      d$disease <- "all"
    }
  }
  d$w <- w_by_draw[as.character(d$draw)]
  groups <- split(d, list(d$outcome, d$disease), drop = TRUE)
  out <- lapply(groups, function(g) {
    w <- g$w / sum(g$w)
    mu <- sum(w * g$x)
    s <- sqrt(sum(w * (g$x - mu)^2))
    q <- weighted_quantile(g$x, w, probs)
    data.frame(outcome = g$outcome[1], disease = g$disease[1],
               mean = mu, sd = s, mcse = s / sqrt(ess),
               t(stats::setNames(q, paste0("q", probs))),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, unname(out))
  rownames(out) <- NULL
  attr(out, "n_draws") <- object$n_draws
  attr(out, "n_rejected") <- object$n_rejected
  attr(out, "ess") <- ess
  out
}
