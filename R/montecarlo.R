#' Monte Carlo propagation of input uncertainty
#'
#' Draws input sets with a user-supplied sampler and adjusts each with
#' [deconflate()]. Because every draw is a complete [cm_model()], disease
#' probabilities and associations are shared across all outcomes within a
#' draw, so the outcomes' uncertainty stays dependent.
#'
#' Draws whose inputs are infeasible (an association incompatible with the
#' sampled marginals, or a jointly infeasible set for the global method) are
#' rejected and counted. Conditioning on feasibility changes the effective
#' input distribution, so the rejection rate is part of the result: report it.
#'
#' @param sampler A function of one argument (the draw index) returning a
#'   [cm_model()].
#' @param n_draws Number of draws.
#' @param method Adjustment method passed to [deconflate()].
#' @param ... Passed to [deconflate()].
#'
#' @return A `cm_mc` object: `draws` (long data frame of raw and adjusted
#'   impacts by draw), `params` (one row per accepted draw with the sampled
#'   inputs, used by [cm_reweight()]), `weights` (initially equal),
#'   `n_draws`, `n_rejected`, `rejections` (draw and reason) and
#'   `sign_changes` (number of adjusted impacts per draw that changed sign).
#' @export
cm_monte_carlo <- function(sampler, n_draws, method = "simultaneous", ...) {
  if (!is.function(sampler)) cm_abort("`sampler` must be a function returning a cm_model().")
  draws <- vector("list", n_draws)
  params <- vector("list", n_draws)
  sign_changes <- integer(n_draws)
  rejected <- data.frame(draw = integer(0), reason = character(0))
  for (d in seq_len(n_draws)) {
    model <- sampler(d)
    res <- tryCatch(
      withCallingHandlers(
        deconflate(model, method = method, warn = FALSE, ...),
        deconflate_nonconvergence = function(w) invokeRestart("muffleWarning")
      ),
      deconflate_infeasible = function(e) e
    )
    if (inherits(res, "condition")) {
      rejected <- rbind(rejected, data.frame(draw = d, reason = conditionMessage(res)))
      next
    }
    draws[[d]] <- cbind(draw = d, res$adjusted[, c("outcome", "disease", "raw", "adjusted")])
    params[[d]] <- cbind(draw = d, flatten_model(model))
    sign_changes[d] <- sum(res$diagnostics$n_sign_changes)
  }
  keep <- !vapply(draws, is.null, logical(1))
  params <- do.call(rbind, params[keep])
  structure(list(
    draws = do.call(rbind, draws[keep]),
    params = params,
    weights = if (any(keep)) rep(1 / sum(keep), sum(keep)) else numeric(0),
    n_draws = n_draws, n_rejected = nrow(rejected), rejections = rejected,
    sign_changes = data.frame(draw = which(keep), n = sign_changes[keep]),
    method = method
  ), class = "cm_mc")
}

# One-row data frame of a model's numeric inputs.
flatten_model <- function(model) {
  v <- stats::setNames(model$diseases$prob, paste0("prob_", model$diseases$id))
  a <- model$associations
  if (!is.null(a) && nrow(a)) {
    v <- c(v, stats::setNames(a$value, paste0("assoc_", a$disease1, "_", a$disease2)))
  }
  i <- model$impacts
  if (!is.null(i) && nrow(i)) {
    v <- c(v, stats::setNames(i$value, paste0("impact_", i$outcome, "_", i$disease)))
  }
  x <- model$interactions
  if (!is.null(x) && nrow(x)) {
    v <- c(v, stats::setNames(x$value, paste0("inter_", x$outcome, "_", x$disease1, "_", x$disease2)))
  }
  as.data.frame(as.list(v), check.names = FALSE)
}

#' Reweight Monte Carlo draws for a scenario (importance sampling)
#'
#' Reweights accepted draws from [cm_monte_carlo()] so that summaries
#' reflect a scenario's input distribution instead of the sampling
#' (proposal) distribution, without re-running the adjustment. Weights are
#' self-normalised, so densities need only be known up to a constant.
#'
#' The proposal must cover the scenario: give the sampler wider spread than
#' the base case (e.g. a defensive mixture). Reweighting cannot recover
#' rejected draws: if the scenario puts weight on regions where draws were
#' infeasible, the result is conditional on feasibility. Check the effective
#' sample size.
#'
#' @param mc A `cm_mc` object.
#' @param log_ratio A function taking `mc$params` (one row per accepted draw)
#'   and returning `log(target density / proposal density)` for each row.
#' @return `mc` with updated `weights` and an `ess` element (effective sample
#'   size, `1 / sum(w^2)`).
#' @export
cm_reweight <- function(mc, log_ratio) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  lw <- log_ratio(mc$params)
  if (length(lw) != nrow(mc$params)) cm_abort("`log_ratio` must return one value per accepted draw.")
  if (anyNA(lw)) cm_abort("`log_ratio` returned missing values.")
  w <- exp(lw - max(lw))
  w <- w / sum(w)
  mc$weights <- w
  mc$ess <- 1 / sum(w^2)
  if (mc$ess < 0.1 * length(w)) {
    cm_warn(sprintf("Effective sample size is low (%.1f of %d draws).", mc$ess, length(w)))
  }
  mc
}

#' Summarise Monte Carlo results
#'
#' @param object A `cm_mc` object.
#' @param probs Quantiles to report.
#' @param ... Unused.
#' @return A data frame of (weighted) means, SDs and quantiles of adjusted
#'   impacts by outcome and disease.
#' @export
summary.cm_mc <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  if (is.null(object$draws)) cm_abort("All draws were rejected.")
  w_by_draw <- stats::setNames(object$weights, object$params$draw)
  d <- object$draws
  d$w <- w_by_draw[as.character(d$draw)]
  groups <- split(d, list(d$outcome, d$disease), drop = TRUE)
  out <- lapply(groups, function(g) {
    w <- g$w / sum(g$w)
    mu <- sum(w * g$adjusted)
    q <- weighted_quantile(g$adjusted, w, probs)
    data.frame(outcome = g$outcome[1], disease = g$disease[1],
               mean = mu, sd = sqrt(sum(w * (g$adjusted - mu)^2)),
               t(stats::setNames(q, paste0("q", probs))),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, unname(out))
  rownames(out) <- NULL
  attr(out, "n_draws") <- object$n_draws
  attr(out, "n_rejected") <- object$n_rejected
  out
}
