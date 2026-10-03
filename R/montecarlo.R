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
#' @section Several methods:
#' With more than one `method`, every accepted draw is adjusted with each
#' method, so the methods are compared on identical inputs. A draw is
#' rejected if any method fails on it. Use [compare_methods()] on the result
#' for a side-by-side table.
#'
#' @section Stabilising the estimates:
#' [summary.cm_mc()] flags unstable estimates and suggests remedies (see
#' [cm_diagnose()]). Two variance-reduction options are available here:
#' * `sampling = "lhs"`: Latin hypercube sampling of the inputs that are not
#'   drawn with an outcome correlation. It stratifies each input's
#'   distribution and usually reduces the Monte Carlo error of means.
#' * `proposal`: importance sampling. Named inputs are drawn from the given
#'   proposal distributions instead of their own, and each draw is weighted
#'   by the ratio of the input densities. A proposal must cover the whole
#'   range of the input's own distribution; a defensive mixture such as the
#'   one built by [cm_suggest_proposal()] does. Importance sampling reduces
#'   the error of means that exist; it cannot fix a mean that does not exist
#'   (see [cm_diagnose()]).
#'
#' @param sampler A function of the draw index returning a [cm_model()],
#'   typically from [cm_sampler()].
#' @param n_draws Number of draws.
#' @param method Adjustment method(s) passed to [deconflate()]: one or more
#'   of `"simultaneous"`, `"published"` and `"global"`.
#' @param economics Optional list with named elements `observed` and
#'   `unit_value` (and optionally `additional`), as in [value_losses()].
#'   Elements may be numbers or `cm_dist` objects (drawn per draw). When
#'   given, productivity gaps and monetary losses are computed for every
#'   draw. Hazard-ratio outcomes cannot be valued here; use
#'   [attributable_risk()] on individual results.
#' @param seed Optional random seed, for reproducibility.
#' @param progress Logical: print progress every 10% of draws?
#' @param sampling `"random"` (default) or `"lhs"` (Latin hypercube; needs a
#'   [cm_sampler()]).
#' @param proposal Optional named list of `cm_dist` objects (importance
#'   sampling proposals), keyed as in `params` (e.g.
#'   `"impact:fertility:SCK"`). Needs a [cm_sampler()].
#' @param ... Passed to [deconflate()].
#'
#' @return A `cm_mc` object with elements:
#'   * `draws`: long data frame of raw and adjusted impacts by draw and
#'     method;
#'   * `losses`: long data frame of gaps and values by draw, method, outcome
#'     and disease (if `economics` was given);
#'   * `params`: one row per accepted draw with the sampled inputs, keyed
#'     `prob:<id>`, `assoc:<d1>:<d2>`, `impact:<outcome>:<disease>`,
#'     `inter:<outcome>:<d1>:<d2>`, `observed:<outcome>`, `unit_value:<outcome>`;
#'   * `weights` (equal, or importance weights), `log_weights`, `ess`,
#'     `n_draws`, `n_rejected`, `rejections` (draw and reason),
#'     `sign_changes`, `method`, `sampling` and `proposal`.
#' @export
#' @examples
#' s <- cm_sampler(example_supplement(),
#'                 impacts = list("yield:d1" = dist_normal(2.5, 0.5)))
#' mc <- cm_monte_carlo(s, 100, method = c("published", "simultaneous"), seed = 1)
#' compare_methods(mc)
cm_monte_carlo <- function(sampler, n_draws, method = "simultaneous",
                           economics = NULL, seed = NULL, progress = FALSE,
                           sampling = c("random", "lhs"), proposal = NULL, ...) {
  if (!is.function(sampler)) cm_abort("`sampler` must be a function returning a cm_model().")
  method <- unique(match.arg(method, c("simultaneous", "published", "global"), several.ok = TRUE))
  sampling <- match.arg(sampling)
  specs <- attr(sampler, "specs")
  correlated <- attr(sampler, "correlated")
  if ((sampling == "lhs" || !is.null(proposal)) && is.null(specs)) {
    cm_abort("Latin hypercube sampling and importance sampling need a sampler from cm_sampler().")
  }
  if (!is.null(proposal)) {
    if (is.null(names(proposal)) || !all(vapply(proposal, inherits, logical(1), "cm_dist"))) {
      cm_abort("`proposal` must be a named list of cm_dist objects.")
    }
    unk <- setdiff(names(proposal), names(specs))
    if (length(unk)) cm_abort(sprintf("No sampled input named: %s.", paste(unk, collapse = ", ")))
    bad <- intersect(names(proposal), correlated)
    if (length(bad)) {
      cm_abort(sprintf("Inputs drawn with outcome correlation cannot have a proposal: %s.",
                       paste(bad, collapse = ", ")))
    }
  }
  if (!is.null(seed)) set.seed(seed)
  U <- NULL
  if (sampling == "lhs") {
    keys <- setdiff(names(specs), correlated)
    if (length(keys)) {
      U <- vapply(keys, function(k) (sample.int(n_draws) - stats::runif(n_draws)) / n_draws,
                  numeric(n_draws))
      U <- matrix(U, nrow = n_draws, dimnames = list(NULL, keys))
    }
  }
  draws <- vector("list", n_draws)
  losses <- vector("list", n_draws)
  params <- vector("list", n_draws)
  logw <- rep(0, n_draws)
  sign_changes <- integer(n_draws)
  rej_draw <- integer(0)
  rej_reason <- character(0)
  step <- max(1L, n_draws %/% 10L)

  for (d in seq_len(n_draws)) {
    if (progress && d %% step == 0L) message(sprintf("Draw %d of %d", d, n_draws))
    vals <- NULL
    if (!is.null(proposal)) {
      vals <- vapply(names(proposal), function(k) {
        if (!is.null(U) && k %in% colnames(U)) proposal[[k]]$q(U[d, k]) else proposal[[k]]$r(1)
      }, numeric(1))
      lw <- sum(vapply(names(proposal), function(k) {
        specs[[k]]$logd(vals[[k]]) - proposal[[k]]$logd(vals[[k]])
      }, numeric(1)))
      if (is.nan(lw) || lw == Inf) {
        cm_abort("A proposal has zero density where the input's own distribution does not; use a defensive mixture (see cm_suggest_proposal()).")
      }
      logw[d] <- lw
    }
    u <- NULL
    if (!is.null(U)) {
      cols <- setdiff(colnames(U), names(proposal))
      if (length(cols)) u <- stats::setNames(as.numeric(U[d, cols]), cols)
    }
    res <- tryCatch(
      withCallingHandlers({
        model <- if (is.null(u) && is.null(vals)) sampler(d) else sampler(d, u = u, values = vals)
        rs <- lapply(method, function(m) deconflate(model, method = m, warn = FALSE, ...))
        eco <- if (!is.null(economics)) draw_economics(economics) else NULL
        lv <- if (!is.null(eco)) {
          lapply(rs, function(r) value_losses(productivity_gap(r, eco$observed),
                                              eco$unit_value, eco$additional))
        } else NULL
        list(model = model, results = rs, eco = eco, lv = lv)
      }, deconflate_nonconvergence = function(w) invokeRestart("muffleWarning")),
      deconflate_infeasible = function(e) e
    )
    if (inherits(res, "condition")) {
      rej_draw <- c(rej_draw, d)
      rej_reason <- c(rej_reason, conditionMessage(res))
      next
    }
    draws[[d]] <- do.call(rbind, lapply(seq_along(method), function(j) {
      a <- res$results[[j]]$adjusted
      cbind(draw = d, method = method[j], a[, c("outcome", "disease", "raw", "adjusted", "scale")],
            stringsAsFactors = FALSE)
    }))
    p <- flatten_model(res$model)
    if (!is.null(res$eco)) {
      p <- cbind(p, as.data.frame(as.list(c(
        stats::setNames(res$eco$observed, paste0("observed:", names(res$eco$observed))),
        stats::setNames(res$eco$unit_value, paste0("unit_value:", names(res$eco$unit_value))))),
        check.names = FALSE))
      losses[[d]] <- do.call(rbind, lapply(seq_along(method), function(j) {
        cbind(draw = d, method = method[j], res$lv[[j]]$by_disease, stringsAsFactors = FALSE)
      }))
    }
    params[[d]] <- cbind(draw = d, p)
    sign_changes[d] <- sum(vapply(res$results, function(r) sum(r$diagnostics$n_sign_changes),
                                  numeric(1)))
  }
  keep <- !vapply(draws, is.null, logical(1))
  lw <- logw[keep]
  w <- if (any(keep)) exp(lw - max(lw)) else numeric(0)
  w <- if (length(w)) w / sum(w) else w
  structure(list(
    draws = do.call(rbind, draws[keep]),
    losses = if (!is.null(economics)) do.call(rbind, losses[keep]) else NULL,
    params = do.call(rbind, params[keep]),
    weights = w,
    log_weights = lw,
    ess = if (length(w)) 1 / sum(w^2) else 0,
    n_draws = n_draws, n_rejected = length(rej_draw),
    rejections = data.frame(draw = rej_draw, reason = rej_reason, stringsAsFactors = FALSE),
    sign_changes = data.frame(draw = which(keep), n = sign_changes[keep]),
    method = method,
    sampling = sampling,
    proposal = proposal,
    specs = specs,
    correlated = correlated,
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
#'   `1 / sum(w^2)`). Importance weights from the original run (see the
#'   `proposal` argument of [cm_monte_carlo()]) are kept.
#' @export
cm_reweight <- function(mc, log_ratio) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  lw <- log_ratio(mc$params)
  if (length(lw) != nrow(mc$params)) cm_abort("`log_ratio` must return one value per accepted draw.")
  if (anyNA(lw)) cm_abort("`log_ratio` returned missing values.")
  # Keep any importance-sampling weights from the original run.
  if (length(mc$log_weights) == length(lw)) lw <- lw + mc$log_weights
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
#' Summarises adjusted impacts (or losses) over the accepted draws, and checks
#' each estimate's stability. When some estimates look unstable, a message
#' explains why and suggests a remedy (see [cm_diagnose()]).
#'
#' @param object A `cm_mc` object.
#' @param what `"adjusted"` (adjusted impacts by outcome and disease),
#'   `"loss"` (monetary losses by outcome and disease), `"total"` (total
#'   losses per outcome and overall) or `"rejections"`.
#' @param probs Quantiles to report.
#' @param trim Fraction trimmed from each tail for `trimmed_mean`.
#' @param diagnose Logical: print suggestions when estimates look unstable?
#' @param ... Unused.
#' @return A data frame of (weighted) means, SDs, the Monte Carlo standard
#'   error of the mean (`mcse`, based on the effective sample size),
#'   quantiles, a trimmed mean and stability diagnostics:
#'   * `rel_mcse`: `mcse` relative to the absolute mean;
#'   * `tail_share`: the share of the variance contributed by the most
#'     extreme 1% of draws;
#'   * `stability`: `"ok"`, `"imprecise"` (`rel_mcse` above 5%),
#'     `"heavy_tail"` (`tail_share` above 60%) or `"no_mean"` (the published
#'     approximation divides by a quantity that changes sign across draws,
#'     so the mean does not exist).
#' @export
summary.cm_mc <- function(object, what = c("adjusted", "loss", "total", "rejections"),
                          probs = c(0.025, 0.5, 0.975), trim = 0.05, diagnose = TRUE, ...) {
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
  } else {
    if (is.null(object$losses)) cm_abort("No losses: run cm_monte_carlo() with `economics`.")
    d <- object$losses
    d$x <- d$value
    d$raw <- NA_real_
    d$adjusted <- NA_real_
    if (is.null(d$method)) d$method <- object$method[1]
    if (what == "total") {
      by_out <- stats::aggregate(x ~ draw + outcome + method, data = d, FUN = sum)
      all_out <- stats::aggregate(x ~ draw + method, data = d, FUN = sum)
      all_out$outcome <- "total"
      d <- rbind(by_out, all_out[, c("draw", "outcome", "method", "x")])
      d$disease <- "all"
      d$raw <- NA_real_
      d$adjusted <- NA_real_
    }
  }
  if (is.null(d$method)) d$method <- object$method[1]
  if (is.null(d$scale)) d$scale <- NA_character_
  d$w <- w_by_draw[as.character(d$draw)]
  grp <- list(factor(d$outcome, levels = unique(d$outcome)),
              factor(d$disease, levels = unique(d$disease)),
              factor(d$method, levels = unique(d$method)))
  groups <- split(d, grp, drop = TRUE)
  stats_rows <- lapply(groups, function(g) mc_stats(g, probs, trim, ess))
  out <- do.call(rbind, lapply(stats_rows, `[[`, "row"))
  info <- do.call(rbind, lapply(stats_rows, `[[`, "info"))
  rownames(out) <- NULL
  rownames(info) <- NULL
  attr(out, "n_draws") <- object$n_draws
  attr(out, "n_rejected") <- object$n_rejected
  attr(out, "ess") <- ess
  diag <- diagnosis_table(out, info, ess)
  attr(out, "diagnosis") <- diag
  if (diagnose && nrow(diag)) message(format_diagnosis(diag))
  out
}

# Weighted summary and stability statistics for one group of draws.
mc_stats <- function(g, probs, trim, ess) {
  n_bad <- sum(!is.finite(g$x))
  g <- g[is.finite(g$x), , drop = FALSE]
  if (!nrow(g)) cm_abort("All draws of an estimate are non-finite.")
  w <- g$w / sum(g$w)
  x <- g$x
  mu <- sum(w * x)
  s <- sqrt(sum(w * (x - mu)^2))
  q <- weighted_quantile(x, w, probs)
  lim <- weighted_quantile(x, w, c(trim, 1 - trim))
  inside <- x >= lim[1] & x <= lim[2]
  tmean <- if (any(inside)) sum(w[inside] * x[inside]) / sum(w[inside]) else NA_real_
  mcse <- s / sqrt(ess)
  rel <- if (isTRUE(abs(mu) > 1e-12)) mcse / abs(mu) else NA_real_
  cdev <- w * (x - mu)^2
  k <- max(1L, ceiling(0.01 * length(x)))
  tail_share <- if (isTRUE(sum(cdev) > 0)) sum(sort(cdev, decreasing = TRUE)[seq_len(k)]) / sum(cdev) else 0
  # The published approximation m^2 / (m + c) has a pole where m + c = 0. If
  # m + c (recovered as m^2 / adjusted) takes both signs across draws, the
  # pole lies inside the input distribution and the mean does not exist.
  flip <- NA_real_
  if (identical(g$method[1], "published") && all(!is.na(g$raw))) {
    hr <- identical(g$scale[1], "hazard_ratio")
    m <- if (hr) g$raw - 1 else g$raw
    a <- if (hr) g$adjusted - 1 else g$adjusted
    ok <- is.finite(m) & is.finite(a) & abs(m) > 1e-12 & abs(a) > 1e-12
    if (sum(ok) > 1) {
      den <- sign(m[ok]^2 / a[ok])
      flip <- min(mean(den > 0), mean(den < 0))
    }
  }
  stability <- if ((!is.na(flip) && flip > 0) || n_bad > 0) {
    "no_mean"
  } else if (length(x) >= 50 && tail_share > 0.6) {
    "heavy_tail"
  } else if (!is.na(rel) && rel > 0.05) {
    "imprecise"
  } else {
    "ok"
  }
  row <- data.frame(outcome = g$outcome[1], disease = g$disease[1], method = g$method[1],
                    mean = mu, sd = s, mcse = mcse,
                    t(stats::setNames(q, paste0("q", probs))),
                    trimmed_mean = tmean, rel_mcse = rel, tail_share = tail_share,
                    stability = stability, check.names = FALSE, stringsAsFactors = FALSE)
  list(row = row, info = data.frame(flip = flip, n_bad = n_bad, stringsAsFactors = FALSE))
}

diagnosis_table <- function(out, info, ess) {
  bad <- out$stability != "ok"
  if (!any(bad)) {
    return(data.frame(outcome = character(0), disease = character(0), method = character(0),
                      stability = character(0), detail = character(0),
                      suggestion = character(0), stringsAsFactors = FALSE))
  }
  o <- out[bad, , drop = FALSE]
  f <- info$flip[bad]
  detail <- character(nrow(o))
  sugg <- character(nrow(o))
  for (r in seq_len(nrow(o))) {
    st <- o$stability[r]
    if (st == "no_mean" && info$n_bad[bad][r] > 0 && (is.na(f[r]) || f[r] == 0)) {
      detail[r] <- sprintf("%d draw(s) gave non-finite values (division by zero), so the mean does not exist", info$n_bad[bad][r])
      sugg[r] <- "report the median or trimmed_mean, or use method = \"simultaneous\" (exact, no division)."
    } else if (st == "no_mean") {
      detail[r] <- sprintf("the published approximation divides by m + c, which changes sign in %.1f%% of draws, so the mean does not exist", 100 * f[r])
      sugg[r] <- "report the median or trimmed_mean, or use method = \"simultaneous\" (exact, no division). More draws, Latin hypercube or importance sampling will not make this mean converge."
    } else if (st == "heavy_tail") {
      detail[r] <- sprintf("the most extreme 1%% of draws contribute %.0f%% of the variance", 100 * o$tail_share[r])
      sugg[r] <- sprintf("importance sampling: prop <- cm_suggest_proposal(mc, \"%s\", \"%s\", method = \"%s\"), then cm_monte_carlo(sampler, n_draws, method = \"%s\", proposal = prop). Also report the median.",
                         o$outcome[r], o$disease[r], o$method[r], o$method[r])
    } else {
      need <- ceiling(ess * (o$rel_mcse[r] / 0.02)^2)
      detail[r] <- sprintf("the Monte Carlo standard error is %.1f%% of the mean", 100 * o$rel_mcse[r])
      sugg[r] <- sprintf("increase n_draws to about %s for a 2%% standard error, or use sampling = \"lhs\" (Latin hypercube).",
                         format(need, big.mark = ",", scientific = FALSE))
    }
  }
  data.frame(outcome = o$outcome, disease = o$disease, method = o$method,
             stability = o$stability, detail = detail, suggestion = sugg,
             stringsAsFactors = FALSE)
}

format_diagnosis <- function(diag) {
  lines <- sprintf("* %s / %s (%s): %s.\n    Suggestion: %s", diag$outcome, diag$disease,
                   diag$method, diag$detail, diag$suggestion)
  paste0(sprintf("%d Monte Carlo estimate(s) may be unstable:\n", nrow(diag)),
         paste(lines, collapse = "\n"),
         "\n(See ?cm_diagnose; use summary(..., diagnose = FALSE) to silence this message.)")
}

#' Diagnose unstable Monte Carlo estimates
#'
#' Lists the estimates that [summary.cm_mc()] flags as unstable, with the
#' reason and a suggested remedy:
#' * `"no_mean"`: with the published approximation, `m^2 / (m + c)` has a
#'   pole where `m + c = 0`. If `m + c` takes both signs across draws, the
#'   adjusted impact's distribution has tails so heavy that its mean does
#'   not exist, and no sampling scheme makes it converge. Report the median
#'   or a trimmed mean, or use the exact (`"simultaneous"`) method. (Strictly,
#'   with unbounded input distributions such as the normal, the pole is
#'   always inside the support; the flag is raised when draws actually reach
#'   it.)
#' * `"heavy_tail"`: a few extreme draws dominate the variance. If they come
#'   from an identifiable region of one input, importance sampling
#'   ([cm_suggest_proposal()] and the `proposal` argument of
#'   [cm_monte_carlo()]) samples that region more often and down-weights it,
#'   which reduces the Monte Carlo error.
#' * `"imprecise"`: the Monte Carlo standard error is large relative to the
#'   mean. Use more draws, or Latin hypercube sampling (`sampling = "lhs"`).
#'
#' @param mc A `cm_mc` object.
#' @param what Passed to [summary.cm_mc()].
#' @return A data frame (class `cm_diagnosis`) with one row per flagged
#'   estimate: outcome, disease, method, stability, detail and suggestion.
#' @export
cm_diagnose <- function(mc, what = c("adjusted", "loss", "total")) {
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
#' one adjusted impact, and builds a defensive mixture proposal for it: half
#' the input's own distribution and half a uniform distribution over the
#' range of that input in the extreme draws. Re-running [cm_monte_carlo()]
#' with this proposal samples the region that produces extreme values more
#' often and down-weights it, which reduces the Monte Carlo error when the
#' mean exists. Weights are bounded by `1 / (1 - weight)`, so the run cannot
#' be dominated by a single draw.
#'
#' It cannot help when the mean does not exist (stability `"no_mean"`; see
#' [cm_diagnose()]).
#'
#' @param mc A `cm_mc` object from a [cm_sampler()]-based run.
#' @param outcome,disease The estimate to stabilise.
#' @param method Adjustment method (default: the first in `mc`).
#' @param top Fraction of draws treated as extreme.
#' @param weight Weight of the uniform component in the mixture.
#' @return A named list with one `cm_dist`, for the `proposal` argument of
#'   [cm_monte_carlo()]. The attribute `"explanation"` describes it.
#' @export
cm_suggest_proposal <- function(mc, outcome, disease, method = NULL, top = 0.02, weight = 0.5) {
  if (!inherits(mc, "cm_mc")) cm_abort("`mc` must come from cm_monte_carlo().")
  if (is.null(mc$specs)) cm_abort("`mc` has no recorded distributions; use cm_sampler().")
  if (!(weight > 0 && weight < 1)) cm_abort("`weight` must be between 0 and 1.")
  method <- method %||% mc$method[1]
  d <- mc$draws
  if (is.null(d$method)) d$method <- mc$method[1]
  d <- d[d$outcome == outcome & d$disease == disease & d$method == method, , drop = FALSE]
  if (nrow(d) < 20L) cm_abort("Too few accepted draws for this estimate.")
  keys <- setdiff(names(mc$specs), mc$correlated)
  keys <- keys[vapply(keys, function(k) mc$specs[[k]]$type != "fixed", logical(1))]
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
  top_v <- p[[key]][p$draw %in% top_draws]
  all_v <- p[[key]]
  span <- diff(range(all_v))
  if (span <= 0) span <- abs(mean(all_v)) + 1
  spec <- mc$specs[[key]]
  bound <- function(u, default) {
    b <- tryCatch(spec$q(u), error = function(e) default)
    if (is.finite(b)) b else default
  }
  lo <- max(min(top_v) - 0.1 * span, bound(0, -Inf))
  hi <- min(max(top_v) + 0.1 * span, bound(1, Inf))
  if (!(hi > lo)) cm_abort("Could not find a range for the proposal.")
  prop <- dist_mixture(spec, dist_uniform(lo, hi), weights = c(1 - weight, weight))
  out <- stats::setNames(list(prop), key)
  attr(out, "explanation") <- sprintf(
    "Proposal for %s: %.0f%% its own distribution, %.0f%% uniform on [%.4g, %.4g], the range of %s in the %d most extreme draws of %s / %s (%s).",
    key, 100 * (1 - weight), 100 * weight, lo, hi, key, length(top_draws), outcome, disease, method)
  message(attr(out, "explanation"))
  out
}
