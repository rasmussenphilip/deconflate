#' Input distributions
#'
#' Constructors for the distributions of uncertain inputs, drawn by
#' [deconflate()] (`n_draws`): for the `distributions` of [cm_model()] and
#' for `overall_risk`. (The input tables give the same distributions in
#' their `dist` and `p1`-`p4` columns; see [cm_read_inputs()].) Each returns
#' a `cm_dist` object with a quantile function (used for sampling, including
#' Latin hypercube sampling), a distribution function and a log-density.
#'
#' * `dist_fixed()`: a constant.
#' * `dist_normal()`: normal, optionally truncated to `[lower, upper]`
#'   (e.g. `lower = 0` for odds ratios, which conditions on positive values).
#' * `dist_lognormal()`: log-normal; `dist_lognormal_ci()` builds one from a
#'   point estimate and confidence interval (e.g. a published odds ratio).
#' * `dist_beta()`: beta, optionally rescaled to `[min, max]`.
#' * `dist_pert()`: PERT with `min`, `mode` and `max`. `dist_pert_mean()` builds
#'   one from a mean instead of a mode. Note: in Rasmussen et al. (2024),
#'   Tables 2-4, the central value of PERT distributions is the mode.
#' * `dist_uniform()`: uniform.
#' * `dist_mixture()`: a finite mixture of distributions (e.g. two pooled
#'   sources of evidence).
#'
#' @param value Constant value.
#' @param mean,sd Normal mean and standard deviation.
#' @param lower,upper Truncation bounds.
#' @param meanlog,sdlog Log-normal parameters (log scale).
#' @param estimate,ci_lower,ci_upper,level Point estimate and confidence
#'   interval for `dist_lognormal_ci()`.
#' @param shape1,shape2 Beta shape parameters.
#' @param min,max Range.
#' @param mode Mode of the PERT distribution.
#' @param lambda PERT shape parameter (4 for the standard PERT).
#' @param ... `cm_dist` components for `dist_mixture()`.
#' @param weights Mixture weights (normalised internally).
#'
#' @return A `cm_dist` object, with the quantile function `q`, distribution
#'   function `p`, log-density `logd`, sampler `r`, `mean`, and its support:
#'   `support` (a two-column matrix of disjoint intervals on which the density
#'   is positive), its range `lower` and `upper`, and `discrete = TRUE` for
#'   point masses. The support of a mixture is the union of the supports of
#'   its components with positive weight, so it can have gaps.
#' @name distributions
#' @examples
#' d <- dist_pert(1.19, 3.30, 10.71)
#' d$mean
#' d$r(5)
#' dist_lognormal_ci(2.7, 1.5, 4.9)
NULL

# `lower`/`upper` give the support; `discrete` marks point masses (fixed
# values), whose `logd` is a probability mass rather than a density.
new_dist <- function(type, params, q, p, logd, mean, lower = -Inf, upper = Inf,
                     discrete = FALSE, support = NULL) {
  support <- support %||% matrix(c(lower, upper), ncol = 2)
  structure(list(type = type, params = params, q = q, p = p, logd = logd,
                 mean = mean, r = function(n) q(stats::runif(n)),
                 lower = min(support[, 1]), upper = max(support[, 2]),
                 discrete = discrete, support = support),
            class = "cm_dist")
}

# Support of a distribution as a matrix of disjoint, sorted intervals.
dist_support <- function(d) {
  d$support %||% matrix(c(d$lower %||% -Inf, d$upper %||% Inf), ncol = 2)
}

# Union of intervals (rows of a two-column matrix), merging overlapping or
# touching ones.
interval_union <- function(m) {
  m <- m[order(m[, 1], m[, 2]), , drop = FALSE]
  out <- m[1, , drop = FALSE]
  for (r in seq_len(nrow(m))[-1]) {
    last <- nrow(out)
    if (m[r, 1] <= out[last, 2]) {
      out[last, 2] <- max(out[last, 2], m[r, 2])
    } else {
      out <- rbind(out, m[r, ])
    }
  }
  unname(out)
}

# Does the support of `cover` contain the support of `target`? Each target
# interval must lie inside one interval of the (merged) cover; a gap of
# positive length anywhere in the target's support fails.
support_covers <- function(cover, target) {
  cs <- interval_union(dist_support(cover))
  ts <- interval_union(dist_support(target))
  all(vapply(seq_len(nrow(ts)), function(r) {
    any(cs[, 1] <= ts[r, 1] & cs[, 2] >= ts[r, 2])
  }, logical(1)))
}

format_support <- function(d) {
  s <- interval_union(dist_support(d))
  paste(sprintf("[%g, %g]", s[, 1], s[, 2]), collapse = " u ")
}

#' @rdname distributions
#' @export
dist_fixed <- function(value) {
  check_numeric(value, "value")
  new_dist("fixed", list(value = value),
           q = function(u) rep(value, length(u)),
           p = function(x) as.numeric(x >= value),
           logd = function(x) ifelse(x == value, 0, -Inf),
           mean = value, lower = value, upper = value, discrete = TRUE)
}

#' @rdname distributions
#' @export
dist_normal <- function(mean, sd, lower = -Inf, upper = Inf) {
  check_numeric(c(mean, sd), "mean/sd")
  if (sd <= 0) cm_abort("`sd` must be positive.")
  if (!is.numeric(c(lower, upper)) || anyNA(c(lower, upper))) {
    cm_abort("`lower` and `upper` must be numbers (-Inf and Inf for no bound).")
  }
  if (lower >= upper) cm_abort("`lower` must be below `upper`.")
  Fl <- stats::pnorm(lower, mean, sd)
  Fu <- stats::pnorm(upper, mean, sd)
  Z <- Fu - Fl
  if (Z <= 0) cm_abort("The truncation range has no probability mass.")
  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  m <- mean + sd * (stats::dnorm(a) - stats::dnorm(b)) / Z
  new_dist("normal", list(mean = mean, sd = sd, lower = lower, upper = upper),
           q = function(u) stats::qnorm(Fl + u * Z, mean, sd),
           p = function(x) pmin(pmax((stats::pnorm(x, mean, sd) - Fl) / Z, 0), 1),
           logd = function(x) ifelse(x < lower | x > upper, -Inf,
                                     stats::dnorm(x, mean, sd, log = TRUE) - log(Z)),
           mean = m, lower = lower, upper = upper)
}

#' @rdname distributions
#' @export
dist_lognormal <- function(meanlog, sdlog) {
  check_numeric(c(meanlog, sdlog), "meanlog/sdlog")
  if (sdlog <= 0) cm_abort("`sdlog` must be positive.")
  new_dist("lognormal", list(meanlog = meanlog, sdlog = sdlog),
           q = function(u) stats::qlnorm(u, meanlog, sdlog),
           p = function(x) stats::plnorm(x, meanlog, sdlog),
           logd = function(x) stats::dlnorm(x, meanlog, sdlog, log = TRUE),
           mean = exp(meanlog + sdlog^2 / 2), lower = 0, upper = Inf)
}

#' @rdname distributions
#' @export
dist_lognormal_ci <- function(estimate, ci_lower, ci_upper, level = 0.95) {
  check_numeric(c(estimate, ci_lower, ci_upper, level), "estimate/ci_lower/ci_upper/level")
  if (!(ci_lower > 0 && ci_lower < estimate && estimate < ci_upper)) {
    cm_abort("Need 0 < ci_lower < estimate < ci_upper.")
  }
  if (!(level > 0 && level < 1)) cm_abort("`level` must be between 0 and 1 (e.g. 0.95).")
  z <- stats::qnorm(1 - (1 - level) / 2)
  dist_lognormal(log(estimate), (log(ci_upper) - log(ci_lower)) / (2 * z))
}

#' @rdname distributions
#' @export
dist_beta <- function(shape1, shape2, min = 0, max = 1) {
  check_numeric(c(shape1, shape2, min, max), "shape1/shape2/min/max")
  if (shape1 <= 0 || shape2 <= 0) cm_abort("Beta shape parameters must be positive.")
  if (min >= max) cm_abort("`min` must be below `max`.")
  w <- max - min
  new_dist("beta", list(shape1 = shape1, shape2 = shape2, min = min, max = max),
           q = function(u) min + w * stats::qbeta(u, shape1, shape2),
           p = function(x) stats::pbeta((x - min) / w, shape1, shape2),
           logd = function(x) stats::dbeta((x - min) / w, shape1, shape2, log = TRUE) - log(w),
           mean = min + w * shape1 / (shape1 + shape2), lower = min, upper = max)
}

#' @rdname distributions
#' @export
dist_pert <- function(min, mode, max, lambda = 4) {
  check_numeric(c(min, mode, max, lambda), "min/mode/max/lambda")
  if (!(min < max && mode >= min && mode <= max)) {
    cm_abort("PERT needs min <= mode <= max and min < max.")
  }
  a <- 1 + lambda * (mode - min) / (max - min)
  b <- 1 + lambda * (max - mode) / (max - min)
  d <- dist_beta(a, b, min, max)
  d$type <- "pert"
  d$params <- list(min = min, mode = mode, max = max, lambda = lambda)
  d
}

#' @rdname distributions
#' @export
dist_pert_mean <- function(min, mean, max, lambda = 4) {
  check_numeric(c(min, mean, max, lambda), "min/mean/max/lambda")
  mode <- ((lambda + 2) * mean - min - max) / lambda
  if (mode < min || mode > max) {
    cm_abort(sprintf("No PERT distribution on [%g, %g] has mean %g.", min, max, mean))
  }
  dist_pert(min, mode, max, lambda)
}

#' @rdname distributions
#' @export
dist_uniform <- function(min, max) {
  check_numeric(c(min, max), "min/max")
  if (min >= max) cm_abort("`min` must be below `max`.")
  new_dist("uniform", list(min = min, max = max),
           q = function(u) stats::qunif(u, min, max),
           p = function(x) stats::punif(x, min, max),
           logd = function(x) stats::dunif(x, min, max, log = TRUE),
           mean = (min + max) / 2, lower = min, upper = max)
}

#' @rdname distributions
#' @export
dist_mixture <- function(..., weights = NULL) {
  comps <- list(...)
  if (!length(comps) || !all(vapply(comps, inherits, logical(1), "cm_dist"))) {
    cm_abort("Mixture components must be cm_dist objects.")
  }
  k <- length(comps)
  w <- weights %||% rep(1 / k, k)
  if (length(w) != k || any(w < 0) || sum(w) <= 0) cm_abort("Invalid mixture weights.")
  w <- w / sum(w)
  pos <- which(w > 0)
  pfun <- function(x) {
    out <- 0
    for (j in seq_len(k)) out <- out + w[j] * comps[[j]]$p(x)
    out
  }
  lo <- min(vapply(comps[pos], function(cmp) cmp$q(1e-10), numeric(1)))
  hi <- max(vapply(comps[pos], function(cmp) cmp$q(1 - 1e-10), numeric(1)))
  d <- new_dist(
    "mixture", list(components = comps, weights = w),
    q = function(u) vapply(u, function(ui) {
      stats::uniroot(function(x) pfun(x) - ui, c(lo, hi), tol = 1e-10,
                     extendInt = "yes")$root
    }, numeric(1)),
    p = pfun,
    logd = function(x) {
      dens <- 0
      for (j in seq_len(k)) dens <- dens + w[j] * exp(comps[[j]]$logd(x))
      log(dens)
    },
    mean = sum(w * vapply(comps, function(cmp) cmp$mean, numeric(1))),
    # Only components with positive weight give support.
    support = interval_union(do.call(rbind, lapply(comps[pos], dist_support))),
    discrete = any(vapply(comps[pos], function(cmp) isTRUE(cmp$discrete), logical(1)))
  )
  # Direct sampling is faster than inverting the mixture distribution function.
  d$r <- function(n) {
    j <- sample.int(k, n, replace = TRUE, prob = w)
    vapply(j, function(jj) comps[[jj]]$r(1), numeric(1))
  }
  d
}

#' @export
print.cm_dist <- function(x, ...) {
  par <- x$params
  if (x$type == "mixture") {
    cat(sprintf("<cm_dist> mixture of %d components (weights %s)\n",
                length(par$components), paste(round(par$weights, 3), collapse = ", ")))
  } else {
    cat(sprintf("<cm_dist> %s(%s), mean %.4g\n", x$type,
                paste(names(par), signif(unlist(par), 4), sep = " = ", collapse = ", "),
                x$mean))
  }
  invisible(x)
}
