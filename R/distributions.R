#' Input distributions for Monte Carlo analysis
#'
#' Constructors for the distributions used by [cm_sampler()]. Each returns a
#' `cm_dist` object with a quantile function (used for sampling and for
#' correlated draws), a distribution function and a log-density (used for
#' importance reweighting in [cm_scenario()]).
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
#' * `dist_mixture()`: a finite mixture, e.g. a defensive mixture of a base
#'   distribution and a wider one, so that scenario reweighting has support.
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
#' @return A `cm_dist` object.
#' @name distributions
#' @examples
#' d <- dist_pert(1.19, 3.30, 10.71)
#' d$mean
#' d$r(5)
#' dist_lognormal_ci(2.7, 1.5, 4.9)
NULL

new_dist <- function(type, params, q, p, logd, mean) {
  structure(list(type = type, params = params, q = q, p = p, logd = logd,
                 mean = mean, r = function(n) q(stats::runif(n))),
            class = "cm_dist")
}

#' @rdname distributions
#' @export
dist_fixed <- function(value) {
  check_numeric(value, "value")
  new_dist("fixed", list(value = value),
           q = function(u) rep(value, length(u)),
           p = function(x) as.numeric(x >= value),
           logd = function(x) ifelse(x == value, 0, -Inf),
           mean = value)
}

#' @rdname distributions
#' @export
dist_normal <- function(mean, sd, lower = -Inf, upper = Inf) {
  check_numeric(c(mean, sd), "mean/sd")
  if (sd <= 0) cm_abort("`sd` must be positive.")
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
           mean = m)
}

#' @rdname distributions
#' @export
dist_lognormal <- function(meanlog, sdlog) {
  if (sdlog <= 0) cm_abort("`sdlog` must be positive.")
  new_dist("lognormal", list(meanlog = meanlog, sdlog = sdlog),
           q = function(u) stats::qlnorm(u, meanlog, sdlog),
           p = function(x) stats::plnorm(x, meanlog, sdlog),
           logd = function(x) stats::dlnorm(x, meanlog, sdlog, log = TRUE),
           mean = exp(meanlog + sdlog^2 / 2))
}

#' @rdname distributions
#' @export
dist_lognormal_ci <- function(estimate, ci_lower, ci_upper, level = 0.95) {
  if (!(ci_lower > 0 && ci_lower < estimate && estimate < ci_upper)) {
    cm_abort("Need 0 < ci_lower < estimate < ci_upper.")
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  dist_lognormal(log(estimate), (log(ci_upper) - log(ci_lower)) / (2 * z))
}

#' @rdname distributions
#' @export
dist_beta <- function(shape1, shape2, min = 0, max = 1) {
  if (shape1 <= 0 || shape2 <= 0) cm_abort("Beta shape parameters must be positive.")
  if (min >= max) cm_abort("`min` must be below `max`.")
  w <- max - min
  new_dist("beta", list(shape1 = shape1, shape2 = shape2, min = min, max = max),
           q = function(u) min + w * stats::qbeta(u, shape1, shape2),
           p = function(x) stats::pbeta((x - min) / w, shape1, shape2),
           logd = function(x) stats::dbeta((x - min) / w, shape1, shape2, log = TRUE) - log(w),
           mean = min + w * shape1 / (shape1 + shape2))
}

#' @rdname distributions
#' @export
dist_pert <- function(min, mode, max, lambda = 4) {
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
  mode <- ((lambda + 2) * mean - min - max) / lambda
  if (mode < min || mode > max) {
    cm_abort(sprintf("No PERT distribution on [%g, %g] has mean %g.", min, max, mean))
  }
  dist_pert(min, mode, max, lambda)
}

#' @rdname distributions
#' @export
dist_uniform <- function(min, max) {
  if (min >= max) cm_abort("`min` must be below `max`.")
  new_dist("uniform", list(min = min, max = max),
           q = function(u) stats::qunif(u, min, max),
           p = function(x) stats::punif(x, min, max),
           logd = function(x) stats::dunif(x, min, max, log = TRUE),
           mean = (min + max) / 2)
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
  pfun <- function(x) {
    out <- 0
    for (j in seq_len(k)) out <- out + w[j] * comps[[j]]$p(x)
    out
  }
  lo <- min(vapply(comps, function(cmp) cmp$q(1e-10), numeric(1)))
  hi <- max(vapply(comps, function(cmp) cmp$q(1 - 1e-10), numeric(1)))
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
    mean = sum(w * vapply(comps, function(cmp) cmp$mean, numeric(1)))
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
