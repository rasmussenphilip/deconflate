test_that("quantile and distribution functions are inverses", {
  ds <- list(
    dist_normal(2, 0.5), dist_normal(2.63, 1.44, lower = 0),
    dist_lognormal(log(2), 0.3), dist_beta(13.95, 441.59),
    dist_pert(1.19, 3.30, 10.71), dist_uniform(-1, 3),
    dist_mixture(dist_normal(2, 0.3), dist_normal(3, 1), weights = c(0.8, 0.2))
  )
  u <- c(0.01, 0.2, 0.5, 0.9, 0.99)
  for (d in ds) {
    expect_equal(d$p(d$q(u)), u, tolerance = 1e-6)
    # The density integrates to one.
    lo <- d$q(1e-9); hi <- d$q(1 - 1e-9)
    expect_equal(stats::integrate(function(x) exp(d$logd(x)), lo, hi)$value, 1,
                 tolerance = 1e-4)
  }
})

test_that("means are correct", {
  expect_equal(dist_pert(1, 4, 13)$mean, (1 + 4 * 4 + 13) / 6)
  expect_equal(dist_pert_mean(1, 5, 13)$mean, 5)
  expect_equal(dist_pert_mean(1, 5, 13)$params$mode, 4)
  expect_equal(dist_beta(2, 6)$mean, 0.25)
  expect_equal(dist_beta(2, 6, min = 1, max = 5)$mean, 2)
  expect_equal(dist_uniform(-1, 3)$mean, 1)
  expect_equal(dist_lognormal(0, 0.5)$mean, exp(0.125))
  expect_equal(dist_mixture(dist_fixed(1), dist_fixed(3), weights = c(3, 1))$mean, 1.5)
  tn <- dist_normal(1, 1, lower = 0)
  num <- stats::integrate(function(x) x * exp(tn$logd(x)), 0, Inf)$value
  expect_equal(tn$mean, num, tolerance = 1e-6)
  expect_error(dist_pert_mean(0, 9.9, 10), "No PERT")
})

test_that("truncation and confidence-interval constructors behave", {
  set.seed(3)
  x <- dist_normal(0.5, 1, lower = 0)$r(2000)
  expect_true(all(x >= 0))
  d <- dist_lognormal_ci(2.7, 1.5, 4.9)
  expect_equal(d$q(0.5), 2.7, tolerance = 1e-12)
  expect_equal(exp((log(d$q(0.975)) + log(d$q(0.025))) / 2), 2.7, tolerance = 1e-12)
  expect_equal(log(d$q(0.975)) - log(d$q(0.025)), log(4.9) - log(1.5), tolerance = 1e-10)
  expect_equal(dist_fixed(3)$r(4), rep(3, 4))
})

test_that("distributions carry their support", {
  sup <- function(d) c(d$lower, d$upper)
  expect_equal(sup(dist_fixed(3)), c(3, 3))
  expect_equal(sup(dist_normal(0, 1)), c(-Inf, Inf))
  expect_equal(sup(dist_normal(2.63, 1.44, lower = 0)), c(0, Inf))
  expect_equal(sup(dist_normal(0, 1, lower = -1, upper = 2)), c(-1, 2))
  expect_equal(sup(dist_lognormal(log(2), 0.3)), c(0, Inf))
  expect_equal(sup(dist_lognormal_ci(2.7, 1.5, 4.9)), c(0, Inf))
  expect_equal(sup(dist_beta(2, 3)), c(0, 1))
  expect_equal(sup(dist_beta(2, 3, min = -1, max = 4)), c(-1, 4))
  expect_equal(sup(dist_pert(1.19, 3.30, 10.71)), c(1.19, 10.71))
  expect_equal(sup(dist_pert_mean(1, 5, 13)), c(1, 13))
  expect_equal(sup(dist_uniform(-1, 3)), c(-1, 3))
  # A mixture's support is the union of its components' supports.
  expect_equal(sup(dist_mixture(dist_uniform(0, 1), dist_uniform(2, 3))), c(0, 3))
  expect_equal(sup(dist_mixture(dist_normal(1, 1, lower = 0), dist_uniform(-2, 1))), c(-2, Inf))

  # Only point masses are discrete; a mixture with a point mass is discrete.
  expect_true(dist_fixed(3)$discrete)
  for (d in list(dist_normal(0, 1), dist_lognormal(0, 1), dist_beta(2, 3),
                 dist_pert(0, 1, 2), dist_uniform(0, 1),
                 dist_mixture(dist_normal(0, 1), dist_uniform(0, 1)))) {
    expect_false(d$discrete)
  }
  expect_true(dist_mixture(dist_fixed(1), dist_normal(0, 1))$discrete)

  # Types and classes.
  expect_s3_class(dist_pert(0, 1, 2), "cm_dist")
  expect_equal(dist_pert(0, 1, 2)$type, "pert")
  expect_equal(dist_pert_mean(1, 5, 13)$type, "pert")
  expect_equal(dist_mixture(dist_normal(0, 1), dist_uniform(0, 1))$type, "mixture")
})

test_that("draws, distribution functions and densities respect the support", {
  set.seed(10)
  ds <- list(dist_normal(2.63, 1.44, lower = 0), dist_normal(0, 1, lower = -1, upper = 2),
             dist_lognormal(log(2), 0.3), dist_beta(2, 3, min = -1, max = 4),
             dist_pert(1.19, 3.30, 10.71), dist_uniform(-1, 3),
             dist_mixture(dist_uniform(0, 1), dist_uniform(2, 3)))
  for (d in ds) {
    x <- d$r(300)
    expect_true(all(x >= d$lower & x <= d$upper))
    if (is.finite(d$lower)) {
      expect_equal(d$p(d$lower), 0, tolerance = 1e-12)
      expect_equal(d$logd(d$lower - 1), -Inf)
    }
    if (is.finite(d$upper)) {
      expect_equal(d$p(d$upper), 1, tolerance = 1e-12)
      expect_equal(d$logd(d$upper + 1), -Inf)
    }
  }
  # The mixture has no density in the gap between its components.
  expect_equal(dist_mixture(dist_uniform(0, 1), dist_uniform(2, 3))$logd(1.5), -Inf)
  # A point mass: probability mass 1 at the value, 0 elsewhere.
  f <- dist_fixed(3)
  expect_equal(f$logd(c(3, 2.9)), c(0, -Inf))
  expect_equal(f$p(c(2.9, 3, 3.1)), c(0, 1, 1))
  expect_equal(f$q(c(0.1, 0.9)), c(3, 3))
  expect_equal(f$mean, 3)
})

test_that("constructors validate their arguments", {
  expect_error(dist_fixed(NA_real_), "finite")
  expect_error(dist_normal(0, 0), "positive")
  expect_error(dist_normal(0, 1, lower = 2, upper = 1), "below")
  expect_error(dist_normal(0, 1, lower = 50), "no probability mass")
  expect_error(dist_lognormal(0, 0), "positive")
  expect_error(dist_lognormal_ci(2, 3, 4), "ci_lower")
  expect_error(dist_lognormal_ci(2, 1.5, 3, level = 95), "level")
  expect_error(dist_normal(0, 1, lower = NA), "lower")
  expect_error(dist_beta(0, 1), "positive")
  expect_error(dist_beta(1, 1, min = 1, max = 1), "below")
  expect_error(dist_pert(3, 1, 2), "PERT")
  expect_error(dist_uniform(1, 1), "below")
  expect_error(dist_mixture(1), "cm_dist")
  expect_error(dist_mixture(dist_normal(0, 1), weights = c(0.5, 0.5)), "Invalid mixture weights")
  expect_error(dist_mixture(dist_normal(0, 1), dist_normal(1, 1), weights = c(-1, 2)),
               "Invalid mixture weights")
})

test_that("distributions print", {
  expect_output(print(dist_pert(1.19, 3.30, 10.71)), "pert\\(min = 1.19, mode = 3.3, max = 10.71")
  expect_output(print(dist_normal(2, 0.5)), "normal")
  expect_output(print(dist_mixture(dist_normal(0, 1), dist_uniform(0, 1))),
                "mixture of 2 components")
})

test_that("a model's distributions are keyed by its inputs", {
  m <- example_supplement()
  d <- cm_model(m, m$impacts,
                distributions = list("impact:d1" = dist_normal(2.5, 0.5),
                                     "assoc:d3:d2" = dist_lognormal_ci(3, 2, 4.5),
                                     "prob:d1" = dist_beta(10, 90)))
  # Keys are put in canonical form (associations as listed in their table).
  expect_setequal(names(d$distributions), c("impact:d1", "assoc:d2:d3", "prob:d1"))
  expect_true(all(vapply(d$distributions, inherits, logical(1), "cm_dist")))
  expect_output(print(d), "Uncertain inputs \\(with a distribution\\): 3")
  expect_null(m$distributions)
  bad <- function(dists) cm_model(m, m$impacts, distributions = dists)
  expect_error(bad(list(dist_normal(0, 1))), "named list")
  expect_error(bad(list("impact:d1" = 1)), "named list")
  expect_error(bad(list("impact:d1" = dist_fixed(1), "impact:d1" = dist_fixed(2))), "more than once")
  expect_error(bad(list("assoc:d2:d3" = dist_fixed(3), "assoc:d3:d2" = dist_fixed(3))), "more than once")
  expect_error(bad(list("assoc:d1" = dist_fixed(1))), "assoc:<d1>:<d2>")
  expect_error(bad(list("inter:d1:d2" = dist_fixed(1))), "no interaction")
  expect_error(bad(list("three:d1:d2:d3" = dist_fixed(1))), "no three-way terms")
  expect_error(bad(list("prob:zz" = dist_fixed(0.1))), "unknown disease")
  expect_error(bad(list("foo:d1" = dist_fixed(1))), "unknown input type")
  # A pair without an association has no association to be uncertain about.
  m2 <- m
  m2$associations <- m$associations[-2, ]
  expect_error(cm_model(m2, m2$impacts, distributions = list("assoc:d1:d3" = dist_fixed(1))),
               "no association row")
  # Interactions and three-way terms are keyed in their tables' order.
  mi <- cm_model(cm_population(m$diseases, m$associations, cm_three_way("d3", "d1", "d2", 2)),
                 m$impacts, cm_interactions("d2", "d1", 0.5),
                 distributions = list("inter:d1:d2" = dist_normal(0.5, 0.1),
                                      "three:d1:d2:d3" = dist_lognormal(log(2), 0.3)))
  expect_setequal(names(mi$distributions), c("inter:d2:d1", "three:d3:d1:d2"))
  # Setting an association to a scenario value drops its distribution.
  d2 <- set_association(d, "d2", "d3", 1.5)
  expect_setequal(names(d2$distributions), c("impact:d1", "prob:d1"))
})

test_that("the overall risk can be a distribution within (0, 1), used at its mean", {
  expect_equal(check_overall_risk(NULL, FALSE), list(value = NULL, dist = NULL))
  r <- check_overall_risk(0.25, TRUE)
  expect_equal(r$value, 0.25)
  expect_null(r$dist)
  b <- dist_beta(250, 750)
  rb <- check_overall_risk(b, TRUE)
  expect_equal(rb$value, 0.25)
  expect_identical(rb$dist, b)
  expect_equal(check_overall_risk(dist_uniform(0.2, 0.3), TRUE)$value, 0.25)
  expect_equal(check_overall_risk(dist_pert(0.1, 0.2, 0.6), TRUE)$value, (0.1 + 4 * 0.2 + 0.6) / 6)
  # The support must lie in [0, 1].
  expect_error(check_overall_risk(dist_normal(0.25, 0.1), TRUE), "between 0 and 1")
  expect_error(check_overall_risk(dist_uniform(0.5, 1.5), TRUE), "between 0 and 1")
  expect_error(check_overall_risk(dist_mixture(dist_uniform(0.1, 0.2), dist_uniform(0.9, 1.2)), TRUE),
               "\\[0.9, 1.2\\]")
  expect_silent(check_overall_risk(dist_normal(0.25, 0.1, lower = 0, upper = 1), TRUE))
  expect_error(check_overall_risk(dist_fixed(1), TRUE), "strictly between 0 and 1")
  expect_error(check_overall_risk(0.25, FALSE), class = "deconflate_unsupported")
})
