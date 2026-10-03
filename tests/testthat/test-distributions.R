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
  expect_equal(dist_beta(2, 6)$mean, 0.25)
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
