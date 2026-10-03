# Threshold searches (cm_threshold). Reference values:
# inst/validation/reference_v030.py.

test_that("a relative change in the aggregate is found for an association", {
  m <- example_supplement()
  # d1:d3 is independent in the example (baseline OR 1).
  th <- cm_threshold(m, "assoc:d1:d3", c(1, 50), conclusion = "change", target = -0.1)
  expect_s3_class(th, "cm_threshold")
  expect_equal(th$baseline, 1)
  expect_true(th$log_scale)
  expect_equal(th$method, "simultaneous")
  t1 <- th$thresholds[th$thresholds$status == "threshold", ]
  expect_equal(nrow(t1), 1)
  expect_equal(t1$threshold, 5.263907242734966, tolerance = 1e-6)
  expect_lt(t1$lower, t1$threshold)
  expect_gt(t1$upper, t1$threshold)
  expect_equal(th$summary$result, "threshold(s) found")
  expect_equal(nrow(th$regions), 0)
  # Results just below and above the threshold are returned.
  expect_equal(length(th$details), 1)
  expect_true(all(c("below", "above") %in% names(th$details[[1]])))
  expect_output(print(th), "Crossings")
  th2 <- cm_threshold(m, "assoc:d2:d3", c(0.2, 20), conclusion = "change", target = 0.1)
  expect_equal(th2$thresholds$threshold[th2$thresholds$status == "threshold"],
               1.4333305414665463, tolerance = 1e-6)
})

test_that("a model-implied sign change of an adjusted impact is found", {
  th <- cm_threshold(example_supplement(), "impact:d1", c(-5, 5), conclusion = "sign",
                     diseases = "d1")
  expect_false(th$log_scale)
  expect_equal(th$thresholds$item, "d1")
  expect_equal(th$thresholds$status, "threshold")
  expect_equal(th$thresholds$threshold, 0.37414199213982113, tolerance = 1e-6)
  expect_match(th$thresholds$description, "model-implied")
  expect_match(th$thresholds$description, "becomes positive")
})

test_that("varying a probability or an impact draws no random numbers", {
  set.seed(1)
  r1 <- stats::runif(1)
  set.seed(1)
  cm_threshold(example_supplement(), "impact:d1", c(-5, 5), conclusion = "sign", diseases = "d1",
               n_grid = 5)
  cm_threshold(example_supplement(), "prob:d1", c(0.05, 0.5), conclusion = "rank", n_grid = 5)
  expect_equal(stats::runif(1), r1)
})

test_that("a pole of the published approximation is a discontinuity, not a threshold", {
  th <- cm_threshold(example_supplement(), "impact:d1", c(-3, 3), conclusion = "sign",
                     diseases = "d1", method = "published")
  st <- th$thresholds
  # Exactly one crossing: the pole. Touching zero at m1 = 0 (the grid point
  # 0 is exact; the published value is 0 there and positive on both sides)
  # is not a crossing.
  expect_equal(nrow(st), 1)
  expect_equal(st$status, "discontinuity")
  expect_lt(st$lower, -0.5266329181546854)
  expect_gt(st$upper, -0.5266329181546854)
  expect_lt(st$upper - st$lower, 1e-6)
  expect_true(is.na(st$threshold))
  expect_lt(st$below, 0)
  expect_gt(st$above, 0)
  expect_equal(th$summary$n_thresholds, 0L)
  expect_equal(th$summary$result, "no threshold; see discontinuities or unresolved crossings")
  expect_equal(th$scan$status[th$scan$value == 0], "ok")
  expect_equal(unname(th$values[th$scan$value == 0, "d1"]), 0)
})

test_that("exact zeros on the grid count only when the sign differs on both sides", {
  f <- c(1, 0, 1, 0, -1, 0, 0, 2, NA, -1, 0)
  ok <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  # 1 -> 0 -> 1 touches zero; 1 -> 0 -> -1 and -1 -> 0 -> 0 -> 2 cross; the
  # sign change across the unusable point 9 and the final zero are not
  # crossings.
  expect_equal(threshold_brackets(f, ok), rbind(c(3L, 5L), c(5L, 8L)))
  expect_equal(threshold_brackets(c(2, -1), c(TRUE, TRUE)), rbind(c(1L, 2L)))
  expect_equal(nrow(threshold_brackets(c(1, 0, NA, -1), c(TRUE, TRUE, FALSE, TRUE))), 0L)
  expect_equal(nrow(threshold_brackets(c(0, 0, 0), rep(TRUE, 3))), 0L)
})

test_that("the aggregate crossing a target is found", {
  th <- cm_threshold(example_supplement(), "impact:d2", c(0, 20), conclusion = "total", target = 3)
  expect_equal(th$thresholds$threshold, 13.03899186228618, tolerance = 1e-6)
  expect_match(th$thresholds$description, "rises above")
  expect_error(cm_threshold(example_supplement(), "impact:d2", c(0, 20), conclusion = "total"),
               "target")
})

test_that("a change in ranking is found, by contribution", {
  th <- cm_threshold(example_supplement(), "assoc:d2:d3", c(0.1, 100), conclusion = "rank",
                     n_grid = 201)
  expect_equal(th$summary$item, c("d1 vs d2", "d1 vs d3", "d2 vs d3"))
  expect_equal(th$summary$n_thresholds, c(1L, 0L, 0L))
  expect_equal(nrow(th$thresholds), 1)
  t1 <- th$thresholds[th$thresholds$item == "d1 vs d2", ]
  expect_equal(t1$threshold, 10.20365117522034, tolerance = 1e-6)
  expect_match(t1$description, "d1 moves above d2")
  expect_equal(th$summary$result, c("threshold(s) found", "none found in the range",
                                    "none found in the range"))
})

test_that("unusable parts of the range are reported, not mistaken for thresholds", {
  # A probability of 1 or more is invalid. Below 1 every grid point passes
  # the triple screen (reference_v030.py port), so the only unusable region
  # is P(d1) >= 1.
  th <- cm_threshold(example_supplement(), "prob:d1", c(0.05, 1.5), conclusion = "rank")
  expect_equal(nrow(th$regions), 1)
  expect_equal(th$regions$status, "infeasible")
  expect_equal(th$regions$from, 1.007, tolerance = 1e-12)
  expect_equal(th$regions$to, 1.5)
  expect_equal(th$regions$n_points, 35L)
  expect_equal(sort(unique(th$scan$status)), c("infeasible", "ok"))
  expect_true(all(th$scan$status[th$scan$value >= 1] == "infeasible"))
  expect_true(all(th$scan$status[th$scan$value < 1] == "ok"))
  expect_equal(th$summary$n_thresholds, c(1L, 1L, 0L))
  expect_equal(th$thresholds$item, c("d1 vs d2", "d1 vs d3"))
  expect_equal(th$thresholds$status, c("threshold", "threshold"))
  expect_equal(th$thresholds$threshold, c(0.22362189686298367, 0.6220804451704026), tolerance = 1e-6)
  expect_equal(th$summary$result, c("threshold(s) found", "threshold(s) found",
                                    "none found where the range could be evaluated"))
  expect_output(print(th), "could not be used")
})

test_that("an invalid value of the input itself is an infeasible point", {
  # Odds ratios <= 0 on a linear grid.
  th <- cm_threshold(example_supplement(), "assoc:d1:d2", c(-1, 5), conclusion = "rank",
                     log_scale = FALSE)
  expect_equal(th$regions$status, "infeasible")
  expect_equal(th$regions$n_points, 17L)
  expect_equal(nrow(th$thresholds), 0)
  expect_true(all(th$summary$result == "none found where the range could be evaluated"))
})

test_that("interaction inputs use the global method", {
  skip_on_cran()
  th <- cm_threshold(example_supplement(), "inter:d1:d2", c(-1, 1), conclusion = "total",
                     target = 2.12, n_grid = 11)
  expect_equal(th$method, "global")
  expect_equal(th$baseline, 0)
  # The aggregate falls as the interaction grows (crude impacts absorb it).
  expect_equal(th$thresholds$status, "threshold")
  expect_equal(th$thresholds$threshold, -0.5347852259874342, tolerance = 1e-5)
  expect_match(th$thresholds$description, "falls below")
  expect_error(cm_threshold(example_supplement(), "inter:d1:d2", c(-1, 1), conclusion = "total",
                            target = 2.12, method = "simultaneous"),
               class = "deconflate_unsupported")
  expect_error(cm_threshold(example_supplement(), "three:d1:d2:d3", c(0.5, 2), conclusion = "total",
                            target = 2.12, method = "simultaneous"),
               class = "deconflate_unsupported")
})

test_that("inputs are validated", {
  m <- example_supplement()
  expect_error(cm_threshold(m, "assoc:d1", c(1, 2)), "wrong number of parts")
  expect_error(cm_threshold(m, "assoc:d1:zz", c(1, 2)), "Unknown diseases")
  expect_error(cm_threshold(m, "foo:d1", c(1, 2)), "Unknown input")
  expect_error(cm_threshold(m, "impact:d1", c(2, 1)), "increasing")
  expect_error(cm_threshold(m, "impact:d1", c(-1, 1), log_scale = TRUE), "positive range")
  expect_error(cm_threshold(m, "impact:d1", c(0, 1), conclusion = "rank", diseases = "d1"),
               "at least two")
  expect_error(cm_threshold(m, "impact:d1", c(0, 1), conclusion = "rank", diseases = c("d1", "d1")),
               "at least two")
  expect_error(cm_threshold(m, "assoc:d1:d1", c(1, 2)), "same disease twice")
  expect_error(cm_threshold(m, c("impact:d1", "impact:d2"), c(1, 2)), "single string")
  expect_error(cm_threshold(m, "impact:d1", c(0, 1), n_grid = NA), "n_grid")
  expect_error(cm_threshold(m, "impact:d1", c(0, 1), log_scale = NA), "log_scale")
  expect_error(cm_threshold(supp_model(c(0, 0, 0)), "impact:d1", c(0, 1), conclusion = "change",
                            target = 0.1), "zero")
  expect_error(cm_threshold(m, "prob:d1", c(0.05, 0.5), joint = fit_joint(m)), "fixed `joint`")
})

test_that("threshold plots draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  th <- cm_threshold(example_supplement(), "assoc:d2:d3", c(0.1, 100), conclusion = "rank")
  expect_invisible(plot(th))
})
