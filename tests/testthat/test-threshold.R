# Threshold searches (cm_threshold). Reference values:
# inst/validation/reference_v030.py.

test_that("a relative change in the aggregate is found for an association", {
  m <- example_supplement()
  # d1:d3 is independent in the example (baseline OR 1).
  th <- cm_threshold(m, "assoc:d1:d3", c(1, 50), conclusion = "change", target = -0.1)
  expect_s3_class(th, "cm_threshold")
  expect_equal(th$baseline, 1)
  expect_true(th$log_scale)
  expect_equal(th$method, "pairwise shortcut")
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
  expect_match(th$thresholds$description, "moves above no effect")
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

# p = 0.5 for a, b and c with phi(a, b) = phi(a, c) = 0.8: the conflation
# matrix is the correlation matrix, singular where phi(b, c) = 2 * 0.8^2 - 1
# = 0.28, i.e. at OR(b, c) = (1.28 / 0.72)^2 = 256 / 81. With raw impacts
# 1, 2, 3 the adjusted impact of a is (s - 3) / (s - 0.28) (s = phi(b, c)),
# which changes sign only at the pole. Below OR(b, c) = 16 (phi 0.6) the
# pairs are jointly infeasible, so the pole lies in the infeasible region.
# Reference values: python emulation of cm_threshold() with R's LINPACK
# rank check (qr tolerance 1e-7).
singular_pole_model <- function() {
  cm_model(cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                         cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(0.8, 0.8, 20),
                                         measure = c("phi", "phi", "OR"))),
           cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
}

test_that("an infeasible part of the range is a region, not a crossing", {
  m <- singular_pole_model()
  # The joint fit reports the region below the pole as infeasible.
  th0 <- cm_threshold(m, "assoc:b:c", c(0.05, 100), conclusion = "sign", diseases = "a",
                      n_grid = 3)
  expect_equal(th0$scan$status, c("infeasible", "infeasible", "ok"))
  expect_equal(nrow(th0$thresholds), 0)
  expect_equal(th0$regions$n_points, 2L)
  expect_true(th0$log_scale)
  expect_equal(th0$baseline, 20)
  # The feasibility argument is gone.
  expect_error(cm_threshold(m, "assoc:b:c", c(0.05, 100), feasibility = "none"),
               class = "deconflate_unsupported")
})

# Crude rows for a, for b and for the combination a;b (all vs none) do not
# tell b_a, b_b and the interaction apart when P(a) + P(b) = 1, whatever the
# association: with P(b) = 0.5, the adjusted impacts have a pole at
# P(a) = 0.5 and change sign across it.
pole_rows_model <- function(pa = 0.3) {
  cm_model(cm_population(cm_diseases(c("a", "b"), c(pa, 0.5)), cm_associations("a", "b", 2)),
           cm_impacts(c("a", "b", "a;b"), c(1, 2, 4)))
}

test_that("a pole of a singular system is a discontinuity, not a threshold", {
  expect_error(deconflate(pole_rows_model(0.5), n_draws = 0), class = "deconflate_singular")
  near <- suppressWarnings(deconflate(pole_rows_model(0.499), n_draws = 0))$adjusted$adjusted
  expect_lt(near[1], -1000)
  # The grid misses the pole; the bisection lands on it.
  th <- cm_threshold(pole_rows_model(), "prob:a", c(0.3, 0.7), conclusion = "sign", diseases = "a",
                     n_grid = 10)
  expect_equal(th$method, "full model")
  expect_false(th$log_scale)
  expect_equal(th$scan$status, rep("ok", 10))
  st <- th$thresholds
  expect_equal(nrow(st), 1)
  expect_equal(st$status, "discontinuity")
  expect_match(st$description, "nearly singular")
  expect_true(is.na(st$threshold))
  expect_lt(st$lower, 0.5)
  expect_gt(st$upper, 0.5)
  expect_lt(st$below, -10)
  expect_gt(st$above, 10)
  expect_equal(th$summary$n_thresholds, 0L)
  expect_equal(th$summary$result, "no threshold; see discontinuities or unresolved crossings")
  # A grid point on the pole is unusable: a region, and no crossing.
  th2 <- cm_threshold(pole_rows_model(), "prob:a", c(0.3, 0.7), conclusion = "sign", diseases = "a",
                      n_grid = 11)
  expect_equal(th2$scan$status[6], "singular")
  expect_equal(th2$regions$status, "singular")
  expect_equal(nrow(th2$thresholds), 0)
  expect_equal(th2$summary$result, "none found where the range could be evaluated")
})

test_that("cm_threshold rejects the retired `method` argument", {
  m <- example_supplement()
  expect_error(cm_threshold(m, "impact:d1", c(-3, 3), conclusion = "sign", diseases = "d1",
                            method = "global"),
               "`method` argument has been removed", class = "deconflate_unsupported")
  # The model stays finite at impact:d1 = -0.5266329181546854, where eq. 16
  # of the 2022 paper divides by zero.
  pole <- -0.5266329181546854
  expect_true(all(is.finite(deconflate(supp_model(c(pole, 5, 7.5)), warn = FALSE)$adjusted$adjusted)))
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

test_that("interaction inputs fit the full model", {
  skip_on_cran()
  th <- cm_threshold(example_supplement(), "inter:d1:d2", c(-1, 1), conclusion = "total",
                     target = 2.12, n_grid = 11)
  # Every point (with or without the interaction row) uses one joint fit.
  expect_equal(th$method, "full model")
  expect_equal(th$baseline, 0)
  # The aggregate falls as the interaction grows (crude impacts absorb it).
  expect_equal(th$thresholds$status, "threshold")
  expect_equal(th$thresholds$threshold, -0.5347852259874342, tolerance = 1e-5)
  expect_match(th$thresholds$description, "falls below")
  # Ranked by contribution (interaction shares included).
  thr <- cm_threshold(example_supplement(), "inter:d1:d2", c(-1, 1), conclusion = "rank", n_grid = 5)
  expect_equal(thr$summary$item, c("d1 vs d2", "d1 vs d3", "d2 vs d3"))
  # An interaction of event impacts is a hazard ratio (an HR row adjusted
  # for all), so the range must be positive.
  expect_error(cm_threshold(supp_hr_model(), "inter:d1:d2", c(-1, 1), conclusion = "total",
                            target = 0.1, event_model = TRUE, overall_risk = 0.25),
               class = "deconflate_invalid_input")
})

test_that("rows are varied by number and by alias", {
  m <- example_supplement()
  a <- cm_threshold(m, "assoc:d2:d3", c(0.2, 20), conclusion = "change", target = 0.1, n_grid = 21)
  b <- cm_threshold(m, "associations:3", c(0.2, 20), conclusion = "change", target = 0.1, n_grid = 21)
  expect_equal(a$thresholds$threshold, b$thresholds$threshold)
  expect_equal(a$baseline, b$baseline)
  i1 <- cm_threshold(m, "impact:d2", c(0, 20), conclusion = "total", target = 3, n_grid = 11)
  i2 <- cm_threshold(m, "impacts:2", c(0, 20), conclusion = "total", target = 3, n_grid = 11)
  expect_equal(i1$thresholds$threshold, 13.03899186228618, tolerance = 1e-6)
  expect_equal(i2$thresholds$threshold, i1$thresholds$threshold)
  # A three-way ratio changes additive results only through combinations:
  # with an interaction it moves the total.
  mi <- set_interaction(m, "d1", "d2", 1)
  t3 <- cm_threshold(mi, "three:d1:d2:d3", c(0.5, 2), conclusion = "total", target = 2.0890878497857703,
                     n_grid = 5)
  expect_true(is.na(t3$baseline))   # no OR_ratio row in the model
  expect_equal(t3$scan$total[c(1, 5)], c(2.0913976662985534, 2.086745618267496), tolerance = 1e-7)
  # Two crude rows of a pair: assoc: is ambiguous; vary a row by number.
  m5 <- read_five()
  expect_error(cm_threshold(m5, "assoc:LAM:MAS", c(1, 3)), "associations:")
  expect_error(cm_threshold(m, "impacts:9", c(0, 1)), "row")
})

test_that("descriptions of relative changes give the level", {
  th <- cm_threshold(example_supplement(), "assoc:d2:d3", c(0.2, 20), conclusion = "change", target = 0.1)
  expect_match(th$thresholds$description, "the aggregate (rises above|falls below) the baseline \\+ 10%")
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
