# cm_threshold_prob() (v0.6 item 12; design note section 19). Reference
# values: inst/validation/reference_v06_threshold_prob.py (a closed form:
# one normal input, losses linear in it).

prob_model <- function() {
  cm_model(example_supplement(), example_supplement()$impacts,
           distributions = list("impacts:1" = dist_normal(2.5, 0.5)))
}

test_that("probabilities match the closed form", {
  res <- deconflate(prob_model(), n_draws = 0)
  p <- cm_threshold_prob(res, "pair:d2:d3", c(1, 9), conclusion = "total", target = 2.15, n_draws = 1000,
                         n_grid = 9, log_scale = FALSE, sampling = "lhs", seed = 1)
  expect_s3_class(p, "cm_threshold_prob")
  cv <- p$curves
  expect_equal(cv$value, 1:9)
  expect_true(all(cv$n == 1000))
  ref <- c(`1` = 0.0000000003, `2` = 0.0562151686, `3` = 0.8128990229, `5` = 0.9998273908, `8` = 0.9999999858)
  # One normal input drawn by Latin hypercube: each probability within 0.01.
  expect_lt(max(abs(cv$probability[match(as.numeric(names(ref)), cv$value)] - ref)), 0.01)
  # The 0.5 crossing lies between 2 and 3: interpolating the reference
  # probabilities gives 2.5865 (the exact crossing, the deterministic
  # threshold, is 2.5810).
  c50 <- p$crossings[p$crossings$probability == 0.5, ]
  expect_equal(nrow(c50), 1L)
  expect_equal(c(c50$lower, c50$upper), c(2, 3))
  expect_equal(c50$direction, "rises")
  expect_lt(abs(c50$threshold - 2.5865), 0.02)
  det <- cm_threshold(res, "pair:d2:d3", c(1, 9), conclusion = "total", target = 2.15)
  expect_equal(det$thresholds$threshold, 2.5810100885, tolerance = 1e-6)
  # A set of diseases: the loss of d1 and d2.
  q <- cm_threshold_prob(res, "pair:d2:d3", c(1, 9), conclusion = "total", diseases = c("d1", "d2"), target = 0.75,
                         n_draws = 1000, n_grid = 9, log_scale = FALSE, sampling = "lhs", seed = 1)
  ref2 <- c(`1` = 0.0000396122, `2` = 0.1256310410, `3` = 0.7321594721, `5` = 0.9984189712, `8` = 0.9999999154)
  expect_lt(max(abs(q$curves$probability[match(as.numeric(names(ref2)), q$curves$value)] - ref2)), 0.01)
})

test_that("the draws are reproducible, and the varied input is not drawn", {
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impacts:1" = dist_normal(2.5, 0.5), "impacts:2" = dist_normal(5, 1)))
  res <- deconflate(m, n_draws = 0)
  a <- cm_threshold_prob(res, "impacts:2", c(3, 7), conclusion = "rank", diseases = c("d1", "d2"),
                         n_draws = 20, n_grid = 3, seed = 4)
  b <- cm_threshold_prob(res, "impacts:2", c(3, 7), conclusion = "rank", diseases = c("d1", "d2"),
                         n_draws = 20, n_grid = 3, seed = 4)
  expect_identical(a$curves, b$curves)
  # The curve moves with the varied input: d1 (2.5 x 0.10) ranks above d2
  # only when d2's impact is low.
  expect_gt(a$curves$probability[1], a$curves$probability[3])
  # The order of `diseases` sets the order of the pair.
  ba <- cm_threshold_prob(res, "impacts:2", c(3, 7), conclusion = "rank", diseases = c("d2", "d1"),
                          n_draws = 20, n_grid = 3, seed = 4)
  expect_equal(unique(ba$curves$item), "d2 vs d1")
  expect_equal(ba$curves$probability, 1 - a$curves$probability)
  # With only the varied input uncertain, there is nothing to draw.
  m1 <- prob_model()
  expect_error(cm_threshold_prob(deconflate(m1, n_draws = 0), "impacts:1", c(2, 3), conclusion = "total",
                                 target = 2, n_draws = 10, n_grid = 3, seed = 1),
               "No input other than impacts:1")
  expect_error(cm_threshold_prob(example_supplement(), "pair:d1:d2", c(1, 3), n_draws = 10, n_grid = 3),
               "No input has a distribution")
})

test_that("crossings are interpolated in both directions", {
  cr <- prob_crossings(c(1, 2, 3, 4), c(0.1, 0.6, 0.9, 0.4), c(0.5, 0.8), FALSE)
  expect_equal(cr$probability, c(0.5, 0.5, 0.8, 0.8))
  expect_equal(cr$direction, c("rises", "falls", "rises", "falls"))
  expect_equal(cr$threshold, c(1 + 0.4 / 0.5, 3 + 0.4 / 0.5, 2 + 0.2 / 0.3, 3 + 0.1 / 0.5))
  # On a log grid, on the log scale.
  cl <- prob_crossings(c(1, 10), c(0, 1), 0.5, TRUE)
  expect_equal(cl$threshold, sqrt(10))
  expect_equal(nrow(prob_crossings(c(1, 2), c(0.2, 0.3), 0.5, FALSE)), 0L)
})

test_that("automatic choice, printing and arguments", {
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impacts:1" = dist_normal(2.5, 0.5), "impacts:2" = dist_normal(5, 1)))
  res <- deconflate(m, n_draws = 0)
  p <- cm_threshold_prob(res, conclusion = "rank", diseases = c("d1", "d2"), top = 2, n_draws = 20, n_grid = 3,
                         seed = 2)
  expect_true(p$auto)
  expect_equal(nrow(p$inputs), 2L)
  out <- capture.output(print(p))
  expect_true(all(nchar(out) <= 80))
  expect_match(paste(out, collapse = " "), "P(d1 ranks above d2 by contribution)", fixed = TRUE)
  expect_error(cm_threshold_prob(res, "impacts:2", c(3, 7), levels = 1.2), "levels")
  expect_error(cm_threshold_prob(res, "impacts:2", c(3, 7), n_draws = 5), "n_draws")
  expect_error(cm_threshold_prob(res, "impacts:2", c(3, 7), n_grid = 2), "n_grid")
  expect_error(cm_threshold_prob(res, "impacts:2", c(3, 7), conclusion = "change", target = 0.1),
               class = "deconflate_unsupported")
})

test_that("inputs that add a row are refitted at each value", {
  # An unknown pair (d1-d3 has no row): the varied input adds a crude row.
  m <- cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5)),
                distributions = list("impacts:1" = dist_normal(2.5, 0.5), "impacts:2" = dist_normal(5, 1)))
  res <- deconflate(m, n_draws = 0)
  p <- cm_threshold_prob(res, "pair:d1:d3", c(0.5, 4), conclusion = "total", target = 1.9, n_draws = 20,
                         n_grid = 3, seed = 3)
  expect_true(all(p$curves$n == 20))
  expect_true(all(p$curves$probability >= 0 & p$curves$probability <= 1))
  # At the model's own value of an entered pair, the deterministic answer
  # is at the median: with no other uncertainty than impacts, P(total < its
  # central value) is near 0.5.
  r2 <- deconflate(prob_model(), n_draws = 0)
  pm <- cm_threshold_prob(r2, "pair:d2:d3", c(3, 3.0001), conclusion = "total",
                          target = r2$totals$adjusted_total, n_draws = 200, n_grid = 3, sampling = "lhs", seed = 5)
  expect_lt(abs(pm$curves$probability[1] - 0.5), 0.02)
})

test_that("relative targets need a positive reference", {
  neg <- supp_model(c(-2.5, -5, -7.5))
  expect_error(cm_threshold(neg, "impacts:1", c(-4, -1), conclusion = "total", target = 0.9,
                            relative_to = "adjusted"), "negative", class = "deconflate_nonfinite")
  expect_error(cm_threshold(example_supplement(), conclusion = "rank", diseases = "d1"), "at least two")
})
