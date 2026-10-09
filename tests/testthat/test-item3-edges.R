# Edge cases of the one-function interface (v0.4 item 3), from a review of
# the draws, the event model and the scenario helpers.

unc_model <- function() {
  cm_model(example_supplement(), example_supplement()$impacts,
           distributions = list("impact:d1" = dist_normal(2.5, 0.5), "assoc:d2:d3" = dist_lognormal_ci(3, 2, 4.5)))
}

test_that("the draw arguments are checked before anything is fitted", {
  m <- unc_model()
  expect_error(deconflate(m, n_draws = 1), "at least 2")
  expect_error(deconflate(m, n_draws = 2.5), "whole number")
  expect_error(deconflate(m, n_draws = -1), "whole number")
  expect_error(deconflate(m, n_draws = 10, seed = 1.5), "seed")
  expect_error(deconflate(m, n_draws = 10, seed = 1e12), "seed")
  expect_error(deconflate(m, n_draws = 10, sampling = "lhs", lhs_replicates = 1), "lhs_replicates")
  expect_error(deconflate(m, n_draws = 3, sampling = "lhs"), "at least 4 draws")
  expect_error(compare_methods(m, n_draws = 1), "at least 2")
  r <- deconflate(m, n_draws = 2, seed = 1)
  expect_equal(r$draws$n_draws, 2)
  expect_output(print(r), "seed 1")
})

test_that("too few accepted draws give no intervals, with a note", {
  m <- example_supplement()
  m$distributions <- list("prob:d1" = dist_uniform(1.1, 1.6))   # every draw is an invalid probability
  r <- deconflate(m, n_draws = 10, seed = 1)
  expect_null(r$adjusted$lower)
  expect_true(any(grepl("All 10 draws were rejected", r$notes)))
  cmp <- compare_methods(m, n_draws = 10, seed = 1)
  expect_true(any(grepl("All 10 draws were rejected", cmp$notes)))
  expect_output(print(cmp), "All 10 draws were rejected")
  # A high rejection share is reported in the notes.
  m$distributions <- list("prob:d1" = dist_uniform(0.05, 1.2))
  r2 <- deconflate(m, n_draws = 40, seed = 2)
  expect_gt(r2$draws$n_rejected, 4)
  expect_true(any(grepl("of the draws were rejected", r2$notes)))
  expect_true(any(grepl("of the draws were rejected", summary(r2)$notes)))
})

test_that("the seed also makes a sampled joint distribution reproducible", {
  skip_on_cran()
  m <- cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5)))
  a <- deconflate(m, backend = "sampled", n_samples = 2000, n_chains = 50, seed = 1, n_draws = 0)
  b <- deconflate(m, backend = "sampled", n_samples = 2000, n_chains = 50, seed = 1, n_draws = 0)
  expect_equal(a$adjusted$adjusted, b$adjusted$adjusted)
  expect_equal(a$unknown_pairs$fitted_or, b$unknown_pairs$fitted_or)
})

test_that("replacing the impacts or interactions of a model drops their distributions", {
  m <- unc_model()
  ev <- cm_model(m, cm_impacts(ids3, c(1.5, 2, 1.3), measure = "HR", estimand = "snapshot_crude"))
  expect_equal(names(ev$distributions), "assoc:d2:d3")
  # The same impacts keep theirs.
  expect_equal(names(cm_model(m, m$impacts)$distributions), names(m$distributions))
  mi <- cm_model(example_supplement(), example_supplement()$impacts, cm_interactions("d1", "d2", 1),
                 distributions = list("inter:d1:d2" = dist_normal(1, 0.2)))
  expect_null(cm_model(mi, mi$impacts)$distributions)
})

test_that("scenario setters fix the value they set", {
  mi <- cm_model(example_supplement(), example_supplement()$impacts, cm_interactions("d1", "d2", 1),
                 distributions = list("inter:d1:d2" = dist_normal(1, 0.2), "impact:d1" = dist_normal(2.5, 0.5)))
  si <- set_interaction(mi, "d1", "d2", 3)
  expect_equal(names(si$distributions), "impact:d1")
  mt <- cm_model(supp_population(three_way = cm_three_way("d1", "d2", "d3", 1.5)),
                 example_supplement()$impacts,
                 distributions = list("three:d1:d2:d3" = dist_lognormal_ci(1.5, 0.8, 2.8)))
  expect_null(set_three_way(mt, "d3", "d1", "d2", 2)$distributions)
})

test_that("negative risk differences give no spurious warnings", {
  m <- supp_hr_model(values = c(1.5, 2, -0.02), measure = c("HR", "HR", "RD"))
  expect_silent(r <- deconflate(m, event_model = TRUE, overall_risk = 0.25, warn = FALSE))
  expect_true(all(is.finite(r$adjusted$adjusted)))
  expect_lt(r$adjusted$adjusted[3], 1)
})

test_that("interactions cannot be combined with event impacts", {
  m <- supp_hr_model()
  m$interactions <- cm_interactions("d1", "d2", 1)
  expect_error(deconflate(m, event_model = TRUE, overall_risk = 0.25), "additive impacts only",
               class = "deconflate_unsupported")
  expect_error(screen_interactions(supp_hr_model(), values = 1), "additive impacts only")
  expect_error(cm_threshold(supp_hr_model(), "inter:d1:d2", c(-1, 1), conclusion = "total", target = 0.05,
                            event_model = TRUE, overall_risk = 0.25), "additive impacts only")
})
