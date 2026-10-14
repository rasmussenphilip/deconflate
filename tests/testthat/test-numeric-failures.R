# Failures inside base R's numerical routines (solve(), uniroot()) come out
# as classed errors, so that draws, comparisons and screens handle them like
# any other failure.

ratio_population <- function(ratio = 1e30) {
  cm_population(cm_diseases(ids3, c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2", "d1"), c("d2", "d3", "d3", "d2"), c(2, 1, 3, ratio),
                                measure = c("OR", "OR", "OR", "OR_ratio"), among = c(NA, NA, NA, "d3")))
}

test_that("a singular system of exact equations is deconflate_singular", {
  expect_error(solve_exact(matrix(1, 2, 2), c(0, 1)), "Could not solve", class = "deconflate_singular")
  expect_equal(solve_exact(diag(2), c(1, 2)), c(1, 2))
})

test_that("the baseline hazard of non-finite multipliers is deconflate_nonfinite", {
  expect_error(baseline_hazard(c(1, Inf), c(0.5, 0.5), 0.2), class = "deconflate_nonfinite")
  expect_error(baseline_hazard(c(1, NaN), c(0.5, 0.5), 0.2), class = "deconflate_nonfinite")
  # A multiplier of a cell with no animals does not matter.
  expect_equal(baseline_hazard(c(1, Inf), c(1, 0), 0.2), -log(0.8), tolerance = 1e-10)
})

test_that("an odds ratio ratio that no distribution can have is deconflate_infeasible", {
  m <- cm_model(ratio_population(), cm_impacts(ids3, c(1, 2, 3)))
  expect_error(deconflate(m, n_draws = 0), "associations, row 4: no joint distribution",
               class = "deconflate_infeasible")
  expect_error(fit_joint(ratio_population()), class = "deconflate_infeasible")
})

test_that("draws that reach such a value are rejected, not fatal", {
  m <- cm_model(ratio_population(1.5), cm_impacts(ids3, c(1, 2, 3)),
                distributions = list("associations:4" = dist_lognormal(log(1.5), 40)))
  res <- suppressWarnings(deconflate(m, n_draws = 20, seed = 1))
  rej <- res$draws$rejections
  expect_gt(nrow(rej), 0)
  # (Some extreme draws also stop the joint fit converging.)
  expect_true(any(rej$type == "infeasible" & grepl("odds ratio ratio", rej$reason)))
  expect_false(any(rej$type == "error"))
  expect_lt(res$draws$n_rejected, res$draws$n_draws)
})
