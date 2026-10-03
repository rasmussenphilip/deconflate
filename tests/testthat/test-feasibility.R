bad_model <- function() {
  cm_model(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
           cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

test_that("the triple screen detects an infeasible triple", {
  f <- check_feasibility(bad_model(), method = "triples")
  expect_false(f$feasible)
  expect_equal(f$triples$gap, 0.2258840035526648, tolerance = 1e-9)
  ok <- check_feasibility(example_supplement(), method = "triples")
  expect_equal(nrow(ok$triples), 0)
})

test_that("the LP check agrees and identifies the conflict", {
  skip_if_not_installed("lpSolve")
  f <- check_feasibility(bad_model(), method = "lp")
  expect_false(f$feasible)
  expect_equal(f$lp$objective, 0.2258840035526648, tolerance = 1e-7)
  expect_true(check_feasibility(example_global_dairy(), method = "lp")$feasible)
  expect_output(print(f), "NOT jointly feasible")
})

test_that("the global method refuses infeasible associations", {
  expect_error(suppressWarnings(deconflate(
    cm_model(bad_model()$diseases, bad_model()$associations,
             cm_impacts(c("a", "b", "c"), c(0.01, 0.02, 0.03))),
    method = "global", max_iter = 500)), class = "deconflate_infeasible")
})
