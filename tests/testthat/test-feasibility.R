# Joint feasibility of pairwise associations. Bob's counterexample: three
# diseases with p = 0.5 and odds ratios 20 (a:b), 20 (a:c) and 0.05 (b:c).
# Each 2x2 table is valid, but no population has all three
# (inst/validation/reference_v020.py).

test_that("the triple screen detects an infeasible triple", {
  f <- check_feasibility(bad_population(), method = "triples")
  expect_false(f$feasible)
  expect_equal(f$triples$gap, 0.2258840035526648, tolerance = 1e-9)
  ok <- check_feasibility(example_supplement(), method = "triples")
  expect_equal(nrow(ok$triples), 0)
  expect_true(is.na(ok$feasible))
})

test_that("the LP check agrees and identifies the conflict", {
  skip_if_not_installed("lpSolve")
  f <- check_feasibility(bad_population(), method = "lp")
  expect_false(f$feasible)
  expect_equal(f$lp$objective, 0.2258840035526648, tolerance = 1e-7)
  expect_true(check_feasibility(example_global_dairy()$population, method = "lp")$feasible)
  expect_output(print(f), "NOT jointly feasible")
})

test_that("the pairwise methods screen triples by default", {
  m <- cm_model(bad_population(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(deconflate(m), class = "deconflate_infeasible")
  expect_error(deconflate(m, method = "published"), class = "deconflate_infeasible")
  expect_warning(res <- deconflate(m, feasibility = "none"), class = "deconflate_sign_change")
  expect_s3_class(res, "cm_result")
  skip_if_not_installed("lpSolve")
  expect_error(deconflate(m, feasibility = "lp"), class = "deconflate_infeasible")
  res <- deconflate(example_supplement(), feasibility = "lp")
  expect_equal(res$diagnostics$feasibility, "LP check passed (jointly feasible)")
})

test_that("the global method stops when the joint distribution does not converge", {
  m <- cm_model(bad_population(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(suppressWarnings(deconflate(m, method = "global", max_iter = 300)),
               class = "deconflate_nonconvergence")
  expect_warning(j <- fit_joint(m, max_iter = 300), class = "deconflate_nonconvergence")
  expect_false(j$converged)
  expect_error(deconflate(m, method = "global", joint = j), class = "deconflate_nonconvergence")
})
