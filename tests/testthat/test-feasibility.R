# Joint feasibility of pairwise associations. Bob's counterexample: three
# diseases with p = 0.5 and odds ratios 20 (a:b), 20 (a:c) and 0.05 (b:c).
# Each 2x2 table is valid, but no population has all three
# (inst/validation/reference_v020.py).

t2_bad_pop <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}
t2_bad_model <- function() cm_model(t2_bad_pop(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))

test_that("the triple screen detects an infeasible triple", {
  f <- check_feasibility(t2_bad_pop(), method = "triples")
  expect_false(f$feasible)
  expect_equal(f$triples$gap, 0.2258840035526648, tolerance = 1e-9)
  expect_equal(c(f$triples$disease1, f$triples$disease2, f$triples$disease3), c("a", "b", "c"))
  ok <- check_feasibility(example_supplement(), method = "triples")
  expect_equal(nrow(ok$triples), 0)
  expect_true(is.na(ok$feasible))
  expect_output(print(ok), "necessary condition only")
})

test_that("triples with an unknown pair are not screened", {
  # Without the b:c association, the triple is unconstrained: no conflict.
  pop <- t2_bad_pop()
  pop$associations <- pop$associations[1:2, ]
  expect_equal(pair_tables(pop)$status, c("specified", "specified", "unknown"))
  f <- check_feasibility(pop, method = "triples")
  expect_equal(nrow(f$triples), 0)
  expect_true(is.na(f$feasible))
  skip_if_not_installed("lpSolve")
  expect_true(check_feasibility(pop, method = "lp")$feasible)
  # No associations at all: nothing to check.
  none <- cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)))
  expect_true(check_feasibility(none, method = "lp")$feasible)
})

test_that("the LP check agrees and identifies the conflict", {
  skip_if_not_installed("lpSolve")
  f <- check_feasibility(t2_bad_pop(), method = "lp")
  expect_false(f$feasible)
  expect_equal(f$lp$objective, 0.2258840035526648, tolerance = 1e-7)
  expect_equal(f$method, "lp + triples")
  # The global dairy example has unknown pairs: only the given pairs are
  # constrained.
  expect_true(check_feasibility(example_global_dairy(), method = "lp")$feasible)
  expect_output(print(f), "NOT jointly feasible")
})

test_that("the pairwise methods screen triples by default", {
  m <- t2_bad_model()
  expect_error(deconflate(m), class = "deconflate_infeasible")
  # The published approximation (internal; compare_methods() and reproduce)
  # is screened too.
  expect_error(adjust_impacts(m, method = "published"), class = "deconflate_infeasible")
  expect_warning(res <- deconflate(m, feasibility = "none"), class = "deconflate_sign_change")
  expect_s3_class(res, "cm_result")
  expect_equal(res$diagnostics$feasibility, "unchecked")
  res_s <- deconflate(example_supplement())
  expect_equal(res_s$diagnostics$feasibility, "triple screen passed (necessary condition only)")
  skip_if_not_installed("lpSolve")
  expect_error(deconflate(m, feasibility = "lp"), class = "deconflate_infeasible")
  res <- deconflate(example_supplement(), feasibility = "lp")
  expect_equal(res$diagnostics$feasibility, "LP check passed (jointly feasible)")
})

test_that("the global method stops when the joint distribution does not converge", {
  m <- t2_bad_model()
  expect_error(suppressWarnings(deconflate(m, method = "global", max_iter = 300)),
               class = "deconflate_nonconvergence")
  expect_warning(j <- fit_joint(m, max_iter = 300), class = "deconflate_nonconvergence")
  expect_false(j$converged)
  expect_error(deconflate(m, method = "global", joint = j), class = "deconflate_nonconvergence")
  # So does the event model, which always uses the joint distribution.
  ev <- cm_model(t2_bad_pop(), cm_impacts(c("a", "b", "c"), c(1.5, 2, 1.3), measure = "HR",
                                          estimand = "snapshot_crude"))
  expect_error(suppressWarnings(deconflate(ev, event_model = TRUE, overall_risk = 0.25, max_iter = 300)),
               class = "deconflate_nonconvergence")
})
