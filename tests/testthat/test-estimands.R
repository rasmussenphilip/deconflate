# Estimands, identifiability and edge cases of the additive engine.
# Reference values: inst/validation/reference_v020.py.

t2_ids <- c("d1", "d2", "d3")
truth <- c(0.02, 0.04, 0.06)

t2_pop <- function(three_way = NULL) {
  cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)),
                three_way = three_way)
}
t2_model <- function(values, estimand = "crude", adjusted_for = NA_character_,
                     interactions = NULL, units = NULL) {
  cm_model(t2_pop(), cm_impacts(t2_ids, values, estimand = estimand, adjusted_for = adjusted_for,
                                units = units),
           interactions)
}
# Bob's infeasible triple: p = 0.5, ORs 20 (a:b), 20 (a:c), 0.05 (b:c).
t2_bad_pop <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

test_that("crude estimates: the simultaneous solution recovers the known truth", {
  pop <- t2_pop()
  raw <- simulate_raw_impacts(pop, stats::setNames(truth, t2_ids))
  expect_equal(raw$estimand, rep("crude", 3))
  res <- deconflate(cm_model(pop, raw))
  expect_equal(res$method, "simultaneous")
  expect_equal(res$adjusted$adjusted, truth, tolerance = 1e-9)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-12)
})

test_that("Bob's counterexample: an adjusted coefficient needs the projection mapping", {
  # d1's estimate is the coefficient of d1 in an additive regression adjusted
  # for d2 (in the same population); d2 and d3 are crude.
  raw <- c(0.018664808346, 0.054064376914, 0.066681750858)
  m <- t2_model(raw, estimand = c("adjusted_linear", "crude", "crude"),
                adjusted_for = c("d2", NA, NA))
  res <- deconflate(m)
  expect_equal(res$adjusted$adjusted, truth, tolerance = 1e-9)
  expect_equal(res$adjusted$estimand, c("adjusted_linear", "crude", "crude"))
  A <- res$conflation$A
  expect_equal(unname(A["d1", "d2"]), 0)
  expect_equal(unname(A["d1", "d1"]), 1)
  # Zeroing the adjusted-for entry of the crude row (the shortcut) is wrong.
  shortcut <- c(0.01866481, 0.04010287, 0.05998282)
  expect_gt(max(abs(res$adjusted$adjusted - shortcut)), 1e-3)
  # The global method gives the same answer.
  expect_equal(deconflate(m, method = "global")$adjusted$adjusted, truth, tolerance = 1e-8)
  # The raw value matches a weighted least-squares fit on the joint
  # distribution (simulate_raw_impacts).
  sim <- simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids),
                              estimand = c("adjusted_linear", "crude", "crude"),
                              adjusted_for = c("d2", NA, NA))
  expect_equal(sim$value, raw, tolerance = 1e-9)
})

test_that("treating an adjusted coefficient as crude gives the wrong answer", {
  raw <- c(0.018664808346, 0.054064376914, 0.066681750858)
  res <- deconflate(t2_model(raw))
  expect_gt(max(abs(res$adjusted$adjusted - truth)), 1e-3)
})

test_that("fully adjusted coefficients need no adjustment", {
  raw <- c(0.021, 0.039, 0.062)
  res <- deconflate(t2_model(raw, estimand = "adjusted_linear", adjusted_for = "all"))
  expect_equal(res$adjusted$adjusted, raw, tolerance = 1e-14)
  expect_equal(unname(res$conflation$A), diag(3))
  # Listing every other disease is the same as "all".
  res2 <- deconflate(t2_model(raw, estimand = "adjusted_linear",
                              adjusted_for = c("d2; d3", "d1; d3", "d1;d2")))
  expect_equal(res2$adjusted$adjusted, raw, tolerance = 1e-14)
})

test_that("mixed estimands recover the truth", {
  raw <- c(0.02, 0.054064376914, 0.066681750858)
  m <- t2_model(raw, estimand = c("adjusted_linear", "crude", "adjusted_linear"),
                adjusted_for = c("all", NA, "d1"))
  expect_equal(deconflate(m)$adjusted$adjusted, truth, tolerance = 1e-9)
  sim <- simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids),
                              estimand = c("adjusted_linear", "crude", "adjusted_linear"),
                              adjusted_for = c("all", NA, "d1"))
  expect_equal(sim$value, raw, tolerance = 1e-9)
})

test_that("unsupported estimands are rejected", {
  expect_error(cm_impacts(t2_ids, truth, adjusted_for = c("d2", NA, NA)),
               class = "deconflate_unsupported")
  expect_error(cm_impacts(t2_ids, truth, estimand = "adjusted_linear"), "adjusted_for")
  expect_error(cm_impacts(t2_ids, truth, estimand = "matched"))
  # Snapshot estimands are for event impacts (with a measure).
  expect_error(cm_impacts(t2_ids, truth, estimand = "snapshot_crude"), class = "deconflate_unsupported")
  expect_error(cm_model(t2_pop(),
                        cm_impacts(t2_ids, truth, estimand = c("adjusted_linear", "crude", "crude"),
                                   adjusted_for = c("d9", NA, NA))),
               "unknown")
  m <- t2_model(truth, estimand = c("adjusted_linear", "crude", "crude"),
                adjusted_for = c("d2", NA, NA))
  # The published approximation (internal adjust_impacts()) is for crude
  # estimates only, and is not a method of deconflate().
  expect_error(adjust_impacts(m, method = "published"), class = "deconflate_unsupported")
  expect_error(adjust_impacts(m, method = "published"), "crude estimates only")
  expect_error(deconflate(t2_model(c(1, 2, 3)), method = "published"), class = "deconflate_unsupported")
  expect_error(deconflate(t2_model(c(1, 2, 3)), method = "published"), "compare_methods")
})

test_that("a non-identifiable adjustment set is rejected", {
  # b is present exactly when a is: a regression on a and b cannot separate
  # them.
  pop <- cm_population(cm_diseases(c("a", "b"), c(0.2, 0.2)),
                       cm_associations("a", "b", 1, measure = "cond_prob"))
  adj <- cm_model(pop, cm_impacts(c("a", "b"), c(1, 2), estimand = c("adjusted_linear", "crude"),
                                  adjusted_for = c("b", NA)))
  expect_error(deconflate(adj), class = "deconflate_singular")
  crude <- cm_model(pop, cm_impacts(c("a", "b"), c(1, 2)))
  expect_error(deconflate(crude), class = "deconflate_singular")
})

test_that("zero impacts give a zero aggregate with undefined shares", {
  res <- deconflate(t2_model(c(0, 0, 0)))
  expect_equal(res$adjusted$adjusted, c(0, 0, 0))
  expect_equal(res$totals$adjusted_total, 0)
  expect_true(all(is.na(res$contributions$share)))
  expect_true(all(is.na(res$adjusted$change)))
  expect_equal(contribution_table(res)$contribution, c(0, 0, 0))
  # The legacy gap (reproduce_rasmussen_2022() only) of a zero aggregate.
  pg <- legacy_gap(res, 100, "decrease", "percent")
  expect_equal(pg$summary$gap, 0)
  expect_equal(pg$attribution$gap, c(0, 0, 0))
  pub <- adjust_impacts(t2_model(c(0, 0, 0)), method = "published")
  expect_equal(pub$adjusted$adjusted, c(0, 0, 0))
})

test_that("cancelling signed impacts are handled without dividing by the total", {
  pop <- cm_population(cm_diseases(c("a", "b"), c(0.3, 0.3)), cm_associations("a", "b", 1))
  res <- deconflate(cm_model(pop, cm_impacts(c("a", "b"), c(5, -5), units = "kg")))
  expect_equal(res$adjusted$adjusted, c(5, -5))
  expect_equal(res$totals$adjusted_total, 0)
  expect_true(all(is.na(res$contributions$share)))
  ct <- contribution_table(res)
  expect_equal(ct$contribution, c(1.5, -1.5))
  expect_true(all(is.na(ct$share)))
  # The legacy gap and value (reproduce_rasmussen_2022() only) attribute
  # the gap without dividing by the zero total.
  pg <- legacy_gap(res, 1000, "decrease", "absolute")
  expect_equal(pg$summary$disease_free, 1000)
  expect_equal(pg$attribution$gap, c(1.5, -1.5))
  expect_equal(legacy_value(pg, 2)$by_disease$value, c(3, -3))
  expect_equal(legacy_value(pg, 2)$value, 0)
})

test_that("a published denominator of zero gives a non-finite result and a warning", {
  pop <- cm_population(cm_diseases(c("a", "b"), c(0.1, 0.15)), cm_associations("a", "b", 2))
  E <- excess_matrix(pop)
  m <- cm_model(pop, cm_impacts(c("a", "b"), c(1, -1 / E["b", "a"])))
  expect_warning(adjust_impacts(m, method = "published"))
  res <- suppressWarnings(adjust_impacts(m, method = "published"))
  x <- res$adjusted$adjusted[1]
  expect_true(!is.finite(x) || abs(x) > 1e10)
})

test_that("jointly infeasible pairs are rejected although A is well conditioned", {
  m <- cm_model(t2_bad_pop(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(deconflate(m), class = "deconflate_infeasible")
  expect_error(adjust_impacts(m, method = "published"), class = "deconflate_infeasible")
  # Unchecked, the infeasible inputs flip the signs of b and c (warn = FALSE).
  res <- deconflate(m, feasibility = "none", warn = FALSE)
  expect_equal(res$diagnostics$sign_changes, "b, c")
  expect_equal(res$diagnostics$condition_number, 6.075710520908865, tolerance = 1e-8)
  expect_equal(res$diagnostics$feasibility, "unchecked")
  expect_error(suppressWarnings(deconflate(m, method = "global", max_iter = 300)),
               class = "deconflate_nonconvergence")
})

test_that("a stale joint distribution is rejected; an impact change is not", {
  m <- example_supplement()
  j <- fit_joint(m)
  expect_silent(deconflate(m, method = "global", joint = j))
  m2 <- set_association(m, "d1", "d2", 3)
  expect_error(deconflate(m2, method = "global", joint = j), class = "deconflate_error")
  m3 <- m
  m3$diseases$prob[1] <- 0.11
  expect_error(deconflate(m3, method = "global", joint = j), class = "deconflate_error")
  m4 <- cm_model(t2_pop(three_way = cm_three_way("d1", "d2", "d3", 2)), m$impacts)
  expect_error(deconflate(m4, method = "global", joint = j), class = "deconflate_error")
  # A pair made unknown changes the constrained pairs.
  m6 <- m
  m6$associations <- m$associations[-2, ]
  expect_error(deconflate(m6, method = "global", joint = j), "constrained pairs")
  m5 <- t2_model(c(1, 2, 3))
  expect_equal(deconflate(m5, method = "global", joint = j)$adjusted$adjusted,
               deconflate(m5)$adjusted$adjusted, tolerance = 1e-8)
})

test_that("interactions: the global method recovers the truth, and is used automatically", {
  # True impacts 0.02, 0.04, 0.06; interactions d1:d2 = 0.01, d2:d3 = 0.015.
  raw <- c(0.027146191901, 0.061369036759, 0.071166773469)
  int <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(0.01, 0.015))
  m <- t2_model(raw, interactions = int)
  expect_error(adjust_impacts(m, method = "simultaneous"), class = "deconflate_unsupported")
  expect_error(adjust_impacts(m, method = "published"), class = "deconflate_unsupported")
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "global")
  expect_equal(res$notes, "The global method was used because of interactions.")
  rs <- deconflate(m, method = "simultaneous", n_draws = 0)
  expect_equal(rs$method, "global")
  expect_match(rs$notes, "instead of the simultaneous method because of interactions")
  expect_equal(rs$adjusted$adjusted, res$adjusted$adjusted)
  expect_equal(res$adjusted$adjusted, truth, tolerance = 1e-8)
  expect_equal(res$totals$adjusted_total, 0.021095698976771753, tolerance = 1e-9)
  expect_equal(res$contributions$total, c(0.002122396963, 0.006547849488, 0.012425452526),
               tolerance = 1e-9)
  expect_equal(sum(res$contributions$total), res$totals$adjusted_total, tolerance = 1e-14)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-10)
  # raw = A b + offset
  expect_equal(as.vector(res$conflation$A %*% truth + res$conflation$offset), raw,
               tolerance = 1e-9)
})

test_that("interactions with adjusted estimands are recovered", {
  int <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(0.01, 0.015))
  est <- c("adjusted_linear", "crude", "adjusted_linear")
  adj <- c("d2", NA, "all")
  raw <- simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids),
                              interactions = int, estimand = est, adjusted_for = adj)
  res <- deconflate(cm_model(t2_pop(), raw, int))
  expect_equal(res$method, "global")
  expect_equal(res$adjusted$adjusted, truth, tolerance = 1e-8)
})

test_that("a three-way term is the log-linear three-way interaction", {
  pop <- t2_pop(three_way = cm_three_way("d1", "d2", "d3", 2.5))
  j <- fit_joint(pop)
  expect_true(j$converged)
  p <- j$prob
  x <- j$cells
  cell <- function(a, b, c) p[x[, 1] == a & x[, 2] == b & x[, 3] == c]
  lr <- log(cell(1, 1, 1) * cell(1, 0, 0) * cell(0, 1, 0) * cell(0, 0, 1) /
              (cell(1, 1, 0) * cell(1, 0, 1) * cell(0, 1, 1) * cell(0, 0, 0)))
  expect_equal(lr, log(2.5), tolerance = 1e-8)
  # Pairwise tables and marginals are still matched.
  j0 <- fit_joint(t2_pop())
  J <- crossprod(x * p, x)
  J0 <- crossprod(j0$cells * j0$prob, j0$cells)
  expect_equal(J, J0, tolerance = 1e-9)
  # Without interactions the additive result does not depend on it. A model
  # with three-way terms is adjusted with the global method.
  m <- cm_model(pop, cm_impacts(t2_ids, c(2.5, 5, 7.5)))
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "global")
  expect_equal(res$notes, "The global method was used because of three-way terms.")
  expect_equal(res$adjusted$adjusted, deconflate(example_supplement())$adjusted$adjusted, tolerance = 1e-8)
  expect_equal(global_reasons(set_interaction(m, "d1", "d2", 1)), c("interactions", "three-way terms"))
  # A ratio of 1 is maximum entropy.
  j1 <- fit_joint(t2_pop(three_way = cm_three_way("d1", "d2", "d3", 1)))
  expect_equal(j1$prob, j0$prob, tolerance = 1e-10)
})

test_that("metadata is kept and impacts are put in disease order", {
  imp <- cm_impacts(c("d3", "d1", "d2"), c(7.5, 2.5, 5), label = "yield", units = "%",
                    source = "test")
  m <- cm_model(t2_pop(), imp)
  expect_equal(m$impacts$disease, t2_ids)
  expect_equal(m$impacts$value, c(2.5, 5, 7.5))
  expect_equal(attr(m$impacts, "label"), "yield")
  expect_equal(attr(m$impacts, "units"), "%")
  expect_equal(attr(m$impacts, "kind"), "additive")
  res <- deconflate(m)
  expect_equal(res$label, "yield")
  expect_equal(res$units, "%")
  expect_output(print(res), "yield")
  expect_output(print(m), "cm_model")
  expect_output(print(m), "additive")
  # Impacts replaced after cm_model() are re-validated and reordered.
  m2 <- example_supplement()
  m2$impacts <- imp
  expect_equal(deconflate(m2)$adjusted$adjusted, res$adjusted$adjusted)
  m2$impacts <- cm_impacts(c("d1", "d2"), c(2.5, 5))
  expect_error(deconflate(m2), "No impact for")
})

test_that("one population serves several outcomes, one impact table each", {
  pop <- t2_pop()
  yield <- cm_model(pop, cm_impacts(t2_ids, c(2.5, 5, 7.5), label = "yield", units = "%"))
  welfare <- cm_model(pop, cm_impacts(t2_ids, c(1, 0, 3), label = "welfare", units = "score"))
  expect_identical(yield$diseases, welfare$diseases)
  expect_identical(yield$associations, welfare$associations)
  ry <- deconflate(yield)
  rw <- deconflate(welfare)
  expect_equal(ry$adjusted$adjusted, deconflate(example_supplement())$adjusted$adjusted)
  expect_equal(rw$units, "score")
  expect_equal(rw$label, "welfare")
  # One joint fit can serve both with the global method.
  j <- fit_joint(pop)
  expect_equal(deconflate(welfare, method = "global", joint = j)$adjusted$adjusted,
               rw$adjusted$adjusted, tolerance = 1e-8)
  expect_identical(deconflate(yield, method = "global", joint = j)$joint, j)
})

test_that("disease ids are validated", {
  expect_error(cm_diseases(c("a:b", "c"), c(0.1, 0.2)), "must not contain")
  expect_error(cm_diseases(c("a;b", "c"), c(0.1, 0.2)), "must not contain")
  expect_error(cm_diseases(c("all", "c"), c(0.1, 0.2)), "reserved")
  expect_error(cm_diseases(c("a", "a"), c(0.1, 0.2)), "unique")
  expect_error(cm_diseases(c("a", "b"), c(0.1, 1.2)), class = "deconflate_infeasible")
})

test_that("covariate-adjusted associations need an explicit choice", {
  d <- cm_diseases(c("a", "b"), c(0.2, 0.3))
  a <- cm_associations("a", "b", 2, adjusted = TRUE, adjusted_for = "parity")
  expect_error(cm_population(d, a), class = "deconflate_unsupported")
  pop <- cm_population(d, a, adjusted_associations = "use_as_marginal")
  expect_s3_class(pop, "cm_population")
  expect_output(print(pop), "used as marginal")
})

test_that("retired association measures are rejected with what to do instead", {
  for (meas in c("independent", "unknown", "table")) {
    expect_error(cm_associations("a", "b", 1, measure = meas), class = "deconflate_unsupported")
  }
  # An unrelated pair is an odds ratio of 1; an unknown pair has no row.
  pop <- cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                       cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)))
  expect_equal(pair_tables(pop)$status, c("specified", "unknown", "specified"))
  expect_equal(pair_tables(t2_pop())$status, rep("specified", 3))
})
