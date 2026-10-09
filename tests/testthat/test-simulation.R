# Validation against simulated populations with known effects
# (R/simulate.R). True impacts 2, 4, 6 (any units; here percent). Reference
# values from reference_v020_tests.py and reference_v020.py (the latter uses
# 0.02, 0.04, 0.06: the values here are 100 times its values).

supp_truth <- c(d1 = 2, d2 = 4, d3 = 6)
supp_raw_crude <- c(2.42130633452375, 5.40643769141281, 6.66817508583985)

test_that("simulated crude raw impacts match the Python reference", {
  pop <- supp_population()
  raw <- simulate_raw_impacts(pop, supp_truth)
  expect_s3_class(raw, "cm_impacts")
  expect_equal(raw$disease, ids3)
  expect_equal(raw$value, supp_raw_crude, tolerance = 1e-9)
  expect_equal(raw$estimand, rep("crude", 3))
  expect_equal(raw$source, rep("simulated (expected)", 3))
  # The order of true_impacts does not matter; a cm_model works too.
  expect_equal(simulate_raw_impacts(example_supplement(), rev(supp_truth))$value, raw$value)
  expect_equal(attr(simulate_raw_impacts(pop, supp_truth, units = "kg"), "units"), "kg")
  expect_error(simulate_raw_impacts(pop, c(2, 4, 6)), "named")
  expect_error(simulate_raw_impacts(pop, c(d1 = 2, d2 = 4)), "named")
})

test_that("simultaneous and global methods recover the truth; published does not", {
  pop <- supp_population()
  m <- cm_model(pop, simulate_raw_impacts(pop, supp_truth))
  expect_equal(deconflate(m)$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
  expect_equal(deconflate(m, method = "global")$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
  expect_lt(deconflate(m)$diagnostics$max_reconstruction_residual, 1e-10)
  pub <- adjust_impacts(m, method = "published")
  expect_gt(max(abs(pub$adjusted$adjusted - supp_truth)), 0.1)
})

test_that("adjusted_linear estimands are simulated and recovered", {
  pop <- supp_population()
  # d1 from a regression adjusted for d2; d2 and d3 crude.
  raw <- simulate_raw_impacts(pop, supp_truth, estimand = c("adjusted_linear", "crude", "crude"),
                              adjusted_for = c("d2", NA, NA))
  expect_equal(raw$estimand, c("adjusted_linear", "crude", "crude"))
  expect_equal(raw$adjusted_for, c("d2", NA, NA))
  expect_equal(raw$value, c(1.8664808346, 5.4064376914, 6.6681750858), tolerance = 1e-9)
  res <- deconflate(cm_model(pop, raw))
  expect_equal(res$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
  expect_equal(res$adjusted$estimand, raw$estimand)
  expect_equal(deconflate(cm_model(pop, raw), method = "global")$adjusted$adjusted, unname(supp_truth),
               tolerance = 1e-8)
  # Treating the adjusted estimate as crude does not recover the truth.
  raw_as_crude <- cm_impacts(ids3, raw$value)
  expect_gt(abs(deconflate(cm_model(pop, raw_as_crude))$adjusted$adjusted[1] - 2), 0.1)

  # Adjusted for all other diseases: the coefficients are the true impacts.
  raw_all <- simulate_raw_impacts(pop, supp_truth, estimand = "adjusted_linear", adjusted_for = "all")
  expect_equal(raw_all$value, unname(supp_truth), tolerance = 1e-9)
  res_all <- deconflate(cm_model(pop, raw_all))
  expect_equal(res_all$adjusted$adjusted, unname(supp_truth), tolerance = 1e-9)
  expect_equal(unname(res_all$conflation$A), diag(3))

  # Mixed: d1 adjusted for all, d2 crude, d3 adjusted for d1.
  raw_mix <- simulate_raw_impacts(pop, supp_truth, estimand = c("adjusted_linear", "crude", "adjusted_linear"),
                                  adjusted_for = c("all", NA, "d1"))
  expect_equal(raw_mix$value, c(2, 5.4064376914, 6.6681750858), tolerance = 1e-9)
  expect_equal(deconflate(cm_model(pop, raw_mix))$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
})

test_that("results are in the units entered and rescale with them", {
  pop <- supp_population()
  raw <- simulate_raw_impacts(pop, supp_truth)
  raw_k <- simulate_raw_impacts(pop, 1000 * supp_truth, units = "kg")
  expect_equal(raw_k$value, 1000 * raw$value, tolerance = 1e-10)
  res <- deconflate(cm_model(pop, raw))
  res_k <- deconflate(cm_model(pop, raw_k))
  expect_equal(res_k$adjusted$adjusted, 1000 * res$adjusted$adjusted, tolerance = 1e-10)
  expect_equal(res_k$adjusted$adjusted, 1000 * unname(supp_truth), tolerance = 1e-8)
  expect_equal(res_k$totals$adjusted_total, 1000 * res$totals$adjusted_total, tolerance = 1e-10)
  expect_equal(res_k$units, "kg")
  pub <- adjust_impacts(cm_model(pop, raw), method = "published")
  pub_k <- adjust_impacts(cm_model(pop, raw_k), method = "published")
  expect_equal(pub_k$adjusted$adjusted, 1000 * pub$adjusted$adjusted, tolerance = 1e-10)
})

test_that("the global method recovers the truth with interactions", {
  pop <- supp_population()
  ints <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(1, 1.5))
  raw <- simulate_raw_impacts(pop, supp_truth, interactions = ints)
  expect_equal(raw$value, c(2.71461919006052, 6.13690367592731, 7.11667734693017), tolerance = 1e-9)
  m_int <- cm_model(pop, raw, ints)
  res <- deconflate(m_int, method = "global")
  expect_equal(res$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-10)
  expect_equal(res$totals$adjusted_total, 2.1095698976771753, tolerance = 1e-8)
  expect_equal(res$contributions$total, c(0.2122396963, 0.6547849488, 1.2425452526), tolerance = 1e-8)
  expect_equal(sum(res$contributions$total), res$totals$adjusted_total)
  expect_error(deconflate(m_int), class = "deconflate_unsupported")
  expect_error(deconflate(m_int, method = "published"), class = "deconflate_unsupported")
  # The internal approximation rejects interactions too.
  expect_error(adjust_impacts(m_int, method = "published"), "Interactions require",
               class = "deconflate_unsupported")

  # With an adjusted estimand as well.
  raw_adj <- simulate_raw_impacts(pop, supp_truth, interactions = ints,
                                  estimand = c("adjusted_linear", "crude", "crude"),
                                  adjusted_for = c("d2", NA, NA))
  res_adj <- deconflate(cm_model(pop, raw_adj, ints), method = "global")
  expect_equal(res_adj$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
})

test_that("Shapley shares sum to the total and split interactions equally", {
  pop <- supp_population()
  ints <- cm_interactions("d1", "d2", 1)
  m_int <- cm_model(pop, simulate_raw_impacts(pop, supp_truth, ints), ints)
  res <- deconflate(m_int, method = "global")
  b <- attribute_burden(res)
  P12 <- res$joint_pairs["d1", "d2"]
  expect_equal(b$interaction, c(0.5 * P12, 0.5 * P12, 0), tolerance = 1e-12)
  L <- sum(supp_truth * c(0.10, 0.15, 0.20)) + P12
  expect_equal(sum(b$total), L, tolerance = 1e-8)
  expect_equal(sum(b$share), 1)
})

test_that("a finite sample adds sampling error only", {
  pop <- supp_population()
  set.seed(11)
  raw_n <- simulate_raw_impacts(pop, supp_truth, n = 5000)
  expect_s3_class(raw_n, "cm_impacts")
  expect_equal(raw_n$source, rep("simulated (n = 5000)", 3))
  expect_true(all(is.finite(raw_n$value)))
  expect_lt(max(abs(raw_n$value - supp_raw_crude)), 1)
  res <- deconflate(cm_model(pop, raw_n))
  expect_lt(max(abs(res$adjusted$adjusted - supp_truth)), 1)
  # Reproducible with a seed; a supplied joint gives the same draws.
  set.seed(11)
  expect_equal(simulate_raw_impacts(pop, supp_truth, n = 5000, joint = fit_joint(pop))$value, raw_n$value)
})

test_that("recovery holds in a denser random 6-disease system", {
  set.seed(42)
  ids <- paste0("x", 1:6)
  pr <- t(utils::combn(ids, 2))
  keep <- stats::runif(nrow(pr)) < 0.6
  pop <- cm_population(cm_diseases(ids, stats::runif(6, 0.03, 0.4)),
                       cm_associations(pr[keep, 1], pr[keep, 2], exp(stats::runif(sum(keep), -0.5, log(6)))))
  truth <- stats::setNames(stats::runif(6, 0.5, 8), ids)
  m <- cm_model(pop, simulate_raw_impacts(pop, truth))
  expect_equal(deconflate(m)$adjusted$adjusted, unname(truth), tolerance = 1e-7)
  adj <- simulate_raw_impacts(pop, truth, estimand = "adjusted_linear", adjusted_for = "x1; x2")
  expect_equal(deconflate(cm_model(pop, adj), method = "global")$adjusted$adjusted, unname(truth),
               tolerance = 1e-7)
})
