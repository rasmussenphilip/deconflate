# Validation against simulated populations with known effects.

supp_true <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)

test_that("simulated raw impacts match the Python reference", {
  raw <- simulate_raw_impacts(example_supplement(), supp_true)
  expect_equal(raw$value, c(0.02421306, 0.05406438, 0.06668175), tolerance = 1e-7)
})

test_that("simultaneous method recovers true additive impacts; published does not", {
  m <- example_supplement()
  m$impacts <- simulate_raw_impacts(m, supp_true, outcome = "yield")
  sim <- deconflate(m)
  expect_equal(sim$adjusted$adjusted, unname(supp_true), tolerance = 1e-8)
  pub <- deconflate(m, method = "published")
  expect_gt(max(abs(pub$adjusted$adjusted - supp_true)), 1e-3)
})

test_that("global method recovers true impacts with interactions", {
  m <- example_supplement()
  ints <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(0.01, 0.015),
                          outcome = "yield")
  raw <- simulate_raw_impacts(m, supp_true, interactions = ints, outcome = "yield")
  expect_equal(raw$value, c(0.02714619, 0.06136904, 0.07116677), tolerance = 1e-7)

  m_int <- cm_model(m$diseases, m$associations, raw, interactions = ints)
  res <- deconflate(m_int, method = "global")
  expect_equal(res$adjusted$adjusted, unname(supp_true), tolerance = 1e-8)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-10)
  expect_error(deconflate(m_int), "Interactions require")
})

test_that("recovery holds in a denser random 6-disease system", {
  set.seed(42)
  ids <- paste0("x", 1:6)
  pr <- t(utils::combn(ids, 2))
  keep <- stats::runif(nrow(pr)) < 0.6
  m <- cm_model(
    cm_diseases(ids, stats::runif(6, 0.03, 0.4)),
    cm_associations(pr[keep, 1], pr[keep, 2], exp(stats::runif(sum(keep), -0.5, log(6))))
  )
  truth <- stats::setNames(stats::runif(6, 0.005, 0.08), ids)
  m$impacts <- simulate_raw_impacts(m, truth)
  expect_equal(deconflate(m)$adjusted$adjusted, unname(truth), tolerance = 1e-7)
})

test_that("Shapley shares sum to the total and split interactions equally", {
  m <- example_supplement()
  ints <- cm_interactions("d1", "d2", 0.01, outcome = "yield")
  m_int <- cm_model(m$diseases, m$associations,
                    simulate_raw_impacts(m, supp_true, ints, outcome = "yield"),
                    interactions = ints)
  res <- deconflate(m_int, method = "global")
  b <- attribute_burden(res)
  P12 <- res$joint_pairs["d1", "d2"]
  expect_equal(b$interaction, c(0.005 * P12, 0.005 * P12, 0), tolerance = 1e-12)
  L <- sum(supp_true * c(0.10, 0.15, 0.20)) + 0.01 * P12
  expect_equal(sum(b$total), L, tolerance = 1e-8)
  expect_equal(sum(b$share), 1)
})
