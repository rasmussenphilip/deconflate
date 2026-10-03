test_that("input constructors validate their arguments", {
  expect_error(cm_diseases(c("a", "a"), c(0.1, 0.2)), "unique")
  expect_error(cm_diseases("a", 1.2), "strictly between")
  expect_error(cm_associations("a", "a", 2), "different")
  expect_error(cm_associations(c("a", "b"), c("b", "a"), c(2, 3)), "Duplicate")
  expect_error(cm_associations("a", "b", -1), "positive")
  expect_error(cm_impacts("a", 1, scale = "absolute"), "units")
  expect_warning(cm_diseases(c("a", "b"), c(0.1, 0.2), time_horizon = c("year", "lactation")),
                 class = "deconflate_warning")
  d <- cm_diseases(c("a", "b"), c(0.1, 0.2))
  expect_error(cm_model(d, impacts = cm_impacts("a", 0.01)), "no impact for: b")
})

test_that("adjusted_for removes the corresponding conflation terms", {
  m <- example_supplement()
  m$impacts <- cm_impacts(c("d1", "d2", "d3"), c(0.025, 0.05, 0.075), outcome = "yield",
                          adjusted_for = c("d2", NA, NA))
  res <- deconflate(m)
  expect_equal(unname(res$conflation$yield$A["d1", "d2"]), 0)
  expect_equal(unname(res$conflation$yield$A["d2", "d1"]),
               unname(excess_matrix(m)["d1", "d2"]))
})

test_that("Monte Carlo propagates draws, rejects infeasible ones and reweights", {
  base <- example_supplement()
  sampler <- function(i) {
    m <- base
    # Draw 3 uses an infeasible conditional probability, P(d1 | d2) = 0.9.
    if (i == 3) {
      m$associations <- cm_associations("d1", "d2", 0.9, measure = "cond_prob")
    } else {
      m$associations <- cm_associations(c("d1", "d2"), c("d2", "d3"),
                                        exp(stats::rnorm(2, log(c(2, 3)), 0.2)))
    }
    m
  }
  set.seed(1)
  mc <- cm_monte_carlo(sampler, 20)
  expect_equal(mc$n_rejected, 1)
  expect_equal(nrow(mc$params), 19)
  expect_equal(length(unique(mc$draws$draw)), 19)
  s <- summary(mc)
  expect_equal(nrow(s), 3)
  expect_true(all(s$`q0.025` <= s$mean & s$mean <= s$`q0.975`))

  # Uniform reweighting leaves the summary unchanged and gives ESS = n.
  mc2 <- cm_reweight(mc, function(p) rep(0, nrow(p)))
  expect_equal(mc2$ess, 19)
  expect_equal(summary(mc2)$mean, s$mean)

  # Shifting the d1:d2 odds ratio upwards lowers d1's adjusted impact.
  mc3 <- suppressWarnings(cm_reweight(mc, function(p) 3 * log(p[["assoc:d1:d2"]])))
  d1 <- function(x) summary(x)$mean[summary(x)$disease == "d1"]
  expect_lt(d1(mc3), d1(mc))
})
