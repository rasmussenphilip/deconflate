# Reference values computed independently in Python (scipy/numpy).

test_that("odds ratio quadratic reproduces the Supplementary File", {
  p11 <- or_to_joint(2, 0.10, 0.15)
  expect_equal(p11 / 0.15, 0.1631959501785622, tolerance = 1e-10)  # P(1 | 2)
  expect_equal(p11 / 0.10, 0.24479392526784327, tolerance = 1e-10) # P(2 | 1)
})

test_that("odds ratio conversion round-trips, including OR < 1 and OR near 1", {
  grid <- expand.grid(or = c(0.2, 0.4, 0.999999, 1 + 1e-9, 1.5, 3, 9.7),
                      p1 = c(0.02, 0.3, 0.6), p2 = c(0.05, 0.38, 0.8))
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    p11 <- or_to_joint(g$or, g$p1, g$p2)
    expect_true(p11 >= max(0, g$p1 + g$p2 - 1) && p11 <= min(g$p1, g$p2))
    expect_equal(joint_to_or(p11, g$p1, g$p2), g$or, tolerance = 1e-6)
  }
})

test_that("all association measures map to the same 2x2 table", {
  p1 <- 0.12; p2 <- 0.3
  p11 <- or_to_joint(2.5, p1, p2)
  pd1_d2 <- p11 / p2
  pd1_nd2 <- (p1 - p11) / (1 - p2)
  phi <- (p11 - p1 * p2) / sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  expect_equal(association_to_joint("RR", pd1_d2 / pd1_nd2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("RD", pd1_d2 - pd1_nd2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("cond_prob", pd1_d2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("phi", phi, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("independent", NA, p1, p2), p1 * p2)
})

test_that("incompatible measures are flagged as infeasible", {
  expect_error(association_to_joint("cond_prob", 0.9, 0.1, 0.5),
               class = "deconflate_infeasible")
  expect_error(association_to_joint("RD", 0.9, 0.1, 0.5),
               class = "deconflate_infeasible")
})

test_that("contingency tables are converted via their odds ratio", {
  a <- cm_associations("a", "b", measure = "table", n11 = 20, n10 = 30, n01 = 40, n00 = 110)
  expect_equal(a$value, (20 * 110) / (30 * 40))
  a0 <- cm_associations("a", "b", measure = "table", n11 = 0, n10 = 30, n01 = 40, n00 = 110)
  expect_equal(a0$value, (0.5 * 110.5) / (30.5 * 40.5))
})

test_that("excess matrix matches the Supplementary File", {
  E <- excess_matrix(example_supplement())
  # E[k, i] = P(k | i) - P(k | not i)
  expect_equal(E["d2", "d1"], 0.10532658, tolerance = 1e-7)  # ep_2,1 = 0.11
  expect_equal(E["d1", "d2"], 0.07434818, tolerance = 1e-7)  # ep_1,2 = 0.07
  expect_equal(E["d3", "d2"], 0.20962356, tolerance = 1e-7)  # ep_3,2 = 0.21
  expect_equal(E["d2", "d3"], 0.16704377, tolerance = 1e-7)  # ep_2,3 = 0.17
  expect_equal(E["d1", "d3"], 0, tolerance = 1e-12)
})

test_that("unknown pairs are rejected by pairwise methods but allowed globally", {
  m <- cm_model(
    cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25)),
    cm_impacts(c("a", "b", "c"), c(3, 4, 5)),
    associations = cm_associations(c("a", "b"), c("b", "c"), c(3, 2)),
    missing_associations = "unknown"
  )
  expect_error(deconflate(m), class = "deconflate_unknown_pairs")
  res <- deconflate(m, method = "global")
  expect_s3_class(res, "cm_result")
  # Unconstrained a:c is implied by the maximum-entropy fit, so it differs
  # from imposing independence.
  ind <- cm_model(m$diseases, m$impacts, associations = m$associations,
                  missing_associations = "independent")
  expect_false(isTRUE(all.equal(res$adjusted$adjusted,
                                deconflate(ind)$adjusted$adjusted)))
})
