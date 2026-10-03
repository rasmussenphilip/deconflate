# Reproduction of the Supplementary File, Rasmussen et al. (2022).
# Reference values computed independently in Python.

test_that("published method reproduces the de-conflated impacts", {
  res <- deconflate(example_supplement(), method = "published")
  # Published (rounded): 2.05, 3.70, 6.74 (%)
  expect_equal(res$adjusted$adjusted, c(0.02065001, 0.03699294, 0.06748473),
               tolerance = 1e-6)
  # The approximation does not reconstruct the raw impacts exactly.
  expect_gt(res$diagnostics$max_reconstruction_residual, 1e-4)
})

test_that("published productivity gap (unrounded) is reproduced", {
  res <- deconflate(example_supplement(), method = "published")
  pg <- productivity_gap(res, c(yield = 10000))
  expect_equal(pg$summary$disease_free, 10215.661683976388, tolerance = 1e-9)
  expect_equal(pg$attribution$gap, c(21.09535158, 56.68609684, 137.88023556),
               tolerance = 1e-7)
  expect_equal(sum(pg$attribution$gap), pg$summary$gap, tolerance = 1e-10)
})

test_that("the published 10,225 comes from rounding adjusted impacts first", {
  res <- deconflate(example_supplement(), method = "published")
  m <- round(res$adjusted$adjusted, 2)
  P <- c(0.10, 0.15, 0.20)
  xh <- 10000 / (1 - sum(m * P))
  expect_equal(xh, 10224.948875255624, tolerance = 1e-9)
  expect_equal(m * P / sum(m * P) * (xh - 10000),
               c(20.44989775, 61.34969325, 143.14928425), tolerance = 1e-7)
})

test_that("simultaneous method solves the additive equations exactly", {
  res <- deconflate(example_supplement())
  expect_equal(res$adjusted$adjusted, c(0.0214325, 0.0338708, 0.06934209),
               tolerance = 1e-6)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-12)
  pg <- productivity_gap(res, c(yield = 10000))
  expect_equal(pg$summary$disease_free, 10215.467592322026, tolerance = 1e-9)
  expect_equal(pg$attribution$gap, c(21.8943058, 51.90090266, 141.67238385),
               tolerance = 1e-7)
})

test_that("global method equals simultaneous without interactions", {
  m <- example_supplement()
  expect_equal(deconflate(m, method = "global")$adjusted$adjusted,
               deconflate(m)$adjusted$adjusted, tolerance = 1e-8)
})

test_that("method scale does not depend on percent vs proportion input", {
  m <- example_supplement()
  m2 <- m
  m2$impacts <- cm_impacts(c("d1", "d2", "d3"), c(0.025, 0.05, 0.075),
                           outcome = "yield", units = "units/animal/year")
  for (meth in c("published", "simultaneous")) {
    expect_equal(deconflate(m, method = meth)$adjusted$adjusted,
                 deconflate(m2, method = meth)$adjusted$adjusted)
  }
})
