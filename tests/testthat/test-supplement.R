# Reproduction of the Supplementary File, Rasmussen et al. (2022).
# Impacts are entered in percent, so results are in percent.
# Reference values computed independently in Python
# (inst/validation/python_reference.py).

test_that("published method reproduces the de-conflated impacts", {
  res <- deconflate(example_supplement(), method = "published")
  # Published (rounded): 2.05, 3.70, 6.74 (%)
  expect_equal(res$adjusted$adjusted, c(2.065001, 3.699294, 6.748473), tolerance = 1e-6)
  # The approximation does not reconstruct the raw impacts exactly.
  expect_gt(res$diagnostics$max_reconstruction_residual, 1e-2)
  expect_equal(res$units, "%")
  expect_equal(res$label, "yield")
})

test_that("published productivity gap (unrounded) is reproduced", {
  res <- deconflate(example_supplement(), method = "published")
  pg <- productivity_gap(res, 10000, direction = "decrease", effect = "percent")
  expect_equal(pg$summary$disease_free, 10215.661683976388, tolerance = 1e-9)
  expect_equal(pg$attribution$gap, c(21.09535158, 56.68609684, 137.88023556),
               tolerance = 1e-7)
  expect_equal(sum(pg$attribution$gap), pg$summary$gap, tolerance = 1e-10)
})

test_that("the published 10,225 comes from rounding adjusted impacts first", {
  res <- deconflate(example_supplement(), method = "published")
  m <- round(res$adjusted$adjusted / 100, 2)
  P <- c(0.10, 0.15, 0.20)
  xh <- 10000 / (1 - sum(m * P))
  expect_equal(xh, 10224.948875255624, tolerance = 1e-9)
  expect_equal(m * P / sum(m * P) * (xh - 10000),
               c(20.44989775, 61.34969325, 143.14928425), tolerance = 1e-7)
})

test_that("simultaneous method solves the additive equations exactly", {
  res <- deconflate(example_supplement())
  expect_equal(res$adjusted$adjusted, c(2.143250478385, 3.387079589184, 6.934209451188),
               tolerance = 1e-9)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-12)
  expect_equal(res$diagnostics$feasibility, "triple screen passed (necessary condition only)")
  pg <- productivity_gap(res, 10000, "decrease", "percent")
  expect_equal(pg$summary$disease_free, 10215.467592322026, tolerance = 1e-9)
  expect_equal(pg$attribution$gap, c(21.8943058, 51.90090266, 141.67238385),
               tolerance = 1e-7)
})

test_that("global method equals simultaneous without interactions", {
  m <- example_supplement()
  expect_equal(deconflate(m, method = "global")$adjusted$adjusted,
               deconflate(m)$adjusted$adjusted, tolerance = 1e-8)
})

test_that("results are in the units supplied and scale with them", {
  kg <- supp_model(c(250, 500, 750), units = "kg/cow/year")
  res <- deconflate(kg)
  expect_equal(res$adjusted$adjusted, c(214.325047838512, 338.707958918446, 693.420945118756),
               tolerance = 1e-9)
  expect_equal(res$units, "kg/cow/year")
  for (meth in c("published", "simultaneous", "global")) {
    a <- deconflate(supp_model(c(2.5, 5, 7.5)), method = meth)$adjusted$adjusted
    b <- deconflate(supp_model(1000 * c(2.5, 5, 7.5)), method = meth)$adjusted$adjusted
    expect_equal(b, 1000 * a, tolerance = 1e-10)
  }
})

test_that("contributions add up to the aggregate", {
  res <- deconflate(example_supplement())
  ct <- attribute_burden(res)
  expect_equal(ct$main, c(0.10, 0.15, 0.20) * res$adjusted$adjusted)
  expect_equal(ct$interaction, c(0, 0, 0))
  expect_equal(sum(ct$total), res$totals$adjusted_total, tolerance = 1e-14)
  expect_equal(sum(ct$share), 1, tolerance = 1e-14)
  expect_equal(res$totals$raw_sum, sum(c(0.10, 0.15, 0.20) * c(2.5, 5, 7.5)))
})
