# Rasmussen et al. (2022): culling pipeline and Tables 8-10.
# Reference values computed independently in Python.

ids <- c("CO", "DA", "DYS", "FAS", "GIN", "LAM", "MAS", "MET", "MF", "NEO", "PTB", "RP", "SCK")

test_that("hazard ratios convert to the published excess culling risks", {
  m <- example_uk_dairy_2022()
  conv <- attr(m, "culling_conversion")
  expect_equal(conv$disease, ids)
  expect_equal(conv$excess,
               c(0, 0.3138414572, 0.1420512995, 0, 0, 0.2556484949, 0.2130688027,
                 0.173857642, 0.2056704258, 0.0988932697, 0.1963730607, 0, 0.1572642738),
               tolerance = 1e-8)
  cull <- m$impacts[m$impacts$outcome == "culling", ]
  expect_equal(cull$value, conv$excess)
})

test_that("published method reproduces the computed 2022 results", {
  m <- example_uk_dairy_2022()
  eco <- uk_dairy_2022_economics()
  res <- deconflate(m, method = "published")
  pg <- productivity_gap(res, eco$observed)
  s <- pg$summary
  expect_equal(s$disease_free[s$outcome == "yield"], 9299.207239070794, tolerance = 1e-9)
  expect_equal(s$disease_free[s$outcome == "fertility"], 375.0932260546869, tolerance = 1e-9)
  expect_equal(s$disease_free[s$outcome == "culling"], 22.545682567951776, tolerance = 1e-9)
  vl <- value_losses(pg, eco$unit_value, eco$additional)
  expect_equal(vl$total, 402.24755302935023, tolerance = 1e-9)
})

test_that("fertility and culling match the published tables", {
  m <- example_uk_dairy_2022()
  eco <- uk_dairy_2022_economics()
  res <- deconflate(m, method = "published")
  pg <- productivity_gap(res, eco$observed)
  vl <- value_losses(pg, eco$unit_value)
  fert <- vl$by_outcome[vl$by_outcome$outcome == "fertility", ]
  expect_equal(pg$summary$disease_free[pg$summary$outcome == "fertility"], 375.09, tolerance = 0.01 / 375)
  expect_equal(fert$value, 101.79, tolerance = 0.02 / 101.79)
  # Table 8 adjusted culling hazard ratios (2 dp).
  hr <- adjusted_hr(res, attr(m, "culling_conversion"), method = "published")
  table8 <- c(CO = 1, DA = 2.68, DYS = 1.60, FAS = 1, GIN = 1, LAM = 3.00, MAS = 2.22,
              MET = 1.55, MF = 1.89, NEO = 1.60, PTB = 1.64, RP = 1, SCK = 1.29)
  expect_lt(max(abs(hr$hr_adjusted - table8[hr$disease])), 0.015)
  expect_equal(hr$hr_adjusted,
               c(1, 2.6714376893, 1.5911278741, 1, 1, 2.9992762358, 2.2076274127,
                 1.5514896989, 1.8904358108, 1.6, 1.6422531356, 1, 1.2906637344),
               tolerance = 1e-8)
  # Culling gap within 0.5% of the published 4.44 percentage points.
  cg <- pg$summary$gap[pg$summary$outcome == "culling"]
  expect_lt(abs(cg - 4.44) / 4.44, 0.005)
})

test_that("SCK yield of 340 kg / 8737 kg reproduces the Table 8 adjusted value", {
  m <- example_uk_dairy_2022(yield_sck = 100 * 340 / 8737)
  res <- deconflate(m, method = "published")
  y <- res$adjusted[res$adjusted$outcome == "yield", ]
  expect_lt(abs(y$adjusted[y$disease == "SCK"] - 0.0265), 0.0002)
  expect_lt(abs(y$adjusted[y$disease == "LAM"] - 0.0476), 0.0002)
})

test_that("simultaneous method on the 2022 inputs matches the reference", {
  m <- example_uk_dairy_2022()
  pg <- productivity_gap(suppressWarnings(deconflate(m)), uk_dairy_2022_economics()$observed)
  expect_equal(pg$summary$disease_free,
               c(9292.915419609979, 376.70813848633463, 22.61083366279789),
               tolerance = 1e-9)
})

test_that("absolute culling scale subtracts excess risk from the observed rate", {
  m <- example_uk_dairy_2022(culling_method = "proportional_hazards", culling_scale = "absolute")
  res <- suppressWarnings(deconflate(m))
  a <- res$adjusted[res$adjusted$outcome == "culling", ]
  P <- m$diseases$prob[match(a$disease, m$diseases$id)]
  pg <- productivity_gap(res, uk_dairy_2022_economics("absolute")$observed)
  expect_equal(pg$summary$disease_free[pg$summary$outcome == "culling"],
               0.27 - sum(a$adjusted * P), tolerance = 1e-12)
})

test_that("proportional-hazards back-conversion recovers the input hazard ratios", {
  d <- cm_diseases(c("a", "b"), c(0.2, 0.1))
  conv <- hr_conversion(d, c(a = 2.5, b = 1.7), overall_risk = 0.25)
  m <- cm_model(d, impacts = as_impacts(conv))   # no associations: no adjustment
  hr <- adjusted_hr(deconflate(m), conv)
  expect_equal(hr$hr_adjusted, c(2.5, 1.7), tolerance = 1e-8)
})

test_that("combine_impacts validates the stacked table", {
  a <- cm_impacts(c("x", "y"), c(1, 2), outcome = "o1", scale = "percent")
  b <- cm_impacts(c("x", "y"), c(3, 4), outcome = "o2")
  ab <- combine_impacts(a, b)
  expect_s3_class(ab, "cm_impacts")
  expect_equal(nrow(ab), 4)
  expect_error(combine_impacts(a, a), "only one impact")
})
