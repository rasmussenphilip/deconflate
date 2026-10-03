# Reproduction of Rasmussen et al. (2022), Tables 8-10, and (2024), Table 5.
# Reference values computed independently in Python
# (inst/validation/reference_2022_shapley.py, python_reference.py).

ids <- c("CO", "DA", "DYS", "FAS", "GIN", "LAM", "MAS", "MET", "MF", "NEO", "PTB", "RP", "SCK")

test_that("hazard ratios treated as odds ratios give the 2022 excess culling risks", {
  uk <- example_uk_dairy_2022()
  cu <- uk$models$culling$impacts
  expect_equal(cu$disease, ids)
  expect_equal(cu$value,
               c(0, 0.3138414572, 0.1420512995, 0, 0, 0.2556484949, 0.2130688027,
                 0.173857642, 0.2056704258, 0.0988932697, 0.1963730607, 0, 0.1572642738),
               tolerance = 1e-8)
  expect_s3_class(attr(uk, "hazard_ratios"), "cm_hazard_ratios")
})

test_that("reproduce_rasmussen_2022() recomputes the published gaps and values", {
  r <- reproduce_rasmussen_2022()
  expect_s3_class(r, "cm_reproduction")
  expect_equal(r$gaps$analysis, c("yield", "fertility", "culling"))
  expect_equal(r$gaps$disease_free, c(9299.207239070794, 375.0932260546869, 22.545682567951776),
               tolerance = 1e-9)
  expect_equal(r$total, 402.24755302935023, tolerance = 1e-9)
  expect_equal(r$additional, c(veterinary = 71.09))
  # Fertility matches Table 9.
  expect_equal(r$gaps$disease_free[2], 375.09, tolerance = 0.01 / 375)
  expect_equal(r$gaps$value[2], 101.79, tolerance = 0.02 / 101.79)
  # Culling gap within 0.5% of the published 4.44 percentage points.
  expect_lt(abs(r$gaps$gap[3] - 4.44) / 4.44, 0.005)
  # Per-disease values add up to the analysis values.
  expect_equal(sum(r$values$total), sum(r$gaps$value), tolerance = 1e-10)
  expect_output(print(r), "Table 9")
})

test_that("adjusted culling hazard ratios (eq. 23) match Table 8", {
  r <- reproduce_rasmussen_2022()
  hc <- r$hr_comparison
  expect_equal(hc$disease, ids)
  expect_lt(max(abs(hc$package - hc$table8)), 0.015)
  expect_equal(hc$package,
               c(1, 2.6714376893, 1.5911278741, 1, 1, 2.9992762358, 2.2076274127,
                 1.5514896989, 1.8904358108, 1.6, 1.6422531356, 1, 1.2906637344),
               tolerance = 1e-8)
})

test_that("SCK yield of 340 kg / 8737 kg reproduces Table 8 adjusted values", {
  r <- reproduce_rasmussen_2022(yield_sck = 100 * 340 / 8737)
  a <- r$adjusted
  expect_lt(abs(a$yield[a$disease == "SCK"] - 2.65), 0.02)
  expect_lt(abs(a$yield[a$disease == "LAM"] - 4.76), 0.02)
})

test_that("the simultaneous method on the 2022 inputs matches the reference", {
  uk <- example_uk_dairy_2022()
  eco <- uk_dairy_2022_economics()
  res <- suppressWarnings(deconflate(uk))
  df <- vapply(names(res), function(nm) {
    v <- eco$valuation[[nm]]
    productivity_gap(res[[nm]], v$observed, v$direction, v$effect)$summary$disease_free
  }, numeric(1))
  expect_equal(unname(df), c(9292.915419609979, 376.70813848633463, 22.61083366279789),
               tolerance = 1e-9)
  val <- sum(vapply(names(res), function(nm) {
    contribution_table(res[[nm]], eco$valuation[[nm]])$value
  }, numeric(length(ids)))) + eco$additional
  expect_equal(unname(val), 393.1318185491835, tolerance = 1e-9)
})

test_that("reproduce_rasmussen_2024() runs the published Monte Carlo", {
  skip_on_cran()
  r <- reproduce_rasmussen_2024(n_draws = 30, seed = 1)
  expect_s3_class(r, "cm_reproduction")
  cmp <- r$comparison
  expect_equal(unique(cmp$analysis), c("yield", "fertility", "culling_hr"))
  expect_equal(nrow(cmp), 36)
  expect_true(all(is.finite(cmp$mean)))
  expect_true(all(cmp$stability %in% c("ok", "imprecise", "heavy_tail", "possible_pole")))
  expect_output(print(r), "Table 5")
})
