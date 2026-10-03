# Rasmussen et al. (2022) Table 6 (culling rate 0.27). Reference values from Python.

test_that("or_approx reproduces the published excess culling probabilities", {
  r <- hr_to_risk(c(3.83, 3.40), c(0.03, 0.30), 0.27, method = "or_approx")
  expect_equal(r$risk_exposed, c(0.5744262135258521, 0.4489539464210607), tolerance = 1e-9)
  expect_equal(r$risk_unexposed, c(0.2605847562827056, 0.19330545153383116), tolerance = 1e-9)
  expect_equal(round(r$excess, 2), c(0.31, 0.26))  # Table 6: 0.31, 0.25 (rounded)
})

test_that("proportional-hazards conversion is consistent", {
  r <- hr_to_risk(c(3.83, 3.40), c(0.03, 0.30), 0.27)
  expect_equal(r$risk_exposed, c(0.6799848806532945, 0.48561141072390346), tolerance = 1e-8)
  expect_equal(r$risk_unexposed, c(0.257320055237527, 0.17759510968975575), tolerance = 1e-8)
  # Population risk is preserved.
  expect_equal((1 - r$prevalence) * r$risk_unexposed + r$prevalence * r$risk_exposed,
               c(0.27, 0.27), tolerance = 1e-10)
  # And the hazard ratio is recovered.
  expect_equal(excess_to_hr(r$excess, r$risk_unexposed), c(3.83, 3.40), tolerance = 1e-8)
})
