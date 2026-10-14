# The worked example of the Supplementary File, Rasmussen et al. (2022).
# Impacts are entered in percent, so results are in percent.
# Reference values computed independently in Python
# (inst/validation/python_reference.py).

test_that("the model solves the additive equations exactly", {
  res <- deconflate(example_supplement())
  expect_equal(res$adjusted$adjusted, c(2.143250478385, 3.387079589184, 6.934209451188),
               tolerance = 1e-9)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-12)
  # The pairwise shortcut gives the answer; the joint fit confirms that the
  # associations can hold together.
  expect_equal(res$diagnostics$method, "pairwise shortcut")
  expect_match(res$diagnostics$feasibility, "joint distribution fitted")
  expect_s3_class(res$joint, "cm_joint")
})

test_that("the full model equals the pairwise shortcut without interactions", {
  m <- example_supplement()
  full <- adjust_additive(m, joint = fit_joint(m))
  expect_equal(full$method, "full model")
  expect_equal(full$adjusted$adjusted, deconflate(m)$adjusted$adjusted, tolerance = 1e-8)
  expect_error(deconflate(m, method = "global"), "has been removed", class = "deconflate_unsupported")
  expect_error(deconflate(m, feasibility = "lp"), "has been removed", class = "deconflate_unsupported")
})

test_that("results are in the units supplied and scale with them", {
  kg <- supp_model(c(250, 500, 750), units = "kg/cow/year")
  res <- deconflate(kg)
  expect_equal(res$adjusted$adjusted, c(214.325047838512, 338.707958918446, 693.420945118756),
               tolerance = 1e-9)
  expect_equal(res$units, "kg/cow/year")
  a <- deconflate(supp_model(c(2.5, 5, 7.5)), n_draws = 0)$adjusted$adjusted
  b <- deconflate(supp_model(1000 * c(2.5, 5, 7.5)), n_draws = 0)$adjusted$adjusted
  expect_equal(b, 1000 * a, tolerance = 1e-10)
})

test_that("contributions add up to the aggregate", {
  res <- deconflate(example_supplement())
  ct <- res$contributions
  expect_equal(ct$main, c(0.10, 0.15, 0.20) * res$adjusted$adjusted)
  expect_equal(ct$interaction, c(0, 0, 0))
  expect_equal(sum(ct$total), res$totals$adjusted_total, tolerance = 1e-14)
  expect_equal(sum(ct$share), 1, tolerance = 1e-14)
  expect_equal(res$totals$raw_total, sum(c(0.10, 0.15, 0.20) * c(2.5, 5, 7.5)))
})
