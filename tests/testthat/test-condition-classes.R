# Every error a user can meet has a class beside deconflate_error (see
# ?deconflate-package, "Conventions").

test_that("invalid arguments, tables and values are deconflate_invalid_input", {
  cls <- "deconflate_invalid_input"
  expect_error(deconflate(example_supplement(), n_draws = -1), class = cls)
  expect_error(deconflate(example_supplement(), event_model = NA), class = cls)
  expect_error(cm_diseases(c("a", "a"), c(0.1, 0.2)), class = cls)
  expect_error(cm_diseases(c("a", "no b"), c(0.1, 0.2)), class = cls)
  expect_error(cm_impacts(ids3, c(1, 2)), class = cls)
  expect_error(dist_normal(0, -1), class = cls)
  expect_error(dist_pert(2, 1, 0), class = cls)
  expect_error(cm_read_inputs(impacts = data.frame(disease = "a", value = 1)), class = cls)
  expect_error(cm_threshold(example_supplement(), "assoc:d2:d3", c(3, 1)), class = cls)
  expect_error(cm_threshold(example_supplement(), "nothing:d2", c(1, 3)), class = cls)
  e <- tryCatch(dist_beta(-1, 1), deconflate_error = function(e) e)
  expect_s3_class(e, c("deconflate_invalid_input", "deconflate_error", "error", "condition"), exact = TRUE)
})

test_that("limits of the package are deconflate_unsupported", {
  expect_error(fit_joint(example_supplement(), max_diseases = 2), class = "deconflate_unsupported")
  expect_error(cm_threshold(example_supplement(), "risk", c(0.1, 0.5)), class = "deconflate_unsupported")
})

test_that("condition labels include invalid_input", {
  e <- tryCatch(dist_normal(0, -1), deconflate_error = function(e) e)
  expect_equal(condition_type(e), "invalid_input")
  expect_equal(type_class("invalid_input"), "deconflate_invalid_input")
  expect_equal(type_class("singular"), "deconflate_singular")
  expect_null(type_class("error"))
  expect_null(type_class("failed"))
})

test_that("errors about a failed baseline or run are of the kind that caused them", {
  expect_error(check_baseline(structure(list(), reason = "singular")), "singular", class = "deconflate_singular")
  e <- tryCatch(check_baseline(structure(list(), reason = "error")), deconflate_error = function(e) e)
  expect_s3_class(e, c("deconflate_error", "error", "condition"), exact = TRUE)
  bad <- cm_model(bad_population(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(compare_methods(bad), "All methods failed", class = "deconflate_infeasible")
  zero <- supp_model(c(0, 0, 0))
  expect_error(cm_threshold(zero, "assoc:d2:d3", c(1, 5), conclusion = "change", target = 0.1), "zero",
               class = "deconflate_nonfinite")
})
