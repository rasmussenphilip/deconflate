# Errors and warnings hold the inputs they are about as fields (see
# ?deconflate, "Errors and warnings").

catch <- function(expr, class = "deconflate_error") {
  tryCatch(expr, condition = function(e) if (inherits(e, class)) e else stop(e))
}

test_that("cm_read_inputs() errors hold every problem found", {
  f <- function(x) extdata("example_with_errors", x)
  e <- catch(cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                            impacts = f("impacts_yield.csv")), "deconflate_input_problems")
  chk <- cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                         impacts = f("impacts_yield.csv"))
  expect_s3_class(e$problems, "cm_problems")
  expect_equal(e$problems, chk)
  expect_true(any(e$problems$severity == "error"))
})

test_that("infeasible pair rows give the conflicting triples", {
  m <- cm_model(bad_population(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  e <- catch(deconflate(m), "deconflate_infeasible")
  expect_equal(e$triples[, c("disease1", "disease2", "disease3")],
               data.frame(disease1 = "a", disease2 = "b", disease3 = "c", stringsAsFactors = FALSE))
  expect_gt(e$triples$gap, 0)
  expect_equal(catch(adjust_published(m), "deconflate_infeasible")$triples, e$triples)
})

test_that("calculations that need every pair name the missing pairs", {
  e <- catch(excess_matrix(supp_unknown_population()), "deconflate_unknown_pairs")
  expect_equal(e$pairs, data.frame(disease1 = "d1", disease2 = "d3", stringsAsFactors = FALSE))
  m <- cm_model(supp_unknown_population(), cm_impacts(ids3, c(1, 2, 3)))
  expect_equal(catch(adjust_published(m), "deconflate_unknown_pairs")$pairs, e$pairs)
})

test_that("sign changes and unidentified terms are named", {
  m <- cm_model(bad_population(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  w <- catch(adjust_additive(m, pt = pair_tables(m)), "deconflate_sign_change")
  expect_s3_class(w, "deconflate_warning")
  expect_equal(w$terms, c("b", "c"))
  e <- catch(identify_terms(matrix(1, 2, 2), c("a", "b")), "deconflate_singular")
  expect_equal(e$terms, c("a", "b"))
})

test_that("contradicting exact rows are given by table and row", {
  a <- supp_population(cm_associations("d1", "d2", 3))
  e <- catch(deconflate(cm_model(a, cm_impacts(ids3, c(1, 2, 3)))), "deconflate_infeasible")
  expect_equal(e$table, "associations")
  expect_equal(e$rows, c(1, 4))
  im <- cm_impacts(c("d1", "d1", "d2", "d3"), c(2, 3, 5, 7.5))
  e <- catch(deconflate(cm_model(supp_population(), im)), "deconflate_infeasible")
  expect_equal(e$table, "impacts")
  # Every row is exact, so all four are reported, largest residual first.
  expect_setequal(e$rows, 1:4)
  expect_setequal(e$rows[1:2], 1:2)
  expect_equal(abs(e$residuals), sort(abs(e$residuals), decreasing = TRUE))
})

test_that("errors about one row give its table and row", {
  e <- catch(supp_population(cm_associations("d1", "d2", 0.1, measure = "RR", adjusted_for = "d3")),
             "deconflate_unsupported")
  expect_equal(e$table, "associations")
  expect_equal(e$row, 4)
})

test_that("non-convergence gives the residual", {
  w <- catch(fit_joint(cm_model(bad_population(), cm_impacts(c("a", "b", "c"), c(1, 2, 3))), max_iter = 50),
             "deconflate_nonconvergence")
  expect_s3_class(w, "deconflate_warning")
  expect_true(is.numeric(w$residual) && w$residual > 0)
})

test_that("an incomplete Shapley allocation gives the amounts", {
  j <- fit_joint(example_supplement())
  w <- catch(shapley_by_cell(j, loss_additive(c(d1 = 1, d2 = 2, d3 = 3)), max_present = 2),
             "deconflate_incomplete_allocation")
  sh <- suppressWarnings(shapley_by_cell(j, loss_additive(c(d1 = 1, d2 = 2, d3 = 3)), max_present = 2))
  expect_equal(w$skipped_mass, attr(sh, "skipped_mass"))
  expect_equal(w$skipped_loss, attr(sh, "skipped_loss"))
  expect_equal(w$total, attr(sh, "total"))
})
