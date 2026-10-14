# An invalid choice for an argument, and the problems found by
# cm_read_inputs(), are deconflate_invalid_input.

test_that("an invalid choice is deconflate_invalid_input and names the argument", {
  cls <- "deconflate_invalid_input"
  expect_error(deconflate(example_supplement(), sampling = "x"), "`sampling` must be one of: 'random', 'lhs'",
               class = cls)
  expect_error(deconflate(example_supplement(), pool = "x"), "`pool`", class = cls)
  expect_error(pair_tables(example_supplement(), pool = "x"), "`pool`", class = cls)
  expect_error(fit_joint(example_supplement(), backend = "x"), "`backend`", class = cls)
  expect_error(check_feasibility(example_supplement(), method = "x"), "`method`", class = cls)
  expect_error(example_global_dairy("x"), "`outcome`", class = cls)
  expect_error(example_uk_dairy_2022("x"), "`outcome`", class = cls)
  expect_error(cm_threshold(example_supplement(), "assoc:d2:d3", c(1, 3), conclusion = "x"), "`conclusion`",
               class = cls)
  expect_error(compare_methods(example_supplement(), methods = "x"),
               "`methods` must be one or more of: 'published', 'deconflate'", class = cls)
  expect_error(sensitivity_oat(example_supplement(), inputs = "x"), "`inputs` must be one or more of", class = cls)
})

test_that("valid choices are matched as by match.arg()", {
  expect_equal(cm_match_arg("l", c("random", "lhs")), "lhs")
  expect_equal(cm_match_arg(c("a", "b"), c("a", "b")), "a")
  expect_equal(cm_match_arg(NULL, c("a", "b")), "a")
  expect_equal(cm_match_arg(c("b", "a"), c("a", "b"), several_ok = TRUE), c("b", "a"))
  f <- function(x = c("one", "two")) cm_match_arg(x)
  expect_equal(f(), "one")
  expect_equal(f("tw"), "two")
  expect_error(f("three"), "`x` must be one of: 'one', 'two'", class = "deconflate_invalid_input")
  expect_equal(pair_tables(example_supplement(), pool = "fix"), pair_tables(example_supplement(), pool = "fixed"))
  expect_equal(compare_methods(example_supplement(), methods = "pub")$methods, "published")
})

test_that("problems found by cm_read_inputs() are also deconflate_invalid_input", {
  read <- function() cm_read_inputs(diseases = data.frame(disease = "a", value = 2),
                                    impacts = data.frame(disease = "a", value = 1))
  expect_error(read(), class = "deconflate_input_problems")
  expect_error(read(), class = "deconflate_invalid_input")
})
