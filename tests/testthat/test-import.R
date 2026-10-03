# Reading inputs from CSV files and data frames.

test_that("the template round-trips through cm_read_inputs()", {
  dir <- file.path(tempdir(), paste0("deconflate-template-", sample.int(1e6, 1)))
  expect_message(paths <- cm_template(dir), "Wrote 5 files")
  expect_true(all(file.exists(paths)))
  expect_error(cm_template(dir), "already exist")
  inp <- cm_read_inputs(dir = dir)
  expect_s3_class(inp, "cm_inputs")
  expect_equal(inp$model$diseases$id, c("LAM", "SCK", "MET"))
  expect_equal(inp$model$diseases$prob[2], 1 - exp(-0.48))
  expect_equal(unique(inp$model$impacts$scale), c("proportion", "hazard_ratio"))
  expect_null(inp$model$interactions)
  expect_equal(names(attr(inp$sampler, "specs")),
               c("prob:LAM", "assoc:LAM:SCK", "impact:yield:LAM", "impact:yield:SCK",
                 "impact:culling:LAM"))
  expect_output(print(inp), "Uncertain inputs: 5")
  res <- deconflate(inp$model, method = "global")
  expect_equal(unique(res$adjusted$outcome), c("yield", "culling"))
  mc <- cm_monte_carlo(inp$sampler, 20, seed = 1)
  expect_equal(mc$n_rejected, 0)
  unlink(dir, recursive = TRUE)
})

test_that("tables typed in R build the same model as the constructors", {
  inp <- cm_read_inputs(
    diseases = data.frame(id = c("d1", "d2", "d3"), value = c(0.10, 0.15, 0.20)),
    associations = data.frame(disease1 = c("d1", "d1", "d2"), disease2 = c("d2", "d3", "d3"),
                              value = c(2, 1, 3)),
    impacts = data.frame(disease = c("d1", "d2", "d3"), outcome = "yield",
                         value = c(2.5, 5, 7.5), scale = "percent")
  )
  expect_equal(deconflate(inp$model)$adjusted$adjusted,
               deconflate(example_supplement())$adjusted$adjusted)
  expect_null(inp$sampler)
})

test_that("all problems are reported at once", {
  pr <- cm_check_inputs(
    diseases = data.frame(id = c("a", "b", "b"), value = c(0.1, 1.5, 0.2)),
    associations = data.frame(disease1 = "a", disease2 = "c", value = -1),
    impacts = data.frame(disease = c("a", "b"), outcome = "yield", value = c("2", "x"),
                         scale = "percent", extra = 1),
    uncertainty = data.frame(key = c("impact:yield:a", "prob:z", "foo"),
                             dist = c("normal", "beta", "pert"), p1 = c(2, 1, 1))
  )
  expect_s3_class(pr, "cm_problems")
  err <- pr$problem[pr$severity == "error"]
  expect_true(any(grepl("Duplicate disease id 'b'", err)))
  expect_true(any(grepl("between 0 and 1", err)))
  expect_true(any(grepl("Unknown disease 'c'", err)))
  expect_true(any(grepl("must be positive", err)))
  expect_true(any(grepl("'x' is not a number", err)))
  expect_true(any(grepl("needs p1 \\(mean\\) and p2 \\(sd\\)", err)))
  expect_true(any(grepl("does not match an input", err)))
  expect_true(any(grepl("must start with", err)))
  expect_true(any(pr$severity == "note" & pr$column %in% "extra"))
  expect_output(print(pr), "problem")
  expect_error(cm_read_inputs(diseases = data.frame(id = "a", value = 2)), "Found 1 problem")
})

test_that("missing impacts and bad files are reported", {
  pr <- cm_check_inputs(
    diseases = data.frame(id = c("a", "b"), value = c(0.1, 0.2)),
    impacts = data.frame(disease = "a", outcome = "yield", value = 0.02)
  )
  expect_true(any(grepl("Outcome 'yield' has no impact for: b", pr$problem)))
  pr2 <- cm_check_inputs(diseases = file.path(tempdir(), "no-such-file.csv"))
  expect_true(any(grepl("File not found", pr2$problem)))
})

test_that("cm_dist_table builds distributions", {
  d <- cm_dist_table(data.frame(key = c("a", "b", "c"), dist = c("pert", "lognormal_ci", "normal"),
                                p1 = c(1, 2, 2.6), p2 = c(2, 1.4, 1.4), p3 = c(4, 2.9, 0)))
  expect_equal(d$a$type, "pert")
  expect_equal(d$a$params$mode, 2)
  expect_equal(d$c$params$lower, 0)
  expect_error(cm_dist_table(data.frame(key = "a", dist = "gamma", p1 = 1)), "Unknown distribution")
})
