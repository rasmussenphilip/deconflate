# Reading inputs from CSV files and data frames (cm_read_inputs() and friends).

# The Supplementary File population as tables typed in R.
supp_tables <- function() {
  list(
    diseases = data.frame(id = ids3, value = c(0.10, 0.15, 0.20)),
    associations = data.frame(disease1 = c("d1", "d1", "d2"), disease2 = c("d2", "d3", "d3"),
                              value = c(2, 1, 3))
  )
}

# Two impact tables (analyses yield and fertility) on the supplement diseases.
two_impacts <- function() {
  list(
    yield = data.frame(disease = ids3, value = c(2.5, 5, 7.5), units = "%"),
    fertility = data.frame(disease = ids3, value = c(1, 2, 0), units = "% increase")
  )
}

# Is there a problem in `table` at `row` and `column` (NA = table level)?
has_problem <- function(pr, table, row = NA, column = NA, severity = "error") {
  any(pr$severity == severity & pr$table == table & pr$row %in% row & pr$column %in% column)
}

# A fresh folder under tempdir() with the given tables written as CSV files.
input_dir <- function(name, tables = list()) {
  dir <- file.path(tempdir(), paste0("deconflate-import-", name))
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]], file.path(dir, paste0(nm, ".csv")), row.names = FALSE, na = "")
  }
  dir
}

test_that("the template round-trips through cm_read_inputs()", {
  dir <- input_dir("template")
  expect_message(paths <- cm_template(dir), "Wrote 8 files")
  expect_true(all(file.exists(paths)))
  expect_setequal(basename(paths),
                  c("diseases.csv", "associations.csv", "three_way.csv", "impacts_yield.csv",
                    "impacts_calving_interval.csv", "interactions_yield.csv",
                    "hazard_ratios.csv", "uncertainty.csv"))
  expect_error(cm_template(dir), "already exist")
  expect_message(cm_template(dir, overwrite = TRUE), "Wrote 8 files")

  inp <- cm_read_inputs(dir = dir)
  expect_s3_class(inp, "cm_inputs")
  expect_equal(nrow(inp$problems), 0)
  expect_s3_class(inp$population, "cm_population")
  expect_equal(inp$population$diseases$id, c("LAM", "SCK", "MET"))
  expect_equal(inp$population$diseases$prob[2], 1 - exp(-0.48))
  # The empty three-way and interaction files are ignored.
  expect_null(inp$population$three_way)

  # Two analyses (in file-name order) on one population.
  expect_s3_class(inp$analyses, "cm_analyses")
  expect_equal(names(inp$analyses$models), c("calving_interval", "yield"))
  expect_null(inp$model)
  expect_identical(inp$analyses$population, inp$population)
  for (m in inp$analyses$models) {
    expect_s3_class(m, "cm_model")
    expect_identical(m$diseases, inp$population$diseases)
    expect_identical(m$associations, inp$population$associations)
    expect_null(m$interactions)
  }
  y <- inp$analyses$models$yield$impacts
  expect_equal(y$value, c(4.81, 8.40, 5.61))
  expect_equal(attr(y, "label"), "milk yield loss")
  expect_equal(attr(y, "units"), "% of yield")
  ci <- inp$analyses$models$calving_interval$impacts
  expect_equal(ci$value, c(12, 4, 18))
  expect_equal(attr(ci, "units"), "days")

  # Hazard ratios on the same population.
  expect_s3_class(inp$hr_model, "cm_hr_model")
  expect_equal(inp$hr_model$hazard_ratios$disease, c("LAM", "SCK", "MET"))
  expect_equal(inp$hr_model$hazard_ratios$value, c(1.74, 1.92, 1.50))
  expect_identical(inp$hr_model$population, inp$population)

  # A batch sampler sharing the population draws.
  s <- inp$sampler
  expect_s3_class(s, "cm_batch_sampler")
  expect_equal(names(s$samplers), c("calving_interval", "yield"))
  expect_equal(s$population_keys, c("prob:LAM", "assoc:LAM:SCK"))
  expect_equal(names(attr(s$samplers$yield, "specs")),
               c("prob:LAM", "assoc:LAM:SCK", "impact:LAM", "impact:SCK"))
  expect_equal(names(attr(s$samplers$calving_interval, "specs")),
               c("prob:LAM", "assoc:LAM:SCK", "impact:MET"))
  expect_equal(attr(s$samplers$calving_interval, "specs")[["impact:MET"]]$params$mode, 18)
  expect_equal(attr(s$samplers$yield, "specs")[["assoc:LAM:SCK"]]$params$lower, 0)
  expect_s3_class(s$samplers$yield(1), "cm_model")

  expect_output(print(inp), "Analyses: calving_interval, yield")
  expect_output(print(inp), "Hazard ratios: yes")
  expect_output(print(inp), "batch sampler over 2 analyses")

  res <- deconflate(inp$analyses)
  expect_s3_class(res, "cm_results")
  expect_equal(res$yield$units, "% of yield")
  expect_equal(res$calving_interval$units, "days")
  unlink(dir, recursive = TRUE)
})

test_that("a single impact table typed in R is the analysis 'impacts'", {
  tabs <- supp_tables()
  inp <- cm_read_inputs(
    diseases = tabs$diseases, associations = tabs$associations,
    impacts = data.frame(disease = ids3, value = c(2.5, 5, 7.5), units = "%")
  )
  expect_equal(names(inp$analyses$models), "impacts")
  expect_s3_class(inp$model, "cm_model")
  expect_identical(inp$model, inp$analyses$models$impacts)
  expect_equal(attr(inp$model$impacts, "label"), "impacts")
  expect_equal(attr(inp$model$impacts, "units"), "%")
  expect_equal(names(inp$tables$impacts), "impacts")
  expect_equal(nrow(inp$problems), 0)
  expect_null(inp$sampler)
  expect_null(inp$hr_model)
  expect_equal(deconflate(inp$model)$adjusted$adjusted,
               deconflate(supp_model(c(2.5, 5, 7.5)))$adjusted$adjusted)
  expect_equal(deconflate(inp$model, method = "published")$adjusted$adjusted,
               deconflate(example_supplement(), method = "published")$adjusted$adjusted)
})

test_that("a named list of impact tables gives one analysis each", {
  tabs <- supp_tables()
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = two_impacts())
  expect_equal(names(inp$analyses$models), c("yield", "fertility"))
  expect_null(inp$model)
  expect_equal(attr(inp$analyses$models$fertility$impacts, "units"), "% increase")
  ref <- cm_analyses(supp_population(),
                     yield = cm_impacts(ids3, c(2.5, 5, 7.5)),
                     fertility = cm_impacts(ids3, c(1, 2, 0)))
  res <- deconflate(inp$analyses)
  ref_res <- deconflate(ref)
  for (nm in c("yield", "fertility")) {
    expect_equal(res[[nm]]$adjusted$adjusted, ref_res[[nm]]$adjusted$adjusted)
  }
})

test_that("adjusted_linear impacts are read with their adjustment sets", {
  tabs <- supp_tables()
  imp <- data.frame(disease = ids3, value = c(2.2, 5, 7.5),
                    estimand = c("adjusted_linear", "crude", "crude"),
                    adjusted_for = c("d2", NA, NA))
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imp)
  ref <- supp_model(c(2.2, 5, 7.5), estimand = c("adjusted_linear", "crude", "crude"),
                    adjusted_for = c("d2", NA, NA))
  expect_equal(deconflate(inp$model)$adjusted$adjusted, deconflate(ref)$adjusted$adjusted)

  bad <- imp
  bad$adjusted_for <- c(NA, "d1", NA)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad)
  expect_true(has_problem(pr, "impacts", 1, "adjusted_for"))  # adjusted_linear without a set
  expect_true(has_problem(pr, "impacts", 2, "adjusted_for"))  # adjusted_for with a crude estimate
  bad <- imp
  bad$adjusted_for[1] <- "d9"
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad)
  expect_true(has_problem(pr, "impacts", 1, "adjusted_for"))
  bad <- imp
  bad$estimand[2] <- "adjusted"
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad)
  expect_true(has_problem(pr, "impacts", 2, "estimand"))
})

test_that("separate impact files reuse one population", {
  tabs <- supp_tables()
  dir <- input_dir("separate", list(
    diseases = tabs$diseases, associations = tabs$associations,
    impacts_yield = data.frame(disease = ids3, value = c(2.5, 5, 7.5), label = "milk yield loss",
                               units = "%"),
    impacts_fertility = data.frame(disease = ids3, value = c(1, 2, 0), units = "% increase")
  ))
  inp <- cm_read_inputs(dir = dir)
  m <- inp$analyses$models
  expect_equal(names(m), c("fertility", "yield"))
  expect_identical(m$yield$diseases, m$fertility$diseases)
  expect_identical(m$yield$associations, m$fertility$associations)
  expect_identical(m$yield$diseases, inp$population$diseases)
  expect_identical(m$yield$associations, inp$population$associations)
  expect_equal(attr(m$yield$impacts, "label"), "milk yield loss")
  expect_equal(attr(m$fertility$impacts, "label"), "fertility")
  expect_equal(attr(m$fertility$impacts, "units"), "% increase")

  # The same files given as paths, in a named list.
  inp2 <- cm_read_inputs(
    diseases = file.path(dir, "diseases.csv"), associations = file.path(dir, "associations.csv"),
    impacts = list(yield = file.path(dir, "impacts_yield.csv"),
                   fertility = file.path(dir, "impacts_fertility.csv"))
  )
  expect_equal(names(inp2$analyses$models), c("yield", "fertility"))
  expect_equal(deconflate(inp2$analyses)$yield$adjusted$adjusted,
               deconflate(inp$analyses)$yield$adjusted$adjusted)
  expect_equal(deconflate(inp$analyses)$yield$adjusted$adjusted,
               deconflate(supp_model(c(2.5, 5, 7.5)))$adjusted$adjusted)

  # Arguments given explicitly take precedence over the folder.
  inp3 <- cm_read_inputs(dir = dir, impacts = list(yield = file.path(dir, "impacts_yield.csv")))
  expect_equal(names(inp3$analyses$models), "yield")
  expect_s3_class(inp3$model, "cm_model")
  unlink(dir, recursive = TRUE)
})

test_that("the shipped global dairy files reproduce example_global_dairy()", {
  dir <- system.file("extdata", "global_dairy_2024", package = "deconflate")
  expect_true(nzchar(dir))
  expect_equal(nrow(cm_check_inputs(dir = dir)), 0)
  inp <- cm_read_inputs(dir = dir)
  ref <- example_global_dairy()
  expect_setequal(names(inp$analyses$models), names(ref$models))
  expect_null(inp$model)
  expect_equal(inp$population$diseases$id, ref$population$diseases$id)
  expect_equal(inp$population$diseases$type, ref$population$diseases$type)
  expect_equal(inp$population$diseases$prob, ref$population$diseases$prob, tolerance = 1e-12)

  res <- deconflate(inp$analyses, method = "published", warn = FALSE)
  ref_res <- deconflate(ref, method = "published", warn = FALSE)
  for (nm in names(ref$models)) {
    expect_equal(res[[nm]]$adjusted$disease, ref_res[[nm]]$adjusted$disease)
    expect_equal(res[[nm]]$adjusted$raw, ref_res[[nm]]$adjusted$raw, tolerance = 1e-12)
    expect_equal(res[[nm]]$adjusted$adjusted, ref_res[[nm]]$adjusted$adjusted, tolerance = 1e-10)
    expect_equal(res[[nm]]$label, ref_res[[nm]]$label)
    expect_equal(res[[nm]]$units, ref_res[[nm]]$units)
  }

  # Culling hazard ratios.
  ref_hr <- example_global_dairy_hr()
  expect_s3_class(inp$hr_model, "cm_hr_model")
  expect_equal(inp$hr_model$hazard_ratios$disease, ref_hr$hazard_ratios$disease)
  expect_equal(inp$hr_model$hazard_ratios$value, ref_hr$hazard_ratios$value, tolerance = 1e-12)
  expect_equal(deconflate_hr(inp$hr_model, method = "published", warn = FALSE)$adjusted$adjusted,
               deconflate_hr(ref_hr, method = "published", warn = FALSE)$adjusted$adjusted,
               tolerance = 1e-10)

  # The uncertainty table matches sampler_global_dairy().
  ref_s <- sampler_global_dairy()
  s <- inp$sampler
  expect_s3_class(s, "cm_batch_sampler")
  expect_setequal(s$population_keys, ref_s$population_keys)
  for (nm in names(ref_s$samplers)) {
    sp <- attr(s$samplers[[nm]], "specs")
    rs <- attr(ref_s$samplers[[nm]], "specs")
    expect_setequal(names(sp), names(rs))
    expect_equal(vapply(sp[names(rs)], function(d) d$type, character(1)),
                 vapply(rs, function(d) d$type, character(1)))
    expect_equal(vapply(sp[names(rs)], function(d) d$mean, numeric(1)),
                 vapply(rs, function(d) d$mean, numeric(1)), tolerance = 1e-10)
  }
})

test_that("the shipped example with errors reports each deliberate error", {
  dir <- system.file("extdata", "example_with_errors", package = "deconflate")
  expect_true(nzchar(dir))
  pr <- cm_check_inputs(dir = dir)
  expect_s3_class(pr, "cm_problems")
  expected <- data.frame(
    table = c("diseases", "diseases",
              "associations", "associations", "associations", "associations",
              "impacts_yield", "impacts_yield", "impacts_yield",
              "impacts_fertility", "hazard_ratios", "interactions_culling",
              rep("uncertainty", 9), "uncertainty_fertility", "uncertainty_fertility"),
    row = c(3, 5,
            3, 4, 5, 6,
            2, 3, NA,
            NA, 3, NA,
            2, 3, 4, 5, 6, 7, 9, 10, 11, 2, 3),
    column = c("value", "id",
               "disease2", NA, "value", "adjusted",
               "value", "adjusted_for", NA,
               "units", "value", NA,
               "dist", "key", "key", "key", "key", "dist", "key", "key", "dist", "key", "key"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(expected))) {
    expect_true(has_problem(pr, expected$table[i], expected$row[i], expected$column[i]),
                info = paste(expected$table[i], expected$row[i], expected$column[i]))
  }
  expect_equal(sum(pr$severity == "error"), nrow(expected))
  expect_true(any(grepl("MET", pr$problem[pr$table == "impacts_yield" & is.na(pr$row)])))

  # Every row documented as an error (source or note column) is reported, and
  # no other row is.
  for (f in list.files(dir, pattern = "\\.csv$")) {
    tab <- utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE)
    doc <- if ("note" %in% names(tab)) tab$note else tab$source
    flagged <- which(grepl("^Error", doc))
    lab <- sub("\\.csv$", "", f)
    for (r in flagged) {
      expect_true(any(pr$table == lab & pr$row %in% c(r, NA)), info = paste(f, "row", r))
    }
    reported <- unique(stats::na.omit(pr$row[pr$table == lab]))
    expect_true(all(reported %in% flagged), info = f)
  }

  expect_output(print(pr), "Found 23 problem")
  expect_error(cm_read_inputs(dir = dir), class = "deconflate_input_problems")
})

test_that("all problems in tables typed in R are reported at once", {
  pr <- cm_check_inputs(
    diseases = data.frame(id = c("a", "b", "b", "c"), value = c(0.1, 1.5, 0.2, 0.3)),
    associations = data.frame(disease1 = "a", disease2 = "z", value = -1),
    impacts = data.frame(disease = c("a", "b"), value = c("2", "x"), extra = 1),
    uncertainty = data.frame(key = c("impact:a", "prob:z", "foo"),
                             dist = c("normal", "beta", "pert"), p1 = c(2, 1, 1))
  )
  expect_s3_class(pr, "cm_problems")
  expect_true(has_problem(pr, "diseases", 2, "value"))
  expect_true(has_problem(pr, "diseases", 3, "id"))
  expect_true(has_problem(pr, "associations", 1, "disease2"))
  expect_true(has_problem(pr, "associations", 1, "value"))
  expect_true(has_problem(pr, "impacts", 2, "value"))
  expect_true(has_problem(pr, "impacts", NA, NA))              # no impact for c
  expect_true(has_problem(pr, "impacts", NA, "extra", severity = "note"))
  expect_true(has_problem(pr, "uncertainty", 1, "dist"))       # normal without sd
  expect_false(has_problem(pr, "uncertainty", 1, "key"))       # one analysis: key may omit it
  expect_true(has_problem(pr, "uncertainty", 2, "key"))
  expect_true(has_problem(pr, "uncertainty", 3, "key"))
  expect_true(has_problem(pr, "uncertainty", 3, "dist"))
  expect_error(cm_read_inputs(diseases = data.frame(id = "a", value = 2)), "Found 1 problem",
               class = "deconflate_input_problems")
})

test_that("missing tables, files and folders are reported", {
  pr <- cm_check_inputs(impacts = data.frame(disease = "a", value = 1))
  expect_true(has_problem(pr, "diseases"))
  pr <- cm_check_inputs(diseases = file.path(tempdir(), "no-such-file.csv"))
  expect_true(any(grepl("File not found", pr$problem[pr$table == "diseases"])))
  pr <- cm_check_inputs(dir = file.path(tempdir(), "no-such-folder"))
  expect_true(has_problem(pr, "inputs"))
  tabs <- supp_tables()
  pr <- cm_check_inputs(diseases = tabs$diseases, impacts = unname(two_impacts()))
  expect_true(has_problem(pr, "impacts"))                      # list without names
  pr <- cm_check_inputs(diseases = tabs$diseases,
                        impacts = list("milk yield" = two_impacts()$yield))
  expect_true(has_problem(pr, "impacts"))                      # space in an analysis name
  pr <- cm_check_inputs(diseases = tabs$diseases,
                        hazard_ratios = data.frame(disease = c("d1", "d2"), value = c(1.5, 2),
                                                   adjusted_for = c("d2", NA)))
  expect_true(has_problem(pr, "hazard_ratios", 1, "adjusted_for"))
  expect_true(has_problem(pr, "hazard_ratios", NA, NA))        # no hazard ratio for d3
})

test_that("an impact table with an outcome column gives the migration error", {
  dis <- data.frame(id = c("a", "b"), value = c(0.1, 0.2))
  old <- data.frame(disease = c("a", "b", "a", "b"), outcome = c("yield", "yield", "culling", "culling"),
                    value = c(2, 3, 1.5, 1.8))
  pr <- cm_check_inputs(diseases = dis, impacts = old)
  expect_true(has_problem(pr, "impacts", NA, "outcome"))
  expect_error(cm_read_inputs(diseases = dis, impacts = old), class = "deconflate_input_problems")
  # scale and direction are dropped with a note only
  pr <- cm_check_inputs(diseases = dis,
                        impacts = data.frame(disease = c("a", "b"), value = c(2, 3),
                                             scale = "percent", direction = "decrease"))
  expect_true(has_problem(pr, "impacts", NA, NA, severity = "note"))
  expect_false(any(pr$severity == "error"))
})

test_that("uncertainty keys name the analysis in uncertainty.csv and may omit it in uncertainty_<analysis>.csv", {
  tabs <- supp_tables()
  imps <- two_impacts()
  dir <- input_dir("uncertainty", list(
    diseases = tabs$diseases, associations = tabs$associations,
    impacts_yield = imps$yield, impacts_fertility = imps$fertility,
    uncertainty = data.frame(
      key = c("prob:d1", "assoc:d1:d2", "impact:yield:d1", "impact:fertility:d3"),
      dist = c("beta", "lognormal_ci", "normal", "normal"),
      p1 = c(10, 2, 2.5, 0), p2 = c(90, 1.4, 0.5, 0.2), p3 = c(NA, 2.9, NA, NA)),
    uncertainty_yield = data.frame(key = c("impact:d2", "impact:yield:d3"), dist = "normal",
                                   p1 = c(5, 7.5), p2 = 1)
  ))
  inp <- cm_read_inputs(dir = dir)
  s <- inp$sampler
  expect_s3_class(s, "cm_batch_sampler")
  expect_equal(s$population_keys, c("prob:d1", "assoc:d1:d2"))
  expect_setequal(names(attr(s$samplers$yield, "specs")),
                  c("prob:d1", "assoc:d1:d2", "impact:d1", "impact:d2", "impact:d3"))
  expect_setequal(names(attr(s$samplers$fertility, "specs")),
                  c("prob:d1", "assoc:d1:d2", "impact:d3"))
  expect_equal(attr(s$samplers$yield, "specs")[["impact:d2"]]$mean, 5)
  expect_equal(attr(s$samplers$fertility, "specs")[["impact:d3"]]$mean, 0)
  expect_equal(attr(s$samplers$yield, "specs")[["prob:d1"]]$type, "beta")

  # The same input in uncertainty.csv and uncertainty_yield.csv.
  utils::write.csv(data.frame(key = c("impact:d2", "impact:d1"), dist = "normal", p1 = c(5, 2.5), p2 = 1),
                   file.path(dir, "uncertainty_yield.csv"), row.names = FALSE, na = "")
  pr <- cm_check_inputs(dir = dir)
  expect_true(has_problem(pr, "uncertainty_yield", 2, "key"))
  expect_equal(sum(pr$severity == "error"), 1)
  expect_error(cm_read_inputs(dir = dir), class = "deconflate_input_problems")

  # Population inputs belong in uncertainty.csv, and an analysis's own file
  # cannot name another analysis.
  utils::write.csv(data.frame(key = c("impact:d2", "prob:d2", "assoc:d1:d2", "impact:fertility:d1"),
                              dist = "normal", p1 = c(5, 0.15, 2, 1), p2 = c(1, 0.01, 0.2, 0.1)),
                   file.path(dir, "uncertainty_yield.csv"), row.names = FALSE, na = "")
  pr <- cm_check_inputs(dir = dir)
  expect_false(has_problem(pr, "uncertainty_yield", 1, "key"))
  expect_true(has_problem(pr, "uncertainty_yield", 2, "key"))
  expect_true(has_problem(pr, "uncertainty_yield", 3, "key"))
  expect_true(has_problem(pr, "uncertainty_yield", 4, "key"))
  unlink(dir, recursive = TRUE)
})

test_that("unqualified impact keys need a single analysis", {
  tabs <- supp_tables()
  imps <- two_impacts()
  unc <- data.frame(key = c("impact:yield:d1", "impact:d2"), dist = "normal", p1 = c(2.5, 5), p2 = 1)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = imps, uncertainty = unc)
  expect_false(has_problem(pr, "uncertainty", 1, "key"))
  expect_true(has_problem(pr, "uncertainty", 2, "key"))
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                              impacts = imps, uncertainty = unc),
               class = "deconflate_input_problems")

  # With one analysis, keys with and without the analysis both work.
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = imps["yield"],
                        uncertainty = data.frame(key = c("impact:d1", "impact:yield:d2"),
                                                 dist = "normal", p1 = c(2.5, 5), p2 = 1))
  expect_s3_class(inp$sampler, "cm_sampler")
  expect_equal(names(attr(inp$sampler, "specs")), c("impact:d1", "impact:d2"))
  expect_s3_class(inp$sampler(1), "cm_model")
  expect_output(print(inp), "Uncertain inputs: 2")

  # The shared table can also be given as list element `shared`.
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = imps,
                        uncertainty = list(
                          shared = data.frame(key = c("assoc:d2:d3", "impact:yield:d1"),
                                              dist = "normal", p1 = c(3, 2.5), p2 = c(0.3, 0.5)),
                          fertility = data.frame(key = "impact:d2", dist = "normal", p1 = 2, p2 = 0.5)))
  expect_equal(inp$sampler$population_keys, "assoc:d2:d3")
  expect_equal(names(attr(inp$sampler$samplers$fertility, "specs")), c("assoc:d2:d3", "impact:d2"))

  # An analysis-specific table for an analysis that does not exist.
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imps,
                        uncertainty = list(culling = data.frame(key = "impact:d1", dist = "fixed", p1 = 1)))
  expect_true(has_problem(pr, "uncertainty_culling"))
})

test_that("interaction tables are matched to their analyses", {
  tabs <- supp_tables()
  imps <- two_impacts()
  int <- data.frame(disease1 = "d1", disease2 = "d2", value = 0.5)
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imps,
                        interactions = list(yield = int),
                        uncertainty = data.frame(key = "inter:yield:d2:d1", dist = "normal",
                                                 p1 = 0.5, p2 = 0.1))
  expect_equal(nrow(inp$analyses$models$yield$interactions), 1)
  expect_null(inp$analyses$models$fertility$interactions)
  expect_equal(names(attr(inp$sampler$samplers$yield, "specs")), "inter:d1:d2")
  expect_error(deconflate(inp$analyses$models$yield), class = "deconflate_unsupported")
  expect_output(print(inp$analyses), "1 interactions")

  # interactions_<analysis> for an analysis without an impact table
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imps,
                        interactions = list(culling = int))
  expect_true(has_problem(pr, "interactions_culling"))
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                              impacts = imps, interactions = list(culling = int)),
               class = "deconflate_input_problems")
  # an unnamed interaction table with several analyses
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imps,
                        interactions = int)
  expect_true(has_problem(pr, "interactions"))
})

test_that("three_way.csv is read into the population", {
  tabs <- supp_tables()
  dir <- input_dir("three-way", list(
    diseases = tabs$diseases, associations = tabs$associations,
    three_way = data.frame(disease1 = "d1", disease2 = "d2", disease3 = "d3", ratio = 2,
                           source = "scenario"),
    impacts = data.frame(disease = ids3, value = c(2.5, 5, 7.5)),
    interactions = data.frame(disease1 = "d1", disease2 = "d2", value = 0.5)
  ))
  inp <- cm_read_inputs(dir = dir)
  expect_s3_class(inp$population$three_way, "cm_three_way")
  expect_equal(inp$population$three_way$ratio, 2)
  expect_equal(inp$population$three_way$source, "scenario")
  expect_equal(names(inp$analyses$models), "impacts")
  expect_identical(inp$model$three_way, inp$population$three_way)
  expect_equal(nrow(inp$model$interactions), 1)
  expect_output(print(inp), "Three-way terms: 1")
  ref <- cm_model(supp_population(cm_three_way("d1", "d2", "d3", ratio = 2)),
                  cm_impacts(ids3, c(2.5, 5, 7.5)), cm_interactions("d1", "d2", 0.5))
  expect_equal(deconflate(inp$model, method = "global")$adjusted$adjusted,
               deconflate(ref, method = "global")$adjusted$adjusted, tolerance = 1e-8)

  pr <- cm_check_inputs(diseases = tabs$diseases,
                        three_way = data.frame(disease1 = "d1", disease2 = "d1", disease3 = "d4",
                                               ratio = -1))
  expect_true(has_problem(pr, "three_way", 1, "disease3"))
  expect_true(has_problem(pr, "three_way", 1, NA))
  expect_true(has_problem(pr, "three_way", 1, "ratio"))
  unlink(dir, recursive = TRUE)
})

test_that("adjusted associations are rejected unless used as marginal", {
  tabs <- supp_tables()
  assoc <- tabs$associations
  assoc$adjusted <- c(TRUE, FALSE, FALSE)
  assoc$adjusted_for <- c("parity", NA, NA)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc)
  expect_true(has_problem(pr, "associations", 1, "adjusted"))
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = assoc),
               class = "deconflate_input_problems")

  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc,
                        adjusted_associations = "use_as_marginal")
  expect_false(any(pr$severity == "error"))
  expect_true(has_problem(pr, "associations", 1, "adjusted", severity = "note"))
  expect_message(inp <- cm_read_inputs(diseases = tabs$diseases, associations = assoc,
                                       adjusted_associations = "use_as_marginal"),
                 "as if it were marginal")
  expect_equal(inp$population$adjusted_associations, "use_as_marginal")
  expect_equal(inp$population$associations$adjusted, c(TRUE, FALSE, FALSE))
  expect_equal(inp$population$associations$adjusted_for[1], "parity")
  expect_null(inp$analyses)
  expect_output(print(inp), "used as marginal: 1")
})

test_that("a zero cell in a contingency table gets a note", {
  tabs <- supp_tables()
  assoc <- data.frame(disease1 = "d1", disease2 = "d2", measure = "table",
                      n11 = 0, n10 = 10, n01 = 20, n00 = 70)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc)
  expect_true(has_problem(pr, "associations", 1, NA, severity = "note"))
  expect_false(any(pr$severity == "error"))
  expect_output(print(pr), "Notes on the inputs")
  expect_message(inp <- cm_read_inputs(diseases = tabs$diseases, associations = assoc), "zero cell")
  a <- inp$population$associations
  expect_true(a$corrected)
  expect_equal(a$value, (0.5 * 70.5) / (10.5 * 20.5))
  expect_output(print(inp), "zero-cell correction: 1")

  empty <- assoc
  empty$n10 <- empty$n01 <- empty$n00 <- 0
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = empty)
  expect_true(has_problem(pr, "associations", 1, NA))
})

test_that("hazard ratios can be read without impacts", {
  tabs <- supp_tables()
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        hazard_ratios = data.frame(disease = ids3, value = c(1.5, 2.0, 1.3)))
  expect_null(inp$analyses)
  expect_null(inp$model)
  expect_null(inp$sampler)
  ref <- supp_hr_model()
  expect_equal(inp$hr_model$hazard_ratios$value, ref$hazard_ratios$value)
  expect_equal(deconflate_hr(inp$hr_model, method = "first_order", warn = FALSE)$adjusted$adjusted,
               deconflate_hr(ref, method = "first_order", warn = FALSE)$adjusted$adjusted)
  expect_output(print(inp), "Hazard ratios: yes")
  expect_output(print(cm_check_inputs(diseases = tabs$diseases)), "No problems found")
})

test_that("cm_dist_table builds distributions", {
  tab <- data.frame(key = c("a", "b", "c", "d"), dist = c("pert", "lognormal_ci", "normal", "fixed"),
                    p1 = c(1, 2, 2.6, 3), p2 = c(2, 1.4, 1.4, NA), p3 = c(4, 2.9, 0, NA))
  d <- cm_dist_table(tab)
  expect_equal(names(d), c("a", "b", "c", "d"))
  expect_true(all(vapply(d, inherits, logical(1), "cm_dist")))
  expect_equal(d$a$type, "pert")
  expect_equal(d$a$params$mode, 2)
  expect_equal(d$b$type, "lognormal")
  expect_equal(d$b$params$meanlog, log(2))
  expect_equal(d$c$params$lower, 0)
  expect_equal(d$d$mean, 3)

  f <- file.path(tempdir(), "deconflate-import-dists.csv")
  utils::write.csv(tab, f, row.names = FALSE, na = "")
  d2 <- cm_dist_table(f)
  expect_equal(names(d2), names(d))
  expect_equal(vapply(d2, function(x) x$mean, numeric(1)), vapply(d, function(x) x$mean, numeric(1)))
  unlink(f)

  expect_error(cm_dist_table(data.frame(key = "a", dist = "gamma", p1 = 1)), "Unknown distribution")
  expect_error(cm_dist_table(data.frame(key = c("a", NA), dist = "normal", p1 = 1, p2 = 1)),
               "missing key")
  expect_error(cm_dist_table(data.frame(key = "a", dist = "pert", p1 = 3, p2 = 2, p3 = 4)), "PERT")
  expect_error(cm_dist_table(data.frame(key = "a", p1 = 1)), "needs columns")
  expect_error(cm_dist_table(data.frame(key = character(0), dist = character(0))), "empty")
  expect_error(cm_dist_table(data.frame(key = c("a", "a"), dist = "fixed", p1 = 1:2)),
               "more than once")
})

test_that("CSV columns are read as text, and empty files are empty tables", {
  dir <- input_dir("text-columns")
  writeLines(c("id,value", "T,0.10", "F,0.15", "01,0.20"), file.path(dir, "diseases.csv"))
  file.create(file.path(dir, "three_way.csv"))
  writeLines(c("disease,value", "T,1", "F,2", "01,3"), file.path(dir, "impacts.csv"))
  inp <- cm_read_inputs(dir = dir)
  expect_equal(inp$population$diseases$id, c("T", "F", "01"))
  expect_null(inp$population$three_way)
  expect_equal(inp$model$impacts$value, c(1, 2, 3))
  unlink(dir, recursive = TRUE)
})

test_that("the same input given twice is reported, in one table or in either pair order", {
  tabs <- supp_tables()
  un <- data.frame(key = c("impact:yield:d1", "impact:yield:d1", "assoc:d1:d2", "assoc:d2:d1"),
                   dist = "normal", p1 = c(2.5, 2.5, 2, 2), p2 = c(0.5, 0.5, 0.2, 0.2))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = list(yield = two_impacts()$yield), uncertainty = un)
  expect_true(has_problem(pr, "uncertainty", 2, "key"))
  expect_true(has_problem(pr, "uncertainty", 4, "key"))
  expect_false(has_problem(pr, "uncertainty", 3, "key"))
  # A key in reverse order is matched to the association row.
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = list(yield = two_impacts()$yield),
                        uncertainty = un[c(1, 4), ])
  expect_true("assoc:d1:d2" %in% names(attr(inp$sampler, "specs")))
})

test_that("two interaction tables for one analysis are reported", {
  tabs <- supp_tables()
  it <- data.frame(disease1 = "d1", disease2 = "d2", value = 0.5)
  dir <- input_dir("two-interactions", list(diseases = tabs$diseases,
                                            associations = tabs$associations,
                                            impacts_yield = two_impacts()$yield,
                                            interactions = it, interactions_yield = it))
  pr <- cm_check_inputs(dir = dir)
  expect_true(any(pr$severity == "error" & grepl("^interactions", pr$table)))
})

test_that("an analysis name that partially matches an argument is read correctly", {
  tabs <- supp_tables()
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                        impacts = list(pop = two_impacts()$yield, inter = two_impacts()$fertility))
  expect_equal(names(inp$analyses$models), c("pop", "inter"))
})

test_that("mixed time horizons give a note", {
  d <- data.frame(id = c("a", "b"), value = c(0.1, 0.2), time_horizon = c("year", "lactation"))
  pr <- cm_check_inputs(diseases = d)
  expect_true(has_problem(pr, "diseases", NA, "time_horizon", severity = "note"))
})
