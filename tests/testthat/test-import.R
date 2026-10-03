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

test_that("the default template is one analysis, read into a model and a sampler", {
  dir <- input_dir("template")
  expect_message(paths <- cm_template(dir), "Wrote 5 files")
  expect_setequal(basename(paths),
                  c("diseases.csv", "associations.csv", "three_way.csv", "impacts.csv",
                    "interactions.csv"))
  expect_error(cm_template(dir), "already exist")
  expect_message(cm_template(dir, overwrite = TRUE), "Wrote 5 files")

  inp <- cm_read_inputs(dir = dir)
  expect_equal(nrow(inp$problems), 0)
  expect_equal(inp$population$diseases$id, c("LAM", "SCK", "MET"))
  expect_equal(inp$population$diseases$prob[2], 1 - exp(-0.48))
  expect_null(inp$population$three_way)
  expect_null(inp$hr_model)
  expect_s3_class(inp$model, "cm_model")
  expect_equal(names(inp$analyses$models), "impacts")
  expect_null(inp$model$interactions)
  expect_equal(inp$model$impacts$value, c(4.81, 8.40, 5.61))
  expect_equal(attr(inp$model$impacts, "units"), "% of yield")

  # A single-model sampler: Latin hypercube and importance sampling work.
  s <- inp$sampler
  expect_s3_class(s, "cm_sampler")
  expect_equal(names(attr(s, "specs")), c("prob:LAM", "assoc:LAM:SCK", "impact:LAM", "impact:SCK"))
  specs <- attr(s, "specs")
  prop <- list("impact:LAM" = dist_mixture(specs[["impact:LAM"]], dist_normal(6, 1.5),
                                           weights = c(0.5, 0.5)))
  mc <- cm_monte_carlo(s, 40, proposal = prop, seed = 1)
  expect_s3_class(mc, "cm_mc")
  expect_lt(mc$ess, 40)
  mc_lhs <- cm_monte_carlo(s, 40, sampling = "lhs", lhs_replicates = 4, seed = 1)
  expect_equal(mc_lhs$n_blocks, 4)
  expect_output(print(inp), "Analyses: impacts")
  unlink(dir, recursive = TRUE)
})

test_that("the analyses template gives two analyses and a batch sampler", {
  dir <- input_dir("template-analyses")
  expect_message(paths <- cm_template(dir, type = "analyses"), "Wrote 6 files")
  expect_setequal(basename(paths),
                  c("diseases.csv", "associations.csv", "three_way.csv", "impacts_yield.csv",
                    "impacts_calving_interval.csv", "interactions_yield.csv"))
  inp <- cm_read_inputs(dir = dir)
  expect_equal(nrow(inp$problems), 0)
  expect_s3_class(inp$analyses, "cm_analyses")
  expect_equal(names(inp$analyses$models), c("calving_interval", "yield"))
  expect_null(inp$model)
  expect_null(inp$hr_model)
  expect_identical(inp$analyses$population, inp$population)
  for (m in inp$analyses$models) {
    expect_identical(m$diseases, inp$population$diseases)
    expect_identical(m$associations, inp$population$associations)
  }
  expect_equal(inp$analyses$models$calving_interval$impacts$value, c(12, 4, 18))
  s <- inp$sampler
  expect_s3_class(s, "cm_batch_sampler")
  expect_equal(s$population_keys, c("prob:LAM", "assoc:LAM:SCK"))
  expect_equal(names(attr(s$samplers$calving_interval, "specs")),
               c("prob:LAM", "assoc:LAM:SCK", "impact:MET"))
  expect_output(print(inp), "batch sampler over 2 analyses")
  res <- deconflate(inp$analyses)
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

  # The distributions in the tables match sampler_global_dairy().
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
    table = c(rep("diseases", 4), rep("associations", 6), rep("impacts_yield", 4),
              "impacts_fertility", "impacts_fertility", "hazard_ratios", "interactions_culling",
              "uncertainty"),
    row = c(2, 3, 4, 5,
            2, 3, 4, 5, 6, 7,
            2, 3, 3, NA,
            NA, 3, 3, NA, NA),
    column = c("dist", "value", "dist", "id",
               "dist", "disease2", NA, "value", "adjusted", "dist",
               "value", "adjusted_for", "dist", NA,
               "units", "p2", "value", NA, NA),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(expected))) {
    expect_true(has_problem(pr, expected$table[i], expected$row[i], expected$column[i]),
                info = paste(expected$table[i], expected$row[i], expected$column[i]))
  }
  expect_equal(sum(pr$severity == "error"), nrow(expected))
  expect_true(any(grepl("MET", pr$problem[pr$table == "impacts_yield" & is.na(pr$row)])))
  expect_true(any(grepl("no longer read", pr$problem[pr$table == "uncertainty"])))
  # A decimal comma in a parameter is reported once, not again as a
  # distribution with a missing parameter.
  expect_false(has_problem(pr, "impacts_fertility", 3, "dist"))
  # Notes: a point value outside its distribution, and distributions given
  # for hazard ratios (ignored).
  expect_true(has_problem(pr, "impacts_fertility", 2, "dist", severity = "note"))
  expect_true(has_problem(pr, "hazard_ratios", NA, "dist", severity = "note"))

  # Every row documented as an error (or note) in the note column is
  # reported, and no other row is.
  for (f in list.files(dir, pattern = "\\.csv$")) {
    tab <- utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE)
    lab <- sub("\\.csv$", "", f)
    for (sev in c("error", "note")) {
      flagged <- which(grepl(if (sev == "error") "^Error" else "^Note", tab$note))
      for (r in flagged) {
        expect_true(any(pr$table == lab & pr$row %in% c(r, NA) & pr$severity == sev),
                    info = paste(f, "row", r))
      }
      reported <- unique(stats::na.omit(pr$row[pr$table == lab & pr$severity == sev]))
      expect_true(all(reported %in% flagged), info = paste(f, sev))
    }
  }

  expect_output(print(pr), "Found 19 problem")
  expect_error(cm_read_inputs(dir = dir), class = "deconflate_input_problems")
})

test_that("all problems in tables typed in R are reported at once", {
  pr <- cm_check_inputs(
    diseases = data.frame(id = c("a", "b", "b", "c"), value = c(0.1, 1.5, 0.2, 0.3)),
    associations = data.frame(disease1 = "a", disease2 = "z", value = -1),
    impacts = data.frame(disease = c("a", "b"), value = c("2", "x"), extra = 1,
                         dist = c("normal", NA), p1 = c(2, 1))
  )
  expect_s3_class(pr, "cm_problems")
  expect_true(has_problem(pr, "diseases", 2, "value"))
  expect_true(has_problem(pr, "diseases", 3, "id"))
  expect_true(has_problem(pr, "associations", 1, "disease2"))
  expect_true(has_problem(pr, "associations", 1, "value"))
  expect_true(has_problem(pr, "impacts", 2, "value"))
  expect_true(has_problem(pr, "impacts", NA, NA))              # no impact for c
  expect_true(has_problem(pr, "impacts", NA, "extra", severity = "note"))
  expect_true(has_problem(pr, "impacts", 1, "dist"))           # normal without sd
  expect_true(has_problem(pr, "impacts", 2, "dist"))           # p1 without a distribution
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
                                                   estimand = "snapshot_crude",
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

test_that("interaction tables are matched to their analyses", {
  tabs <- supp_tables()
  imps <- two_impacts()
  int <- data.frame(disease1 = "d1", disease2 = "d2", value = 0.5)
  inp <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imps,
                        interactions = list(yield = cbind(int, dist = "normal", p1 = 0.5, p2 = 0.1)))
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
                        hazard_ratios = data.frame(disease = ids3, value = c(1.5, 2.0, 1.3),
                                                   estimand = "snapshot_crude"))
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

# ---- Uncertainty in the tables (deconflate 0.3.0) -----------------------------

test_that("distributions are given in the rows of the values they describe", {
  tabs <- supp_tables()
  dir <- input_dir("dists", list(
    diseases = cbind(tabs$diseases, dist = c("beta", NA, NA), p1 = c(10, NA, NA),
                     p2 = c(90, NA, NA)),
    associations = cbind(tabs$associations, dist = c("lognormal_ci", NA, "pert"),
                         p1 = c(2, NA, 2), p2 = c(1.4, NA, 3), p3 = c(2.9, NA, 4.5)),
    three_way = data.frame(disease1 = "d3", disease2 = "d1", disease3 = "d2", ratio = 1.5,
                           dist = "lognormal_ci", p1 = 1.5, p2 = 0.8, p3 = 2.8),
    impacts_yield = cbind(two_impacts()$yield, dist = c("normal", "pert", NA),
                          p1 = c(2.5, 3, NA), p2 = c(0.5, 5, NA), p3 = c(NA, 8, NA)),
    impacts_fertility = cbind(two_impacts()$fertility, dist = c(NA, NA, "fixed"),
                              p1 = c(NA, NA, 0)),
    interactions_yield = data.frame(disease1 = "d2", disease2 = "d1", value = 0.5,
                                    dist = "normal", p1 = 0.5, p2 = 0.1)
  ))
  pr <- cm_check_inputs(dir = dir)
  expect_equal(nrow(pr), 0)
  inp <- cm_read_inputs(dir = dir)
  # Point values are used by the deterministic methods.
  expect_equal(inp$population$diseases$value, c(0.10, 0.15, 0.20))
  expect_equal(inp$population$three_way$ratio, 1.5)
  s <- inp$sampler
  expect_s3_class(s, "cm_batch_sampler")
  expect_setequal(s$population_keys, c("prob:d1", "assoc:d1:d2", "assoc:d2:d3", "three:d3:d1:d2"))
  expect_setequal(names(attr(s$samplers$yield, "specs")),
                  c(s$population_keys, "impact:d1", "impact:d2", "inter:d2:d1"))
  expect_setequal(names(attr(s$samplers$fertility, "specs")), c(s$population_keys, "impact:d3"))
  sp <- attr(s$samplers$yield, "specs")
  expect_equal(sp[["prob:d1"]]$type, "beta")
  expect_equal(sp[["assoc:d1:d2"]]$type, "lognormal")
  expect_equal(sp[["assoc:d1:d2"]]$params$meanlog, log(2))
  expect_equal(sp[["three:d3:d1:d2"]]$params$meanlog, log(1.5))
  expect_equal(sp[["impact:d2"]]$type, "pert")
  expect_equal(attr(s$samplers$fertility, "specs")[["impact:d3"]]$type, "fixed")
  # The three-way ratio is drawn with the population inputs.
  m <- s$samplers$yield(1, values = c("three:d3:d1:d2" = 2.5))
  expect_equal(m$three_way$ratio, 2.5)
  mc <- cm_monte_carlo(s, 10, method = "global", seed = 1)
  p <- mc$analyses$yield$params
  expect_true("three:d3:d1:d2" %in% names(p))
  expect_equal(p[["three:d3:d1:d2"]], mc$analyses$fertility$params[["three:d3:d1:d2"]][
    match(p$draw, mc$analyses$fertility$params$draw)])
  expect_true(all(mc$analyses$fertility$params[["impact:d3"]] == 0))

  # One analysis typed in R: a single sampler.
  inp1 <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                         impacts = data.frame(disease = ids3, value = c(2.5, 5, 7.5),
                                              dist = c("normal", NA, NA), p1 = c(2.5, NA, NA),
                                              p2 = c(0.5, NA, NA)))
  expect_s3_class(inp1$sampler, "cm_sampler")
  expect_equal(names(attr(inp1$sampler, "specs")), "impact:d1")
  expect_output(print(inp1), "Uncertain inputs: 1")
  unlink(dir, recursive = TRUE)
})

test_that("uncertainty files from deconflate 0.2 give a migration error", {
  tabs <- supp_tables()
  dir <- input_dir("old-uncertainty", list(
    diseases = tabs$diseases,
    impacts_yield = two_impacts()$yield,
    uncertainty = data.frame(key = "impact:yield:d1", dist = "normal", p1 = 2.5, p2 = 0.5),
    uncertainty_yield = data.frame(key = "impact:d2", dist = "normal", p1 = 5, p2 = 1)
  ))
  pr <- cm_check_inputs(dir = dir)
  expect_true(has_problem(pr, "uncertainty"))
  expect_true(has_problem(pr, "uncertainty_yield"))
  expect_match(pr$problem[pr$table == "uncertainty"], "no longer read")
  expect_error(cm_read_inputs(dir = dir), "no longer read", class = "deconflate_input_problems")
  expect_false("uncertainty" %in% names(formals(cm_read_inputs)))
  unlink(dir, recursive = TRUE)
})

test_that("problems with distributions are reported per row", {
  dis <- data.frame(id = c("a", "b", "c"), value = c(0.1, 0.2, 0.3),
                    dist = c(NA, "uniform", NA), p1 = c(0.05, 0.25, NA), p2 = c(NA, 0.4, NA))
  assoc <- data.frame(disease1 = c("a", "a", "b"), disease2 = c("b", "c", "c"),
                      value = c(NA, NA, 2), measure = c("table", "independent", "OR"),
                      n11 = c(10, NA, NA), n10 = c(20, NA, NA), n01 = c(30, NA, NA),
                      n00 = c(40, NA, NA), dist = c("normal", "normal", "lognormal"),
                      p1 = c(1, 1, 0.7), p2 = c(0.1, 0.1, 0.2))
  imp <- data.frame(disease = c("a", "b", "c"), value = c(1, 2, 3),
                    dist = c("gamma", "normal", "pert"), p1 = c(1, 2, 1), p2 = c(1, "0,5", 5),
                    p3 = c(NA, NA, 2))
  hr <- data.frame(disease = c("a", "b", "c"), value = c(1.5, 1.2, 1), estimand = "snapshot_crude",
                   dist = c("lognormal_ci", NA, NA), p1 = c(1.5, NA, NA), p2 = c(1.2, NA, NA),
                   p3 = c(1.9, NA, NA))
  pr <- cm_check_inputs(diseases = dis, associations = assoc, impacts = imp, hazard_ratios = hr)
  expect_true(has_problem(pr, "diseases", 1, "dist"))                     # p1 without dist
  expect_true(has_problem(pr, "diseases", 2, "dist", severity = "note"))  # 0.2 outside (0.25, 0.4)
  expect_true(has_problem(pr, "associations", 1, "dist"))                 # a table of counts
  expect_true(has_problem(pr, "associations", 2, "dist"))                 # independent
  expect_false(has_problem(pr, "associations", 3, "dist"))
  expect_true(has_problem(pr, "impacts", 1, "dist"))                      # unknown distribution
  expect_true(has_problem(pr, "impacts", 2, "p2"))                        # decimal comma
  expect_false(has_problem(pr, "impacts", 2, "dist"))                     # reported once
  expect_true(has_problem(pr, "impacts", 3, "dist"))                      # pert with max < mode
  expect_true(has_problem(pr, "hazard_ratios", NA, "dist", severity = "note"))
  expect_equal(sum(pr$severity == "error"), 6)

  # Without impact tables, distributions are noted (Monte Carlo adjusts impacts).
  pr <- cm_check_inputs(diseases = dis[2:3, ])
  expect_true(has_problem(pr, "inputs", NA, "dist", severity = "note"))
})

test_that("the five-disease example reads without problems and matches the reference", {
  # Reference values: inst/validation/reference_five_diseases.py
  dir <- system.file("extdata", "five_diseases", package = "deconflate")
  expect_true(nzchar(dir))
  expect_equal(nrow(cm_check_inputs(dir = dir)), 0)
  inp <- cm_read_inputs(dir = dir)
  ids <- c("LAM", "MAS", "MET", "SCK", "RP")
  expect_equal(inp$population$diseases$id, ids)
  expect_equal(inp$population$diseases$prob, c(0.25, 1 - exp(-0.3), 0.10, 0.35, 0.06),
               tolerance = 1e-12)
  expect_equal(nrow(pair_tables(inp$population)), 10)
  expect_equal(nrow(inp$population$three_way), 1)
  expect_equal(names(inp$analyses$models), c("calving_interval", "welfare", "yield"))
  expect_equal(attr(inp$analyses$models$yield$impacts, "units"), "% of yield")
  expect_s3_class(inp$hr_model, "cm_hr_model")
  expect_equal(inp$hr_model$hazard_ratios$estimand,
               c("snapshot_crude", "snapshot_stratified", "snapshot_crude", "snapshot_crude",
                 "snapshot_stratified"))

  # Distributions from every table, mixed with point values.
  s <- inp$sampler
  expect_s3_class(s, "cm_batch_sampler")
  expect_setequal(s$population_keys,
                  c("prob:LAM", "prob:MAS", "prob:MET", "prob:SCK",
                    "assoc:LAM:MAS", "assoc:LAM:SCK", "assoc:MET:RP", "assoc:MET:SCK",
                    "assoc:MAS:SCK", "assoc:MAS:MET", "assoc:LAM:MET", "three:LAM:MAS:SCK"))
  own <- function(nm) setdiff(names(attr(s$samplers[[nm]], "specs")), s$population_keys)
  expect_setequal(own("yield"), c("impact:LAM", "impact:MAS", "impact:MET", "impact:SCK"))
  expect_setequal(own("calving_interval"), c("impact:LAM", "impact:MAS", "impact:MET", "impact:RP"))
  expect_setequal(own("welfare"), c("impact:LAM", "impact:MAS", "impact:SCK", "impact:RP",
                                    "inter:LAM:MAS", "inter:MET:RP"))
  types <- vapply(attr(s$samplers$yield, "specs"), function(d) d$type, character(1))
  expect_true(all(c("beta", "pert", "uniform", "lognormal", "normal", "fixed") %in% types))

  # Deterministic results.
  m <- inp$analyses$models
  r <- deconflate(m$yield)
  expect_equal(r$adjusted$adjusted,
               c(4.13927424, 2.35580211, 4.5140575, 1.36629496, 2.66542027), tolerance = 1e-7)
  expect_equal(r$totals$adjusted_total, 2.7349337441102657, tolerance = 1e-9)
  rp <- deconflate(m$yield, method = "published")
  expect_equal(rp$adjusted$adjusted,
               c(3.97507962, 2.37379544, 4.36087399, 1.58917725, 2.84352728), tolerance = 1e-7)
  expect_equal(deconflate(m$yield, method = "global")$adjusted$adjusted, r$adjusted$adjusted,
               tolerance = 1e-7)
  rc <- deconflate(m$calving_interval)
  expect_equal(rc$adjusted$adjusted, c(10.27962506, 4.40963526, 15.61622693, 4, 5.06713035),
               tolerance = 1e-7)
  expect_error(deconflate(m$calving_interval, method = "published"), class = "deconflate_unsupported")
  rw <- deconflate(m$welfare, method = "global")
  expect_equal(rw$adjusted$adjusted, c(8.77541245, 4.00171284, 3.18388605, 1.02183469, 2.14413033),
               tolerance = 1e-7)
  expect_equal(rw$totals$adjusted_total, 4.182930749703578, tolerance = 1e-8)
  expect_error(deconflate(m$welfare), class = "deconflate_unsupported")
  expect_error(deconflate_hr(inp$hr_model, method = "published"), class = "deconflate_unsupported")
  expect_s3_class(deconflate_hr(inp$hr_model), "cm_hr_result")

  # Thresholds used in run_all_features.R.
  th <- cm_threshold(m$yield, "impact:SCK", c(0, 2.5), conclusion = "sign", method = "simultaneous")
  t1 <- th$thresholds[th$thresholds$status == "threshold", ]
  expect_equal(nrow(t1), 1)
  expect_equal(t1$threshold, 1.1935832731450715, tolerance = 1e-6)
  th <- cm_threshold(m$yield, "impact:MET", c(0, 15), conclusion = "total", target = 3,
                     method = "simultaneous")
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 9.730627246165994,
               tolerance = 1e-6)
  th <- cm_threshold(m$yield, "prob:SCK", c(0.1, 0.6), conclusion = "rank",
                     diseases = c("MET", "SCK"), method = "simultaneous")
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 0.33002218486278845,
               tolerance = 1e-6)
  th <- cm_threshold(m$yield, "impact:SCK", c(-4, 2.5), conclusion = "sign", diseases = "SCK",
                     method = "published")
  expect_equal(th$thresholds$status, "discontinuity")
  expect_lt(th$thresholds$lower, -1.4328526785573266)
  expect_gt(th$thresholds$upper, -1.4328526785573266)
  skip_on_cran()
  th <- cm_threshold(m$welfare, "inter:LAM:MAS", c(-5, 20), conclusion = "change", target = -0.1)
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 7.607122195588001,
               tolerance = 1e-6)
  # Monte Carlo over all analyses: almost every draw is usable.
  mc <- cm_monte_carlo(s, 30, method = "global", seed = 1)
  expect_s3_class(mc, "cm_mc_batch")
  expect_lte(max(vapply(mc$analyses, function(a) a$n_rejected, numeric(1))), 3)
})
