# Reading inputs from CSV files and data frames (cm_read_inputs() and friends).
# Reference values: inst/validation/reference_v040.py (sections 4 and 6a) and
# reference_five_diseases.py.

# The Supplementary File population as tables typed in R.
supp_tables <- function() {
  list(
    diseases = data.frame(id = ids3, value = c(0.10, 0.15, 0.20)),
    associations = data.frame(disease1 = c("d1", "d1", "d2"), disease2 = c("d2", "d3", "d3"),
                              value = c(2, 1, 3)),
    yield = data.frame(disease = ids3, value = c(2.5, 5, 7.5), units = "%")
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

test_that("cm_read_inputs() takes named tables and returns a model", {
  expect_equal(names(formals(cm_read_inputs)),
               c("diseases", "associations", "impacts", "interactions", "three_way",
                 "adjusted_associations"))
  expect_equal(names(formals(cm_check_inputs)), names(formals(cm_read_inputs)))
  tabs <- supp_tables()
  m <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                      impacts = tabs$yield)
  expect_s3_class(m, "cm_model")
  expect_s3_class(attr(m, "problems"), "cm_problems")
  expect_equal(nrow(attr(m, "problems")), 0)
  expect_equal(m$diseases$id, ids3)
  expect_equal(attr(m$impacts, "units"), "%")
  expect_null(attr(m$impacts, "label"))
  expect_equal(attr(m$impacts, "kind"), "additive")
  expect_null(m$distributions)
  expect_null(m$interactions)
  expect_null(m$three_way)
  expect_equal(deconflate(m, n_draws = 0)$adjusted$adjusted,
               deconflate(supp_model(c(2.5, 5, 7.5)), n_draws = 0)$adjusted$adjusted)
  expect_equal(adjust_impacts(m, method = "published")$adjusted$adjusted,
               adjust_impacts(example_supplement(), method = "published")$adjusted$adjusted)

  # diseases and impacts are required; the arguments of 0.3 are gone.
  expect_error(cm_read_inputs(impacts = tabs$yield), "`diseases` is required")
  expect_error(cm_read_inputs(diseases = tabs$diseases), "`impacts` is required")
  expect_error(cm_read_inputs(dir = tempdir()), "unused argument")
  expect_error(cm_check_inputs(dir = tempdir()), "unused argument")
  expect_error(cm_read_inputs(diseases = tabs$diseases, impacts = tabs$yield,
                              hazard_ratios = tabs$yield), "unused argument")
})

test_that("the template writes six files, read back as yield (additive) and culling (event)", {
  expect_false("type" %in% names(formals(cm_template)))
  dir <- input_dir("template")
  expect_message(paths <- cm_template(dir), "Wrote 6 files")
  expect_setequal(basename(paths),
                  c("diseases.csv", "associations.csv", "yield.csv", "yield_interactions.csv",
                    "culling.csv", "three_way.csv"))
  expect_error(cm_template(dir), "already exist")
  expect_message(cm_template(dir, overwrite = TRUE), "Wrote 6 files")
  f <- function(x) file.path(dir, x)
  expect_equal(nrow(cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                                    impacts = f("yield.csv"), interactions = f("yield_interactions.csv"),
                                    three_way = f("three_way.csv"))), 0)

  # Yield: additive impacts, in percent of yield.
  y <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                      impacts = f("yield.csv"))
  expect_equal(nrow(attr(y, "problems")), 0)
  expect_equal(y$diseases$id, c("LAM", "SCK", "MET"))
  expect_equal(y$diseases$prob, c(0.25, 1 - exp(-0.48), 0.10))
  expect_equal(y$impacts$value, c(4.81, 8.40, 5.61))
  expect_equal(attr(y$impacts, "units"), "% of yield")
  expect_equal(attr(y$impacts, "label"), "milk yield loss")
  expect_equal(impact_kind(y$impacts), "additive")
  expect_equal(names(y$distributions), c("prob:LAM", "assoc:LAM:SCK", "impact:LAM", "impact:SCK"))
  expect_equal(y$distributions[["prob:LAM"]]$type, "beta")
  expect_equal(y$distributions[["assoc:LAM:SCK"]]$params$lower, 0)
  expect_output(print(y), "Uncertain inputs \\(with a distribution\\): 4")
  # The empty interaction and three-way files are empty tables.
  y2 <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                       impacts = f("yield.csv"), interactions = f("yield_interactions.csv"),
                       three_way = f("three_way.csv"))
  expect_null(y2$interactions)
  expect_null(y2$three_way)
  ry <- deconflate(y, n_draws = 0)
  expect_equal(ry$method, "simultaneous")
  expect_length(ry$notes, 0)
  expect_null(ry$draws)
  # reference_v040.py, 6a
  expect_equal(ry$adjusted$adjusted, c(2.8705956676, 7.8191139087, 3.1698811769), tolerance = 1e-9)
  expect_equal(ry$totals$raw_sum, 4.965719508828417, tolerance = 1e-12)
  expect_equal(ry$totals$adjusted_total, 4.015413117941525, tolerance = 1e-9)
  expect_equal(deconflate(y2, n_draws = 0)$adjusted$adjusted, ry$adjusted$adjusted)

  # Culling: event impacts (hazard ratios and a risk ratio).
  cu <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                       impacts = f("culling.csv"))
  expect_equal(impact_kind(cu$impacts), "event")
  expect_equal(cu$impacts$measure, c("HR", "HR", "RR"))
  expect_equal(cu$impacts$estimand, rep("snapshot_crude", 3))
  expect_equal(cu$impacts$value, c(1.74, 1.92, 1.45))
  expect_equal(attr(cu$impacts, "label"), "culling")
  expect_null(attr(cu$impacts, "units"))
  expect_equal(names(cu$distributions), c("prob:LAM", "assoc:LAM:SCK", "impact:LAM"))
  expect_equal(cu$distributions[["impact:LAM"]]$type, "lognormal")
  expect_equal(cu$distributions[["impact:LAM"]]$params$meanlog, log(1.74))
  expect_output(print(cu), "event impacts: use event_model = TRUE")
  expect_error(deconflate(cu, n_draws = 0), "event_model = TRUE", class = "deconflate_unsupported")
  expect_warning(rc <- deconflate(cu, event_model = TRUE, overall_risk = 0.25, n_draws = 0), NA)
  expect_s3_class(rc, "cm_event_result")
  expect_equal(rc$method, "snapshot")
  expect_equal(rc$adjusted$adjusted, c(1.5143397997, 1.7819332524, 1.2012868074), tolerance = 1e-7)
  expect_equal(rc$attributable$summary$attributable, 0.0732820406, tolerance = 1e-8)
  unlink(dir, recursive = TRUE)
})

test_that("tables without associations are read; deconflate() then needs one", {
  tabs <- supp_tables()
  m <- cm_read_inputs(diseases = tabs$diseases, impacts = tabs$yield)
  expect_s3_class(m, "cm_model")
  expect_null(m$associations)
  expect_equal(pair_tables(m)$status, rep("unknown", 3))
  expect_true(all(is.na(pair_tables(m)$measure)))
  expect_output(print(m), "3 unknown")
  expect_error(deconflate(m), "No association estimates", class = "deconflate_unsupported")
  expect_error(deconflate(m), "screen_associations")
  # An associations table without rows is the same as none.
  m0 <- cm_read_inputs(diseases = tabs$diseases, impacts = tabs$yield,
                       associations = tabs$associations[0, ])
  expect_null(m0$associations)
  expect_equal(nrow(cm_check_inputs(diseases = tabs$diseases, impacts = tabs$yield)), 0)
  # The sensitivity tools run without associations.
  sc <- screen_associations(m)
  expect_s3_class(sc, "cm_screen")
  expect_true(all(sc$status == "unknown"))
})

test_that("adjusted_linear impacts are read with their adjustment sets", {
  tabs <- supp_tables()
  imp <- data.frame(disease = ids3, value = c(2.2, 5, 7.5),
                    estimand = c("adjusted_linear", "crude", "crude"),
                    adjusted_for = c("d2", NA, NA))
  m <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imp)
  ref <- supp_model(c(2.2, 5, 7.5), estimand = c("adjusted_linear", "crude", "crude"),
                    adjusted_for = c("d2", NA, NA))
  expect_equal(deconflate(m, n_draws = 0)$adjusted$adjusted,
               deconflate(ref, n_draws = 0)$adjusted$adjusted)

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

test_that("the shipped global dairy files match example_global_dairy()", {
  f <- function(x) extdata("global_dairy_2024", x)
  expect_true(file.exists(f("culling.csv")))
  expect_false(file.exists(f("hazard_ratios.csv")))
  files <- c(yield = "impacts_yield.csv", fertility = "impacts_fertility.csv", culling = "culling.csv")
  for (outcome in names(files)) {
    expect_equal(nrow(cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                                      impacts = f(files[[outcome]]))), 0, info = outcome)
    m <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                        impacts = f(files[[outcome]]))
    ref <- example_global_dairy(outcome)
    expect_equal(m$diseases$id, ref$diseases$id)
    expect_equal(m$diseases$type, ref$diseases$type)
    expect_equal(m$diseases$prob, ref$diseases$prob, tolerance = 1e-12)
    km <- pair_key(m$associations$disease1, m$associations$disease2)
    kr <- pair_key(ref$associations$disease1, ref$associations$disease2)
    expect_setequal(km, kr)
    expect_equal(m$associations$value[match(kr, km)], ref$associations$value, tolerance = 1e-12)
    expect_equal(m$impacts$disease, ref$impacts$disease)
    expect_equal(m$impacts$value, ref$impacts$value, tolerance = 1e-12)
    expect_equal(m$impacts$measure, ref$impacts$measure)
    expect_equal(m$impacts$estimand, ref$impacts$estimand)
    expect_equal(attr(m$impacts, "label"), attr(ref$impacts, "label"))
    expect_equal(attr(m$impacts, "units"), attr(ref$impacts, "units"))
    expect_equal(impact_kind(m$impacts), if (outcome == "culling") "event" else "additive")
    # The distributions in the tables are those of the 2024 analysis.
    dm <- m$distributions
    dr <- ref$distributions
    expect_setequal(names(dm), names(dr))
    expect_equal(vapply(dm[names(dr)], function(d) d$type, character(1)),
                 vapply(dr, function(d) d$type, character(1)))
    expect_equal(vapply(dm[names(dr)], function(d) d$mean, numeric(1)),
                 vapply(dr, function(d) d$mean, numeric(1)), tolerance = 1e-10)
  }
  # With the unlisted pairs independent (as in the 2024 analysis), the
  # published approximation gives the analysis values.
  m <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                      impacts = f("impacts_yield.csv"))
  expect_equal(adjust_impacts(with_independent_pairs(m), method = "published")$adjusted$adjusted,
               adjust_impacts(global_dairy_analyses("analysis")$yield, method = "published")$adjusted$adjusted,
               tolerance = 1e-10)
})

test_that("the shipped example with errors reports each deliberate error", {
  f <- function(x) extdata("example_with_errors", x)
  expect_setequal(list.files(extdata("example_with_errors"), pattern = "\\.csv$"),
                  c("diseases.csv", "associations.csv", "impacts_yield.csv", "culling.csv"))
  pr <- cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                        impacts = f("impacts_yield.csv"))
  expect_s3_class(pr, "cm_problems")
  expect_named(pr, c("table", "row", "column", "severity", "problem"))
  expected <- data.frame(
    table = c(rep("diseases", 4), rep("associations", 10), rep("impacts", 4)),
    row = c(2, 3, 4, 5,
            NA, 2, 3, 4, 5, 6, 7, 8, 9, 9,
            2, 3, 3, NA),
    column = c("dist", "value", "dist", "id",
               NA, "dist", "disease2", NA, "value", "adjusted", "measure", "measure", "measure",
               "disease2",
               "value", "adjusted_for", "dist", "units"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(expected))) {
    expect_true(has_problem(pr, expected$table[i], expected$row[i], expected$column[i]),
                info = paste(expected$table[i], expected$row[i], expected$column[i]))
  }
  # Each error is reported once (checked row by row against R/import.R).
  expect_equal(sum(pr$severity == "error"), nrow(expected))
  # The note: a point value outside its distribution (MET, uniform from 6 to 8).
  expect_equal(sum(pr$severity == "note"), 1)
  expect_true(has_problem(pr, "impacts", 4, "dist", severity = "note"))
  # The decimal comma is a value problem, not a missing value as well.
  expect_equal(sum(pr$table == "impacts" & pr$row %in% 2), 1)
  # The count columns are reported once, for the table; the retired
  # measures once per row, with what to do instead.
  counts <- pr$problem[pr$table == "associations" & is.na(pr$row)]
  expect_length(counts, 1)
  expect_match(counts, "n11 \\* n00")
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 7 & pr$column %in% "measure"],
               "measure OR")
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 8 & pr$column %in% "measure"],
               "odds ratio of 1")
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 9 & pr$column %in% "measure"],
               "leave the pair out")
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 5], "Leave the pair out")
  expect_match(pr$problem[pr$table == "diseases" & pr$row %in% 2], "p2 \\(sd\\)")
  # Every row documented as an error in the note column is reported (in its
  # row, or at table level), and no undocumented row is.
  files <- c(diseases = "diseases.csv", associations = "associations.csv", impacts = "impacts_yield.csv")
  for (lab in names(files)) {
    tab <- utils::read.csv(f(files[[lab]]), stringsAsFactors = FALSE)
    flagged <- which(grepl("^Error", tab$note))
    expect_gt(length(flagged), 0)
    for (r in flagged) {
      expect_true(any(pr$table == lab & pr$row %in% c(r, NA) & pr$severity == "error"),
                  info = paste(files[[lab]], "row", r))
    }
    reported <- unique(stats::na.omit(pr$row[pr$table == lab & pr$severity == "error"]))
    expect_true(all(reported %in% flagged), info = files[[lab]])
  }
  expect_output(print(pr), "Found 18 problem")
  expect_error(cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                              impacts = f("impacts_yield.csv")),
               "diseases, row 3, column 'value'", class = "deconflate_input_problems")

  # culling.csv: event impacts with three mistakes.
  pc <- cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                        impacts = f("culling.csv"))
  expect_true(has_problem(pc, "impacts", 2, "estimand"))       # no estimand
  expect_true(has_problem(pc, "impacts", 3, "value"))          # HR of 0
  expect_true(has_problem(pc, "impacts", 4, "adjusted_for"))   # stratified without a set
  expect_equal(sum(pc$severity == "error" & pc$table == "impacts"), 3)
  expect_match(pc$problem[pc$table == "impacts" & pc$row %in% 2], "snapshot_crude or snapshot_stratified")
  # The other tables report the same problems as before.
  expect_equal(sum(pc$severity == "error"), 17)
  expect_false(any(pc$severity == "note"))
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
  expect_error(cm_read_inputs(diseases = data.frame(id = "a", value = 2),
                              impacts = data.frame(disease = "a", value = 1)),
               "Found 1 problem", class = "deconflate_input_problems")
})

test_that("missing tables and files are reported", {
  pr <- cm_check_inputs(impacts = data.frame(disease = "a", value = 1))
  expect_true(has_problem(pr, "diseases"))
  pr <- cm_check_inputs(diseases = file.path(tempdir(), "no-such-file.csv"))
  expect_true(any(grepl("File not found", pr$problem[pr$table == "diseases"])))
  expect_true(has_problem(pr, "impacts"))                      # no impact table
  tabs <- supp_tables()
  pr <- cm_check_inputs(diseases = tabs$diseases)
  expect_true(any(grepl("impact table is required", pr$problem[pr$table == "impacts"])))
  pr <- cm_check_inputs(diseases = list(1, 2), impacts = tabs$yield)
  expect_true(any(grepl("data frame or the path", pr$problem[pr$table == "diseases"])))
  expect_error(cm_read_inputs(diseases = file.path(tempdir(), "no-such-file.csv"), impacts = tabs$yield),
               "File not found", class = "deconflate_input_problems")
  # Event impacts: an adjustment set for a crude estimate, and a missing row.
  pr <- cm_check_inputs(diseases = tabs$diseases,
                        impacts = data.frame(disease = c("d1", "d2"), value = c(1.5, 2), measure = "HR",
                                             estimand = "snapshot_crude", adjusted_for = c("d2", NA)))
  expect_true(has_problem(pr, "impacts", 1, "adjusted_for"))
  expect_true(has_problem(pr, "impacts", NA, NA))              # no impact for d3
  expect_match(pr$problem[pr$table == "impacts" & is.na(pr$row)], "1 for ratios")
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

test_that("retired association measures and count columns are errors", {
  tabs <- supp_tables()
  assoc <- data.frame(disease1 = c("d1", "d1", "d2"), disease2 = c("d2", "d3", "d3"),
                      value = c(NA, 1, NA), measure = c("table", "independent", "unknown"),
                      n11 = c(10, NA, NA), n10 = c(20, NA, NA), n01 = c(30, NA, NA), n00 = c(40, NA, NA))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield)
  expect_true(has_problem(pr, "associations", NA, NA))         # the count columns, once
  expect_equal(sum(pr$table == "associations" & is.na(pr$row)), 1)
  expect_true(has_problem(pr, "associations", 1, "measure"))
  expect_true(has_problem(pr, "associations", 2, "measure"))
  expect_true(has_problem(pr, "associations", 3, "measure"))
  # A retired row is not also reported as a missing value.
  expect_false(has_problem(pr, "associations", 1, "value"))
  expect_false(has_problem(pr, "associations", 3, "value"))
  expect_equal(sum(pr$severity == "error"), 4)
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield),
               "no longer", class = "deconflate_input_problems")
  # Count columns alone (with an odds ratio in `value`) are still an error.
  assoc2 <- data.frame(disease1 = "d1", disease2 = "d2", value = 2, n11 = 10, n10 = 20, n01 = 30, n00 = 40)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc2, impacts = tabs$yield)
  expect_equal(sum(pr$severity == "error"), 1)
  expect_true(has_problem(pr, "associations", NA, NA))

  # The same in R.
  expect_error(cm_associations("d1", "d2", 1, measure = "independent"), "odds ratio of 1",
               class = "deconflate_unsupported")
  expect_error(cm_associations("d1", "d2", 1, measure = "unknown"), "leave the pair out",
               class = "deconflate_unsupported")
  expect_error(cm_associations("d1", "d2", 1, measure = "table"), "n11 \\* n00",
               class = "deconflate_unsupported")
  expect_error(cm_associations("d1", "d2", NA), "missing for some rows")
  expect_error(cm_associations("d1", "d2"), "Give each association a `value`")
  expect_false(any(c("n11", "n10", "n01", "n00", "zero_cell") %in% names(formals(cm_associations))))
  expect_false("missing_associations" %in% c(names(formals(cm_population)), names(formals(cm_model))))
  # An odds ratio of 1 states independence.
  pt <- pair_tables(supp_population())
  expect_equal(pt$p11[pt$disease1 == "d1" & pt$disease2 == "d3"], 0.10 * 0.20)
  expect_equal(pt$status, rep("specified", 3))
})

test_that("event impact tables have a measure column and are checked", {
  tabs <- supp_tables()
  ev <- data.frame(disease = ids3, value = c(1.5, 1.3, 0.04), measure = c("HR", "RR", "RD"),
                   estimand = "snapshot_crude", label = "culling")
  m <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = ev)
  expect_equal(attr(m$impacts, "kind"), "event")
  expect_equal(m$impacts$measure, c("HR", "RR", "RD"))
  expect_output(print(m), "Measures: HR: 1; RD: 1; RR: 1")
  expect_s3_class(deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0), "cm_event_result")
  # Additive impact tables have no measures.
  y <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = tabs$yield)
  expect_true(all(is.na(y$impacts$measure)))

  # Mistakes in event tables.
  bad <- data.frame(disease = ids3, value = c(1.5, -1, 1.2), measure = c("HR", NA, "RD"),
                    estimand = c(NA, "snapshot_crude", "snapshot_stratified"),
                    adjusted_for = c(NA, NA, "d9"))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad)
  expect_true(has_problem(pr, "impacts", 1, "estimand"))       # estimand required
  expect_true(has_problem(pr, "impacts", 2, "measure"))        # a measure in every row
  expect_true(has_problem(pr, "impacts", 3, "value"))          # RD outside (-1, 1)
  expect_true(has_problem(pr, "impacts", 3, "adjusted_for"))   # unknown disease
  bad2 <- data.frame(disease = ids3, value = c(0, 1.3, 1.2), measure = c("HR", "HZ", "OR"),
                     estimand = c("snapshot_crude", "snapshot_crude", "crude"),
                     adjusted_for = c("d2", NA, NA))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad2)
  expect_true(has_problem(pr, "impacts", 1, "value"))          # a ratio must be positive
  expect_true(has_problem(pr, "impacts", 1, "adjusted_for"))   # set given for snapshot_crude
  expect_true(has_problem(pr, "impacts", 2, "measure"))        # unknown measure
  expect_true(has_problem(pr, "impacts", 3, "estimand"))       # an additive estimand
  expect_match(pr$problem[pr$table == "impacts" & pr$row %in% 3 & pr$column %in% "estimand"],
               "for event impacts")

  # Event estimands in an additive table, and interactions with event impacts.
  mix <- data.frame(disease = ids3, value = c(1.5, 2, 1.3), estimand = "snapshot_crude")
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = mix)
  expect_true(has_problem(pr, "impacts", 1, "estimand"))
  expect_match(pr$problem[pr$table == "impacts" & pr$row %in% 1], "add a measure column")
  int <- data.frame(disease1 = "d1", disease2 = "d2", value = 0.5)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = ev,
                        interactions = int)
  expect_true(has_problem(pr, "interactions", NA, NA))
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = ev,
                              interactions = int),
               "additive impacts only", class = "deconflate_input_problems")
  # The same mistakes made in R.
  expect_error(cm_impacts(ids3, c(1.5, 2, 1.3), estimand = "snapshot_crude"), "measure",
               class = "deconflate_unsupported")
  expect_error(cm_impacts(ids3, c(1.5, 2, 1.3), measure = "HR"), "State the estimand")
  expect_error(cm_impacts(ids3, c(1.5, 2, 1.3), measure = "HR", estimand = "crude"), "Invalid `estimand`")
  expect_error(cm_impacts(ids3, c(1.5, 0, 1.3), measure = "HR", estimand = "snapshot_crude"),
               "must be positive")
  expect_error(cm_impacts(ids3, c(1.5, 2, 1.3), measure = "HZ", estimand = "snapshot_crude"),
               "Invalid `measure`")
  expect_error(cm_model(supp_population(),
                        cm_impacts(ids3, c(1.5, 2, 1.3), measure = "HR", estimand = "snapshot_crude"),
                        cm_interactions("d1", "d2", 0.5)),
               "additive impacts only", class = "deconflate_unsupported")
})

test_that("interactions are read for additive impacts", {
  tabs <- supp_tables()
  int <- data.frame(disease1 = "d1", disease2 = "d2", value = 0.5, dist = "normal", p1 = 0.5, p2 = 0.1)
  m <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = tabs$yield,
                      interactions = int)
  expect_equal(nrow(m$interactions), 1)
  expect_equal(names(m$distributions), "inter:d1:d2")
  expect_output(print(m), "Interactions: 1")
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "global")
  expect_match(res$notes, "because of interactions")
  ref <- supp_model(c(2.5, 5, 7.5), interactions = cm_interactions("d1", "d2", 0.5))
  expect_equal(res$adjusted$adjusted, deconflate(ref, n_draws = 0)$adjusted$adjusted, tolerance = 1e-10)
  res_s <- deconflate(m, method = "simultaneous", n_draws = 0)
  expect_equal(res_s$method, "global")
  expect_match(res_s$notes, "instead of the simultaneous method")
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = tabs$yield,
                        interactions = data.frame(disease1 = c("d1", "d2"), disease2 = c("d1", "d1"),
                                                  value = c(1, 2)))
  expect_true(has_problem(pr, "interactions", 1, "disease2"))   # the same disease twice
})

test_that("three_way is read into the model", {
  tabs <- supp_tables()
  dir <- input_dir("three-way", list(
    diseases = tabs$diseases, associations = tabs$associations,
    three_way = data.frame(disease1 = "d1", disease2 = "d2", disease3 = "d3", ratio = 2,
                           source = "scenario"),
    yield = data.frame(disease = ids3, value = c(2.5, 5, 7.5)),
    yield_interactions = data.frame(disease1 = "d1", disease2 = "d2", value = 0.5)
  ))
  f <- function(x) file.path(dir, x)
  m <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                      impacts = f("yield.csv"), interactions = f("yield_interactions.csv"),
                      three_way = f("three_way.csv"))
  expect_s3_class(m$three_way, "cm_three_way")
  expect_equal(m$three_way$ratio, 2)
  expect_equal(m$three_way$source, "scenario")
  expect_equal(nrow(m$interactions), 1)
  expect_output(print(m), "Three-way terms: 1")
  ref <- cm_model(supp_population(cm_three_way("d1", "d2", "d3", ratio = 2)),
                  cm_impacts(ids3, c(2.5, 5, 7.5)), cm_interactions("d1", "d2", 0.5))
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "global")
  expect_match(res$notes, "interactions and three-way terms")
  expect_equal(res$adjusted$adjusted, deconflate(ref, n_draws = 0)$adjusted$adjusted, tolerance = 1e-8)

  pr <- cm_check_inputs(diseases = tabs$diseases, impacts = tabs$yield,
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
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield)
  expect_true(has_problem(pr, "associations", 1, "adjusted"))
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield),
               class = "deconflate_input_problems")

  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield,
                        adjusted_associations = "use_as_marginal")
  expect_false(any(pr$severity == "error"))
  expect_true(has_problem(pr, "associations", 1, "adjusted", severity = "note"))
  expect_message(m <- cm_read_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield,
                                     adjusted_associations = "use_as_marginal"),
                 "as if it were marginal")
  expect_equal(nrow(attr(m, "problems")), 1)
  expect_equal(m$adjusted_associations, "use_as_marginal")
  expect_equal(m$associations$adjusted, c(TRUE, FALSE, FALSE))
  expect_equal(m$associations$adjusted_for[1], "parity")
  expect_output(print(m), "used as marginal: 1")
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

  # The keys can be given to cm_model() as distributions.
  dd <- cm_dist_table(data.frame(key = c("impact:d1", "assoc:d2:d1"), dist = c("normal", "lognormal_ci"),
                                 p1 = c(2.5, 2), p2 = c(0.5, 1.4), p3 = c(NA, 2.9)))
  m <- cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5)), distributions = dd)
  expect_equal(names(m$distributions), c("impact:d1", "assoc:d1:d2"))

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
  writeLines(c("disease,value", "T,1", "F,2", "01,3"), file.path(dir, "losses.csv"))
  m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"), impacts = file.path(dir, "losses.csv"),
                      three_way = file.path(dir, "three_way.csv"))
  expect_equal(m$diseases$id, c("T", "F", "01"))
  expect_null(m$three_way)
  expect_equal(m$impacts$value, c(1, 2, 3))
  unlink(dir, recursive = TRUE)
})

test_that("mixed time horizons give a note", {
  d <- data.frame(id = c("a", "b"), value = c(0.1, 0.2), time_horizon = c("year", "lactation"))
  pr <- cm_check_inputs(diseases = d, impacts = data.frame(disease = c("a", "b"), value = c(1, 2)))
  expect_true(has_problem(pr, "diseases", NA, "time_horizon", severity = "note"))
  expect_false(any(pr$severity == "error"))
})

# ---- Uncertainty in the tables ---------------------------------------------------

test_that("distributions are given in the rows of the values they describe", {
  tabs <- supp_tables()
  dis <- cbind(tabs$diseases, dist = c("beta", NA, NA), p1 = c(10, NA, NA), p2 = c(90, NA, NA))
  assoc <- cbind(tabs$associations, dist = c("lognormal_ci", NA, "pert"),
                 p1 = c(2, NA, 2), p2 = c(1.4, NA, 3), p3 = c(2.9, NA, 4.5))
  tw <- data.frame(disease1 = "d3", disease2 = "d1", disease3 = "d2", ratio = 1.5,
                   dist = "lognormal_ci", p1 = 1.5, p2 = 0.8, p3 = 2.8)
  imp <- cbind(tabs$yield, dist = c("normal", "pert", "fixed"),
               p1 = c(2.5, 3, 7.5), p2 = c(0.5, 5, NA), p3 = c(NA, 8, NA))
  int <- data.frame(disease1 = "d2", disease2 = "d1", value = 0.5, dist = "normal", p1 = 0.5, p2 = 0.1)
  pr <- cm_check_inputs(diseases = dis, associations = assoc, impacts = imp, interactions = int,
                        three_way = tw)
  expect_equal(nrow(pr), 0)
  m <- cm_read_inputs(diseases = dis, associations = assoc, impacts = imp, interactions = int,
                      three_way = tw)
  # Point values are the central inputs.
  expect_equal(m$diseases$value, c(0.10, 0.15, 0.20))
  expect_equal(m$three_way$ratio, 1.5)
  expect_equal(m$impacts$value, c(2.5, 5, 7.5))
  d <- m$distributions
  expect_equal(names(d), c("prob:d1", "assoc:d1:d2", "assoc:d2:d3", "three:d3:d1:d2",
                           "impact:d1", "impact:d2", "impact:d3", "inter:d2:d1"))
  expect_equal(d[["prob:d1"]]$type, "beta")
  expect_equal(d[["assoc:d1:d2"]]$type, "lognormal")
  expect_equal(d[["assoc:d1:d2"]]$params$meanlog, log(2))
  expect_equal(d[["three:d3:d1:d2"]]$params$meanlog, log(1.5))
  expect_equal(d[["impact:d2"]]$type, "pert")
  expect_equal(d[["impact:d3"]]$type, "fixed")
  # The draws use them (the global method, because of the interaction).
  res <- deconflate(m, n_draws = 20, seed = 1)
  expect_equal(res$method, "global")
  expect_equal(colnames(res$draws$params), names(d))
  expect_gt(length(unique(res$draws$params[, "three:d3:d1:d2"])), 1)
  expect_true(all(res$draws$params[, "impact:d3"] == 7.5))
  expect_equal(nrow(res$draws$params), 20 - res$draws$n_rejected)

  # Event impacts can have distributions too.
  ev <- data.frame(disease = ids3, value = c(1.5, 2, 1.3), measure = "HR", estimand = "snapshot_crude",
                   dist = c("lognormal_ci", NA, NA), p1 = c(1.5, NA, NA), p2 = c(1.2, NA, NA),
                   p3 = c(1.9, NA, NA))
  me <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = ev)
  expect_equal(names(me$distributions), "impact:d1")
  expect_equal(nrow(attr(me, "problems")), 0)
})

test_that("problems with distributions are reported per row", {
  dis <- data.frame(id = c("a", "b", "c"), value = c(0.1, 0.2, 0.3),
                    dist = c(NA, "uniform", NA), p1 = c(0.05, 0.25, NA), p2 = c(NA, 0.4, NA))
  assoc <- data.frame(disease1 = c("a", "b"), disease2 = c("b", "c"), value = c(2, 2),
                      dist = c("normal", "lognormal"), p1 = c(2, 0.7), p2 = c(0.1, 0.2))
  imp <- data.frame(disease = c("a", "b", "c"), value = c(1, 2, 3),
                    dist = c("gamma", "normal", "pert"), p1 = c(1, 2, 1), p2 = c(1, "0,5", 5),
                    p3 = c(NA, NA, 2))
  pr <- cm_check_inputs(diseases = dis, associations = assoc, impacts = imp)
  expect_true(has_problem(pr, "diseases", 1, "dist"))                     # p1 without dist
  expect_true(has_problem(pr, "diseases", 2, "dist", severity = "note"))  # 0.2 outside (0.25, 0.4)
  expect_false(has_problem(pr, "associations", 1, "dist"))
  expect_false(has_problem(pr, "associations", 2, "dist"))
  expect_true(has_problem(pr, "impacts", 1, "dist"))                      # unknown distribution
  expect_true(has_problem(pr, "impacts", 2, "p2"))                        # decimal comma
  expect_false(has_problem(pr, "impacts", 2, "dist"))                     # reported once
  expect_true(has_problem(pr, "impacts", 3, "dist"))                      # pert with max < mode
  expect_equal(sum(pr$severity == "error"), 4)
  expect_equal(sum(pr$severity == "note"), 1)
  expect_output(print(pr), "Found 4 problem")
  # Notes alone do not stop reading; they are printed as a message.
  ok <- dis[2:3, ]
  expect_message(m <- cm_read_inputs(diseases = ok, impacts = data.frame(disease = c("b", "c"), value = 1:2)),
                 "Notes on the inputs")
  expect_equal(nrow(attr(m, "problems")), 1)
  expect_equal(names(m$distributions), "prob:b")
  expect_output(print(cm_check_inputs(diseases = ok[2, ], impacts = data.frame(disease = "c", value = 1))),
                "No problems found")
})

test_that("the five-disease example reads without problems and matches the reference", {
  # Reference values: inst/validation/reference_five_diseases.py and
  # reference_v040.py (section 4).
  f <- function(x) extdata("five_diseases", x)
  expect_false(file.exists(f("hazard_ratios.csv")))
  for (imp in c("impacts_yield.csv", "impacts_calving_interval.csv", "impacts_welfare.csv", "culling.csv")) {
    pr <- cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                          impacts = f(imp), three_way = f("three_way.csv"),
                          interactions = if (imp == "impacts_welfare.csv") f("interactions_welfare.csv"))
    expect_equal(nrow(pr), 0, info = imp)
  }
  ids <- c("LAM", "MAS", "MET", "SCK", "RP")
  y <- read_five("impacts_yield.csv")
  expect_equal(y$diseases$id, ids)
  expect_equal(y$diseases$prob, c(0.25, 1 - exp(-0.3), 0.10, 0.35, 0.06), tolerance = 1e-12)
  pt <- pair_tables(y)
  expect_equal(nrow(pt), 10)
  expect_true(all(pt$status == "specified"))
  expect_setequal(unique(pt$measure), c("OR", "RR", "RD", "cond_prob", "phi"))
  expect_equal(attr(y$impacts, "units"), "% of yield")
  expect_null(y$three_way)

  # Distributions from every table, mixed with point values.
  y3 <- read_five("impacts_yield.csv", three_way = TRUE)
  expect_equal(nrow(y3$three_way), 1)
  expect_setequal(names(y3$distributions),
                  c("prob:LAM", "prob:MAS", "prob:MET", "prob:SCK",
                    "assoc:LAM:MAS", "assoc:LAM:SCK", "assoc:MET:RP", "assoc:MET:SCK",
                    "assoc:MAS:SCK", "assoc:MAS:MET", "assoc:LAM:MET", "three:LAM:MAS:SCK",
                    "impact:LAM", "impact:MAS", "impact:MET", "impact:SCK"))
  types <- vapply(y3$distributions, function(d) d$type, character(1))
  expect_true(all(c("beta", "pert", "uniform", "lognormal", "normal", "fixed") %in% types))
  ci <- read_five("impacts_calving_interval.csv")
  expect_setequal(setdiff(names(ci$distributions), names(y$distributions)[!startsWith(names(y$distributions), "impact:")]),
                  c("impact:LAM", "impact:MAS", "impact:MET", "impact:RP"))
  w <- read_five("impacts_welfare.csv", three_way = TRUE, interactions = "interactions_welfare.csv")
  expect_true(all(c("impact:LAM", "impact:MAS", "impact:SCK", "impact:RP", "inter:LAM:MAS",
                    "inter:MET:RP") %in% names(w$distributions)))

  # Deterministic results.
  r <- deconflate(y, n_draws = 0)
  expect_equal(r$method, "simultaneous")
  expect_equal(r$adjusted$adjusted,
               c(4.13927424, 2.35580211, 4.5140575, 1.36629496, 2.66542027), tolerance = 1e-7)
  expect_equal(r$totals$adjusted_total, 2.7349337441102657, tolerance = 1e-9)
  rp <- adjust_impacts(y, method = "published")
  expect_equal(rp$adjusted$adjusted,
               c(3.97507962, 2.37379544, 4.36087399, 1.58917725, 2.84352728), tolerance = 1e-7)
  # With the three-way term the global method is used; the pairwise tables,
  # and so the additive results, are the same.
  r3 <- deconflate(y3, n_draws = 0)
  expect_equal(r3$method, "global")
  expect_match(r3$notes, "three-way terms")
  expect_equal(r3$adjusted$adjusted, r$adjusted$adjusted, tolerance = 1e-7)
  rc <- deconflate(ci, n_draws = 0)
  expect_equal(rc$method, "simultaneous")
  expect_equal(rc$adjusted$adjusted, c(10.27962506, 4.40963526, 15.61622693, 4, 5.06713035),
               tolerance = 1e-7)
  expect_equal(rc$totals$adjusted_total, 6.978453892271114, tolerance = 1e-9)
  expect_error(adjust_impacts(ci, method = "published"), class = "deconflate_unsupported")
  rw <- deconflate(w, n_draws = 0)
  expect_equal(rw$method, "global")
  expect_equal(rw$adjusted$adjusted, c(8.77541245, 4.00171284, 3.18388605, 1.02183469, 2.14413033),
               tolerance = 1e-7)
  expect_equal(rw$totals$adjusted_total, 4.182930749703578, tolerance = 1e-8)
  expect_error(adjust_impacts(w, method = "simultaneous"), class = "deconflate_unsupported")

  # Culling: event impacts with mixed measures and estimands (and the
  # three-way term in the joint distribution).
  cu <- read_five("culling.csv", three_way = TRUE)
  expect_equal(cu$impacts$measure, c("HR", "HR", "RR", "OR", "RD"))
  expect_equal(cu$impacts$estimand, c("snapshot_crude", "snapshot_stratified", "snapshot_crude",
                                      "snapshot_crude", "snapshot_stratified"))
  expect_equal(cu$impacts$adjusted_for[c(2, 5)], c("LAM;SCK", "all"))
  expect_true(all(c("impact:LAM", "impact:MET", "impact:RP") %in% names(cu$distributions)))
  expect_error(deconflate(cu, n_draws = 0), class = "deconflate_unsupported")
  ev <- deconflate(cu, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_s3_class(ev, "cm_event_result")
  expect_equal(ev$adjusted$adjusted, c(1.58344491, 1.57873827, 1.08072624, 1.19114596, 1.50176712),
               tolerance = 1e-6)
  expect_equal(ev$adjusted$change[c(1, 2)], ev$adjusted$adjusted[1:2] / c(1.74, 1.6) - 1)
  expect_true(all(is.na(ev$adjusted$change[3:5])))
  expect_equal(ev$attributable$summary$attributable, 0.0697502441, tolerance = 1e-8)
  expect_equal(ev$attributable$by_disease$attributable,
               c(0.0248718025, 0.0256191332, 0.0015698203, 0.0123310667, 0.0053584214), tolerance = 1e-8)
  expect_error(adjust_event(cu, method = "first_order", overall_risk = 0.25), class = "deconflate_unsupported")

  # Thresholds used in run_all_features.R.
  th <- cm_threshold(y, "impact:SCK", c(0, 2.5), conclusion = "sign", method = "simultaneous")
  t1 <- th$thresholds[th$thresholds$status == "threshold", ]
  expect_equal(nrow(t1), 1)
  expect_equal(t1$threshold, 1.1935832731450715, tolerance = 1e-6)
  th <- cm_threshold(y, "impact:MET", c(0, 15), conclusion = "total", target = 3,
                     method = "simultaneous")
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 9.730627246165994,
               tolerance = 1e-6)
  th <- cm_threshold(y, "prob:SCK", c(0.1, 0.6), conclusion = "rank",
                     diseases = c("MET", "SCK"), method = "simultaneous")
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 0.33002218486278845,
               tolerance = 1e-6)
  # The published approximation is not a method of cm_threshold(); its pole
  # for SCK (reference_five_diseases.py) is checked on the internal
  # adjust_impacts() instead.
  expect_error(cm_threshold(y, "impact:SCK", c(-4, 2.5), conclusion = "sign", diseases = "SCK",
                            method = "published"),
               class = "deconflate_unsupported")
  pole <- -1.4328526785573266
  A <- rp$conflation$A
  expect_equal(-as.vector((A - diag(nrow(A))) %*% rp$adjusted$raw)[rp$adjusted$disease == "SCK"], pole,
               tolerance = 1e-8)
  sck_published <- function(v) {
    my <- y
    my$impacts$value[my$impacts$disease == "SCK"] <- v
    r <- adjust_impacts(my, method = "published", warn = FALSE)
    r$adjusted$adjusted[r$adjusted$disease == "SCK"]
  }
  expect_lt(sck_published(pole - 1e-3), -1000)
  expect_gt(sck_published(pole + 1e-3), 1000)
  skip_on_cran()
  th <- cm_threshold(w, "inter:LAM:MAS", c(-5, 20), conclusion = "change", target = -0.1)
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 7.607122195588001,
               tolerance = 1e-6)
  # Draws of every input: almost every draw is usable.
  dr <- deconflate(w, n_draws = 30, seed = 1)$draws
  expect_lte(dr$n_rejected, 3)
  expect_equal(colnames(dr$params), names(w$distributions))
})
