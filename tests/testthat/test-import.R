# Reading inputs from CSV files and data frames (cm_read_inputs() and friends).
# Reference values: inst/validation/reference_v040.py (sections 4 and 6a),
# reference_five_diseases.py (the ten pair rows of the five-disease example)
# and reference_v041.py (the five-disease example of 0.4).

# The Supplementary File population as tables typed in R.
supp_tables <- function() {
  list(
    diseases = data.frame(disease = ids3, value = c(0.10, 0.15, 0.20)),
    associations = data.frame(disease = c("d1", "d1", "d2"), associated_with = c("d2", "d3", "d3"),
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

test_that("cm_read_inputs() takes three named tables and returns a model", {
  expect_equal(names(formals(cm_read_inputs)), c("diseases", "associations", "impacts", "..."))
  expect_equal(names(formals(cm_check_inputs)), names(formals(cm_read_inputs)))
  tabs <- supp_tables()
  m <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations,
                      impacts = tabs$yield)
  expect_s3_class(m, "cm_model")
  expect_s3_class(attr(m, "problems"), "cm_problems")
  expect_equal(nrow(attr(m, "problems")), 0)
  expect_equal(m$diseases$disease, ids3)
  expect_equal(m$associations$row, 1:3)
  expect_equal(m$impacts$row, 1:3)
  expect_equal(attr(m$impacts, "units"), "%")
  expect_null(attr(m$impacts, "label"))
  expect_equal(attr(m$impacts, "kind"), "additive")
  expect_null(m$distributions)
  expect_equal(deconflate(m, n_draws = 0)$adjusted$adjusted,
               deconflate(supp_model(c(2.5, 5, 7.5)), n_draws = 0)$adjusted$adjusted)

  # diseases and impacts are required; the arguments of earlier versions
  # give messages.
  expect_error(cm_read_inputs(impacts = tabs$yield), "`diseases` is required")
  expect_error(cm_read_inputs(diseases = tabs$diseases), "`impacts` is required")
  expect_error(cm_read_inputs(dir = tempdir()), "Unknown argument")
  expect_error(cm_check_inputs(dir = tempdir()), "Unknown argument")
  expect_error(cm_read_inputs(diseases = tabs$diseases, impacts = tabs$yield, three_way = "x"),
               "OR_ratio", class = "deconflate_unsupported")
  expect_error(cm_read_inputs(diseases = tabs$diseases, impacts = tabs$yield, interactions = "x"),
               "rows of the impact table", class = "deconflate_unsupported")
  expect_error(cm_check_inputs(diseases = tabs$diseases, impacts = tabs$yield,
                               adjusted_associations = "use_as_marginal"), "no longer used")
})

test_that("the template writes four files, read back as yield (additive) and culling (event)", {
  dir <- input_dir("template")
  expect_message(paths <- cm_template(dir), "Wrote 4 files")
  expect_setequal(basename(paths), c("diseases.csv", "associations.csv", "yield.csv", "culling.csv"))
  expect_error(cm_template(dir), "already exist")
  expect_message(cm_template(dir, overwrite = TRUE), "Wrote 4 files")
  f <- function(x) file.path(dir, x)
  # The columns of each file, in the documented order.
  expect_equal(names(utils::read.csv(f("diseases.csv")))[1:4], c("disease", "measure", "value", "label"))
  expect_equal(names(utils::read.csv(f("associations.csv")))[1:6],
               c("disease", "associated_with", "among", "adjusted_for", "measure", "value"))
  expect_equal(names(utils::read.csv(f("yield.csv")))[1:4], c("disease", "among", "adjusted_for", "value"))
  expect_equal(nrow(cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                                    impacts = f("yield.csv"))), 0)

  # Yield: additive impacts, in percent of yield.
  y <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                      impacts = f("yield.csv"))
  expect_equal(nrow(attr(y, "problems")), 0)
  expect_equal(y$diseases$disease, c("LAM", "SCK", "MET"))
  expect_equal(y$diseases$label, c("lameness", "subclinical ketosis", "metritis"))
  expect_equal(y$diseases$prob, c(0.25, 1 - exp(-0.48), 0.10))
  expect_equal(y$impacts$value, c(4.81, 8.40, 5.61))
  expect_equal(attr(y$impacts, "units"), "% of yield")
  expect_equal(attr(y$impacts, "label"), "milk yield loss")
  expect_equal(impact_kind(y$impacts), "additive")
  expect_equal(names(y$distributions), c("prob:LAM", "associations:1", "impacts:1", "impacts:2"))
  expect_equal(y$distributions[["prob:LAM"]]$type, "beta")
  expect_equal(y$distributions[["associations:1"]]$params$lower, 0)
  expect_output(print(y), "Uncertain inputs \\(with a distribution\\): 4")
  ry <- deconflate(y, n_draws = 0)
  expect_equal(ry$method, "pairwise shortcut")
  expect_length(ry$notes, 0)
  expect_null(ry$draws)
  # reference_v040.py, 6a
  expect_equal(ry$adjusted$adjusted, c(2.8705956676, 7.8191139087, 3.1698811769), tolerance = 1e-9)
  expect_equal(ry$totals$raw_total, 4.965719508828417, tolerance = 1e-12)
  expect_equal(ry$totals$adjusted_total, 4.015413117941525, tolerance = 1e-9)

  # Culling: event impacts (hazard ratios and a risk ratio).
  cu <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                       impacts = f("culling.csv"))
  expect_equal(impact_kind(cu$impacts), "event")
  expect_equal(cu$impacts$measure, c("HR", "HR", "RR"))
  expect_equal(cu$impacts$value, c(1.74, 1.92, 1.45))
  expect_equal(attr(cu$impacts, "label"), "culling")
  expect_null(attr(cu$impacts, "units"))
  expect_equal(names(cu$distributions), c("prob:LAM", "associations:1", "impacts:1"))
  expect_equal(cu$distributions[["impacts:1"]]$type, "lognormal")
  expect_equal(cu$distributions[["impacts:1"]]$params$meanlog, log(1.74))
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
  expect_error(deconflate(m), "cm_sensitivity")
  # An associations table without rows is the same as none.
  m0 <- cm_read_inputs(diseases = tabs$diseases, impacts = tabs$yield,
                       associations = tabs$associations[0, ])
  expect_null(m0$associations)
  expect_equal(nrow(cm_check_inputs(diseases = tabs$diseases, impacts = tabs$yield)), 0)
  # The sensitivity tools run without associations.
  sc <- cm_sensitivity(m, inputs = "unknown")
  expect_s3_class(sc, "cm_sensitivity")
  expect_equal(sum(startsWith(sc$inputs$input, "pair:")), 3)
})

test_that("adjusted impacts are read with their adjustment sets", {
  tabs <- supp_tables()
  imp <- data.frame(disease = ids3, value = c(2.2, 5, 7.5), adjusted_for = c("d2", NA, NA))
  m <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imp)
  ref <- supp_model(c(2.2, 5, 7.5), adjusted_for = c("d2", NA, NA))
  expect_equal(deconflate(m, n_draws = 0)$adjusted$adjusted,
               deconflate(ref, n_draws = 0)$adjusted$adjusted)
  # A name that is not a disease of the model is a factor outside the model:
  # a note, and the row is used as crude.
  bad <- imp
  bad$adjusted_for[1] <- "d9"
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad)
  expect_true(has_problem(pr, "impacts", 1, "adjusted_for", severity = "note"))
  expect_false(any(pr$severity == "error"))
  expect_match(pr$problem[1], "the row is used as crude")
  # The estimand column of earlier versions is ignored, with a note.
  old <- cbind(imp, estimand = c("adjusted_linear", "crude", "crude"))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = old)
  expect_true(has_problem(pr, "impacts", NA, "estimand", severity = "note"))
  expect_false(any(pr$severity == "error"))
})

test_that("rows about combinations and subgroups are read", {
  tabs <- supp_tables()
  # (The d1-d3 pair is described by the prevalence of the combination.)
  assoc <- data.frame(disease = c("d1", "d2", "d1", "d1;d3"),
                      associated_with = c("d2", "d3", "d2", NA),
                      among = c(NA, NA, "d3", NA),
                      measure = c("OR", "OR", "OR_ratio", "prevalence"),
                      value = c(2, 3, 1.5, 0.015))
  imp <- data.frame(disease = c(ids3, "d1;d2", "d1"), among = c(NA, NA, NA, NA, "no others"),
                    adjusted_for = c(NA, NA, NA, "all", NA), value = c(2.5, 5, 7.5, 1, 2.4),
                    dist = c(NA, NA, NA, NA, "normal"), p1 = c(NA, NA, NA, NA, 2.4), p2 = c(NA, NA, NA, NA, 0.5))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc, impacts = imp)
  expect_equal(nrow(pr), 0)
  m <- cm_read_inputs(diseases = tabs$diseases, associations = assoc, impacts = imp)
  expect_equal(m$associations$among, c(NA, NA, "d3", NA))
  expect_equal(pair_tables(m)$status, c("specified", "other rows", "specified"))
  expect_equal(m$impacts$disease, c(ids3, "d1;d2", "d1"))
  expect_output(print(m), "interactions: d1;d2")
  expect_output(print(m), "Association rows: 4")
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "full model")
  expect_equal(res$adjusted$disease, c(ids3, "d1;d2"))
  expect_equal(res$adjusted$rows, c(2, 1, 1, 1))
  # The d1 row among animals with no other disease is fitted by weighted
  # least squares beside the exact rows.
  expect_equal(res$estimates$fitted[1:4], c(2.5, 5, 7.5, 1), tolerance = 1e-10)
  expect_true(is.finite(res$estimates$std_residual[5]))
  # Problems in such rows are reported with their row and column.
  bad_a <- data.frame(disease = c("d1", "d1", "d1", "d1;d2", "d1"),
                      associated_with = c("d2", "d9", "d2", NA, "d1"),
                      among = c("d1", NA, "d3;no d3", NA, NA),
                      measure = c("OR", "OR", "OR", "prevalence", "OR"), value = c(2, 2, 2, 0.1, 2))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = bad_a, impacts = tabs$yield)
  expect_true(has_problem(pr, "associations", 1, "among"))            # about d1 and among d1
  expect_true(has_problem(pr, "associations", 2, "associated_with"))  # unknown disease
  expect_true(has_problem(pr, "associations", 3, "among"))            # d3 twice
  expect_false(has_problem(pr, "associations", 4))                    # a combination's prevalence
  expect_true(has_problem(pr, "associations", 5, "disease"))          # d1 on both sides
  bad_i <- data.frame(disease = c(ids3, "d1;d9", "d1"), among = c(NA, NA, NA, NA, "d1"), value = 1:5)
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad_i)
  expect_true(has_problem(pr, "impacts", 4, "disease"))
  expect_true(has_problem(pr, "impacts", 5, "among"))
  # A whole-population prevalence of one disease belongs in the diseases table.
  pr <- cm_check_inputs(diseases = tabs$diseases, impacts = tabs$yield,
                        associations = data.frame(disease = "d1", measure = "prevalence", value = 0.1))
  expect_match(pr$problem[pr$table == "associations"], "belongs in the diseases table")
})

test_that("the shipped global dairy files match example_global_dairy()", {
  f <- function(x) extdata("global_dairy_2024", x)
  expect_true(file.exists(f("culling.csv")))
  files <- c(yield = "impacts_yield.csv", fertility = "impacts_fertility.csv", culling = "culling.csv")
  for (outcome in names(files)) {
    expect_equal(nrow(cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                                      impacts = f(files[[outcome]]))), 0, info = outcome)
    m <- cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                        impacts = f(files[[outcome]]))
    ref <- example_global_dairy(outcome)
    expect_equal(m$diseases$disease, ref$diseases$disease)
    expect_equal(m$diseases$measure, ref$diseases$measure)
    expect_equal(m$diseases$prob, ref$diseases$prob, tolerance = 1e-12)
    km <- pair_key(m$associations$disease, m$associations$associated_with)
    kr <- pair_key(ref$associations$disease, ref$associations$associated_with)
    expect_setequal(km, kr)
    expect_equal(m$associations$value[match(kr, km)], ref$associations$value, tolerance = 1e-12)
    expect_equal(m$impacts$disease, ref$impacts$disease)
    expect_equal(m$impacts$value, ref$impacts$value, tolerance = 1e-12)
    expect_equal(m$impacts$measure, ref$impacts$measure)
    expect_equal(attr(m$impacts, "label"), attr(ref$impacts, "label"))
    expect_equal(attr(m$impacts, "units"), attr(ref$impacts, "units"))
    expect_equal(impact_kind(m$impacts), if (outcome == "culling") "event" else "additive")
    # The distributions in the tables are those of the 2024 analysis (the
    # rows are in the same order in the files and in R).
    dm <- m$distributions
    dr <- ref$distributions
    expect_setequal(names(dm), names(dr))
    expect_equal(vapply(dm[names(dr)], function(d) d$type, character(1)),
                 vapply(dr, function(d) d$type, character(1)))
    expect_equal(vapply(dm[names(dr)], function(d) d$mean, numeric(1)),
                 vapply(dr, function(d) d$mean, numeric(1)), tolerance = 1e-10)
  }
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
    table = c(rep("diseases", 5), rep("associations", 8), rep("impacts", 4)),
    row = c(2, 3, 4, 5, 6,
            2, 3, 4, 6, 7, 8, 9, 10,
            2, 3, 5, NA),
    column = c("dist", "value", "dist", "disease", "disease",
               "dist", "associated_with", "value", "measure", "measure", "among", NA, "among",
               "value", "dist", "disease", "units"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(expected))) {
    expect_true(has_problem(pr, expected$table[i], expected$row[i], expected$column[i]),
                info = paste(expected$table[i], expected$row[i], expected$column[i]))
  }
  # Each error is reported once.
  expect_equal(sum(pr$severity == "error"), nrow(expected))
  # The notes: an adjustment for factors outside the model, and a point
  # value outside its distribution (MET, uniform from 6 to 8).
  expect_equal(sum(pr$severity == "note"), 2)
  expect_true(has_problem(pr, "associations", 5, "adjusted_for", severity = "note"))
  expect_true(has_problem(pr, "impacts", 4, "dist", severity = "note"))
  # The decimal comma is a value problem, not a missing value as well.
  expect_equal(sum(pr$table == "impacts" & pr$row %in% 2), 1)
  # Retired measures, with what to do instead.
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 6 & pr$column %in% "measure"],
               "now 'prevalence'")
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 7 & pr$column %in% "measure"],
               "odds ratio of 1")
  expect_match(pr$problem[pr$table == "associations" & pr$row %in% 4], "Leave the pair out")
  expect_match(pr$problem[pr$table == "diseases" & pr$row %in% 2], "p2 \\(sd\\)")
  expect_match(pr$problem[pr$table == "diseases" & pr$row %in% 6], "belongs in the associations table")
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
  expect_output(print(pr), "Found 17 problem")
  expect_error(cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                              impacts = f("impacts_yield.csv")),
               "diseases, row 3, column 'value'", class = "deconflate_input_problems")

  # culling.csv: event impacts with three mistakes, and the estimand column
  # of earlier versions (a note).
  pc <- cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                        impacts = f("culling.csv"))
  expect_true(has_problem(pc, "impacts", 2, "measure"))        # no measure
  expect_true(has_problem(pc, "impacts", 3, "value"))          # HR of 0
  expect_true(has_problem(pc, "impacts", 4, "value"))          # RD outside (-1, 1)
  expect_true(has_problem(pc, "impacts", NA, "estimand", severity = "note"))
  expect_equal(sum(pc$severity == "error" & pc$table == "impacts"), 3)
})

test_that("files in the layout of earlier versions say what to change", {
  dis <- data.frame(id = c("a", "b"), value = c(0.1, 0.2), type = "prevalence")
  pr <- cm_check_inputs(diseases = dis, impacts = data.frame(disease = c("a", "b"), value = 1:2))
  expect_true(has_problem(pr, "diseases", NA, "id"))
  expect_true(has_problem(pr, "diseases", NA, "type"))
  expect_match(pr$problem[pr$column %in% "id"], "now called `disease`")
  dis2 <- data.frame(disease = c("a", "b"), value = c(0.1, 0.2))
  old_a <- data.frame(disease1 = "a", disease2 = "b", value = 2, adjusted = FALSE)
  pr <- cm_check_inputs(diseases = dis2, associations = old_a,
                        impacts = data.frame(disease = c("a", "b"), value = 1:2))
  expect_true(has_problem(pr, "associations", NA, "disease1"))
  expect_true(has_problem(pr, "associations", NA, "disease2"))
  expect_match(pr$problem[pr$column %in% "disease2"], "associated_with")
  expect_true(has_problem(pr, "associations", NA, "adjusted", severity = "note"))
  # A three-way file given as the associations.
  tw <- data.frame(disease1 = "a", disease2 = "b", disease3 = "c", ratio = 2)
  pr <- cm_check_inputs(diseases = dis2, associations = tw,
                        impacts = data.frame(disease = c("a", "b"), value = 1:2))
  expect_match(pr$problem[pr$table == "associations"][1], "three-way file")
  # An impact table with an outcome column.
  old <- data.frame(disease = c("a", "b", "a", "b"), outcome = c("yield", "yield", "culling", "culling"),
                    value = c(2, 3, 1.5, 1.8))
  pr <- cm_check_inputs(diseases = dis2, impacts = old)
  expect_true(has_problem(pr, "impacts", NA, "outcome"))
  expect_error(cm_read_inputs(diseases = dis2, impacts = old), class = "deconflate_input_problems")
  # scale and direction are dropped with a note only.
  pr <- cm_check_inputs(diseases = dis2,
                        impacts = data.frame(disease = c("a", "b"), value = c(2, 3),
                                             scale = "percent", direction = "decrease"))
  expect_true(has_problem(pr, "impacts", NA, NA, severity = "note"))
  expect_false(any(pr$severity == "error"))
})

test_that("all problems in tables typed in R are reported at once", {
  pr <- cm_check_inputs(
    diseases = data.frame(disease = c("a", "b", "b", "c"), value = c(0.1, 1.5, 0.2, 0.3)),
    associations = data.frame(disease = "a", associated_with = "z", value = -1),
    impacts = data.frame(disease = c("a", "b"), value = c("2", "x"), extra = 1,
                         dist = c("normal", NA), p1 = c(2, 1))
  )
  expect_s3_class(pr, "cm_problems")
  expect_true(has_problem(pr, "diseases", 2, "value"))
  expect_true(has_problem(pr, "diseases", 3, "disease"))
  expect_true(has_problem(pr, "associations", 1, "value"))
  expect_true(has_problem(pr, "impacts", 2, "value"))
  expect_true(has_problem(pr, "impacts", NA, NA))              # no impact for c
  expect_true(has_problem(pr, "impacts", NA, "extra", severity = "note"))
  expect_true(has_problem(pr, "impacts", 1, "dist"))           # normal without sd
  expect_true(has_problem(pr, "impacts", 2, "dist"))           # p1 without a distribution
  expect_error(cm_read_inputs(diseases = data.frame(disease = "a", value = 2),
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
  # Event impacts without a row for d3.
  pr <- cm_check_inputs(diseases = tabs$diseases,
                        impacts = data.frame(disease = c("d1", "d2"), value = c(1.5, 2), measure = "HR"))
  expect_true(has_problem(pr, "impacts", NA, NA))
  expect_match(pr$problem[pr$table == "impacts" & is.na(pr$row)], "1 for ratios")
})

test_that("retired association measures and count columns are errors", {
  tabs <- supp_tables()
  assoc <- data.frame(disease = c("d1", "d1", "d2", "d1"), associated_with = c("d2", "d3", "d3", "d2"),
                      value = c(NA, 1, NA, 0.4), measure = c("table", "independent", "unknown", "cond_prob"),
                      n11 = c(10, NA, NA, NA), n10 = c(20, NA, NA, NA), n01 = c(30, NA, NA, NA),
                      n00 = c(40, NA, NA, NA))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield)
  expect_true(has_problem(pr, "associations", NA, NA))         # the count columns, once
  expect_equal(sum(pr$table == "associations" & is.na(pr$row)), 1)
  for (r in 1:4) expect_true(has_problem(pr, "associations", r, "measure"), info = r)
  # A retired row is not also reported as a missing value.
  expect_false(has_problem(pr, "associations", 1, "value"))
  expect_false(has_problem(pr, "associations", 3, "value"))
  expect_equal(sum(pr$severity == "error"), 5)
  expect_error(cm_read_inputs(diseases = tabs$diseases, associations = assoc, impacts = tabs$yield),
               "no longer", class = "deconflate_input_problems")
  # Count columns alone (with an odds ratio in `value`) are still an error.
  assoc2 <- data.frame(disease = "d1", associated_with = "d2", value = 2, n11 = 10, n10 = 20, n01 = 30, n00 = 40)
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
  # An odds ratio of 1 states independence.
  pt <- pair_tables(supp_population())
  expect_equal(pt$p11[pt$disease1 == "d1" & pt$disease2 == "d3"], 0.10 * 0.20)
  expect_equal(pt$status, rep("specified", 3))
})

test_that("event impact tables have a measure column and are checked", {
  tabs <- supp_tables()
  ev <- data.frame(disease = ids3, value = c(1.5, 1.3, 0.04), measure = c("HR", "RR", "RD"), label = "culling")
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
                    adjusted_for = c(NA, NA, "d9"))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad)
  expect_true(has_problem(pr, "impacts", 2, "measure"))        # a measure in every row
  expect_true(has_problem(pr, "impacts", 3, "value"))          # RD outside (-1, 1)
  expect_true(has_problem(pr, "impacts", 3, "adjusted_for", severity = "note"))   # not a disease
  bad2 <- data.frame(disease = ids3, value = c(0, 1.3, 1.2), measure = c("HR", "HZ", "OR"))
  pr <- cm_check_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = bad2)
  expect_true(has_problem(pr, "impacts", 1, "value"))          # a ratio must be positive
  expect_true(has_problem(pr, "impacts", 2, "measure"))        # unknown measure
  # The same mistakes made in R.
  expect_error(cm_impacts(ids3, c(1.5, 0, 1.3), measure = "HR"), "must be positive")
  expect_error(cm_impacts(ids3, c(1.5, 2, 1.3), measure = "HZ"), "Invalid `measure`")
})

test_that("several rows per item are read and pooled", {
  tabs <- supp_tables()
  assoc <- rbind(tabs$associations, data.frame(disease = "d2", associated_with = "d1", value = 2.4))
  assoc$dist <- c("lognormal_ci", NA, NA, "lognormal_ci")
  assoc$p1 <- c(2, NA, NA, 2.4)
  assoc$p2 <- c(1.4, NA, NA, 1.5)
  assoc$p3 <- c(2.9, NA, NA, 3.8)
  imp <- data.frame(disease = c(ids3, "d1"), value = c(2.5, 5, 7.5, 2.1), adjusted_for = c(NA, NA, NA, "d2"),
                    dist = c("normal", NA, NA, "normal"), p1 = c(2.5, NA, NA, 2.1), p2 = c(0.5, NA, NA, 0.4))
  m <- cm_read_inputs(diseases = tabs$diseases, associations = assoc, impacts = imp)
  expect_equal(nrow(attr(m, "problems")), 0)
  expect_equal(names(m$distributions), c("associations:1", "associations:4", "impacts:1", "impacts:4"))
  pt <- pair_tables(m)
  expect_equal(pt$measure[1], "pooled")
  expect_equal(pt$rows[1], "1, 4")
  res <- deconflate(m, n_draws = 0)
  h <- res$heterogeneity
  expect_equal(h$table, c("associations", "impacts"))
  expect_equal(h$k, c(2, 4))
  expect_true(all(c("Q", "df", "p", "I2", "tau2") %in% names(h)))
  expect_equal(res$adjusted$rows, c(2, 1, 1))
  expect_equal(res$totals$raw_total, 0.10 * (2.5 / 0.5^2 + 2.1 / 0.4^2) / (1 / 0.5^2 + 1 / 0.4^2) + 0.15 * 5 + 0.20 * 7.5,
               tolerance = 1e-10)
  # Exact rows that disagree are an error.
  imp2 <- data.frame(disease = c(ids3, "d1"), value = c(2.5, 5, 7.5, 2.1))
  m2 <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = imp2)
  expect_error(deconflate(m2, n_draws = 0), "cannot all hold", class = "deconflate_infeasible")
})

test_that("distribution keys given to cm_model() name the rows they refer to", {
  dd <- list("impact:d1" = dist_normal(2.5, 0.5), "assoc:d2:d1" = dist_lognormal_ci(2, 1.4, 2.9))
  m <- cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5)), distributions = dd)
  expect_equal(names(m$distributions), c("impacts:1", "associations:1"))
})

test_that("CSV columns are read as text, and empty files are empty tables", {
  dir <- input_dir("text-columns")
  writeLines(c("disease,value", "T,0.10", "F,0.15", "01,0.20"), file.path(dir, "diseases.csv"))
  file.create(file.path(dir, "associations.csv"))
  writeLines(c("disease,value", "T,1", "F,2", "01,3"), file.path(dir, "losses.csv"))
  m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"), impacts = file.path(dir, "losses.csv"),
                      associations = file.path(dir, "associations.csv"))
  expect_equal(m$diseases$disease, c("T", "F", "01"))
  expect_null(m$associations)
  expect_equal(m$impacts$value, c(1, 2, 3))
  unlink(dir, recursive = TRUE)
})

test_that("mixed time horizons give a note", {
  d <- data.frame(disease = c("a", "b"), value = c(0.1, 0.2), time_horizon = c("year", "lactation"))
  pr <- cm_check_inputs(diseases = d, impacts = data.frame(disease = c("a", "b"), value = c(1, 2)))
  expect_true(has_problem(pr, "diseases", NA, "time_horizon", severity = "note"))
  expect_false(any(pr$severity == "error"))
})

# ---- Uncertainty in the tables ---------------------------------------------------

test_that("distributions are given in the rows of the values they describe", {
  tabs <- supp_tables()
  dis <- cbind(tabs$diseases, dist = c("beta", NA, NA), p1 = c(10, NA, NA), p2 = c(90, NA, NA))
  assoc <- rbind(cbind(tabs$associations, among = NA, measure = "OR"),
                 data.frame(disease = "d3", associated_with = "d1", value = 1.5, among = "d2", measure = "OR_ratio"))
  assoc <- cbind(assoc, dist = c("lognormal_ci", NA, "pert", "lognormal_ci"),
                 p1 = c(2, NA, 2, 1.5), p2 = c(1.4, NA, 3, 0.8), p3 = c(2.9, NA, 4.5, 2.8))
  imp <- data.frame(disease = c(ids3, "d1;d2"), value = c(2.5, 5, 7.5, 0.5), adjusted_for = c(NA, NA, NA, "all"),
                    units = "%", dist = c("normal", "pert", "fixed", "normal"),
                    p1 = c(2.5, 3, 7.5, 0.5), p2 = c(0.5, 5, NA, 0.1), p3 = c(NA, 8, NA, NA))
  pr <- cm_check_inputs(diseases = dis, associations = assoc, impacts = imp)
  expect_equal(nrow(pr), 0)
  m <- cm_read_inputs(diseases = dis, associations = assoc, impacts = imp)
  # Point values are the central inputs.
  expect_equal(m$diseases$value, c(0.10, 0.15, 0.20))
  expect_equal(m$associations$value[4], 1.5)
  expect_equal(m$impacts$value, c(2.5, 5, 7.5, 0.5))
  d <- m$distributions
  expect_equal(names(d), c("prob:d1", "associations:1", "associations:3", "associations:4",
                           "impacts:1", "impacts:2", "impacts:3", "impacts:4"))
  expect_equal(d[["prob:d1"]]$type, "beta")
  expect_equal(d[["associations:1"]]$type, "lognormal")
  expect_equal(d[["associations:1"]]$params$meanlog, log(2))
  expect_equal(d[["associations:4"]]$params$meanlog, log(1.5))
  expect_equal(d[["impacts:2"]]$type, "pert")
  expect_equal(d[["impacts:3"]]$type, "fixed")
  # The draws use them (the full model, because of the interaction).
  res <- deconflate(m, n_draws = 20, seed = 1)
  expect_equal(res$method, "full model")
  expect_equal(colnames(res$draws$params), names(d))
  expect_gt(length(unique(res$draws$params[, "associations:4"])), 1)
  expect_true(all(res$draws$params[, "impacts:3"] == 7.5))
  expect_equal(nrow(res$draws$params), 20 - res$draws$n_rejected)

  # Event impacts can have distributions too.
  ev <- data.frame(disease = ids3, value = c(1.5, 2, 1.3), measure = "HR",
                   dist = c("lognormal_ci", NA, NA), p1 = c(1.5, NA, NA), p2 = c(1.2, NA, NA),
                   p3 = c(1.9, NA, NA))
  me <- cm_read_inputs(diseases = tabs$diseases, associations = tabs$associations, impacts = ev)
  expect_equal(names(me$distributions), "impacts:1")
  expect_equal(nrow(attr(me, "problems")), 0)
})

test_that("problems with distributions are reported per row", {
  dis <- data.frame(disease = c("a", "b", "c"), value = c(0.1, 0.2, 0.3),
                    dist = c(NA, "uniform", NA), p1 = c(0.05, 0.25, NA), p2 = c(NA, 0.4, NA))
  assoc <- data.frame(disease = c("a", "b"), associated_with = c("b", "c"), value = c(2, 2),
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
  f <- function(x) extdata("five_diseases", x)
  expect_setequal(list.files(extdata("five_diseases"), pattern = "\\.csv$"),
                  c("diseases.csv", "associations.csv", "impacts_yield.csv", "impacts_calving_interval.csv",
                    "impacts_welfare.csv", "culling.csv"))
  for (imp in c("impacts_yield.csv", "impacts_calving_interval.csv", "impacts_welfare.csv", "culling.csv")) {
    pr <- cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                          impacts = f(imp))
    expect_equal(nrow(pr), 0, info = imp)
  }
  ids <- c("LAM", "MAS", "MET", "SCK", "RP")
  y <- read_five("impacts_yield.csv")
  expect_equal(y$diseases$disease, ids)
  expect_equal(y$diseases$prob, c(0.25, 1 - exp(-0.3), 0.10, 0.35, 0.06), tolerance = 1e-12)
  expect_equal(y$diseases$label[1], "lameness")
  pt <- pair_tables(y)
  expect_equal(nrow(pt), 10)
  expect_true(all(pt$status == "specified"))
  expect_setequal(unique(pt$measure), c("pooled", "OR", "RR", "RD", "prevalence", "phi"))
  expect_equal(attr(y$impacts, "units"), "% of yield")

  # Distributions from every table, mixed with point values.
  expect_setequal(names(y$distributions),
                  c("prob:LAM", "prob:MAS", "prob:MET", "prob:SCK", paste0("associations:", c(1:8, 12)),
                    paste0("impacts:", 1:4)))
  types <- vapply(y$distributions, function(d) d$type, character(1))
  expect_true(all(c("beta", "pert", "uniform", "lognormal", "normal", "fixed") %in% types))

  # The ten pair rows alone (the example of earlier versions): the pairwise
  # shortcut and the references of reference_five_diseases.py.
  y10 <- read_five_pairs("impacts_yield.csv")
  r10 <- deconflate(y10, n_draws = 0)
  expect_equal(r10$method, "pairwise shortcut")
  expect_equal(r10$adjusted$adjusted,
               c(4.13927424, 2.35580211, 4.5140575, 1.36629496, 2.66542027), tolerance = 1e-7)
  expect_equal(r10$totals$adjusted_total, 2.7349337441102657, tolerance = 1e-9)

  # The full example (reference_v041.py): two pooled LAM-MAS rows and an
  # OR_ratio row, so the full model.
  r <- deconflate(y, n_draws = 0)
  expect_equal(r$method, "full model")
  expect_equal(r$adjusted$adjusted, c(4.1036014117, 2.2841850068, 4.5202762922, 1.380768132, 2.668536736),
               tolerance = 1e-7)
  expect_equal(r$totals$adjusted_total, 2.7133281669, tolerance = 1e-8)
  h <- r$heterogeneity
  expect_equal(h$pooled[h$table == "associations"], 1.9609933673, tolerance = 1e-8)
  expect_equal(h$Q[h$table == "associations"], 0.617707938, tolerance = 1e-8)
  j <- r$joint
  expect_equal(sum(j$prob[j$cells[, "LAM"] == 1 & j$cells[, "MAS"] == 1 & j$cells[, "SCK"] == 1]),
               0.055709677209, tolerance = 1e-8)
  ar <- r$estimates[r$estimates$table == "associations" & r$estimates$row == 12, ]
  expect_equal(ar$fitted, 1.5, tolerance = 1e-6)
  # Calving interval: two MAS rows (weighted least squares), SCK adjusted for
  # all (exact), RP only against healthy cows.
  ci <- read_five("impacts_calving_interval.csv")
  rc <- deconflate(ci, n_draws = 0)
  expect_equal(rc$adjusted$adjusted, c(10.2063061612, 4.7682049994, 14.7741251171, 4, 10.0020415569),
               tolerance = 1e-7)
  expect_equal(rc$totals$adjusted_total, 7.2649434013, tolerance = 1e-8)
  ei <- rc$estimates[rc$estimates$table == "impacts", ]
  expect_equal(ei$std_residual, c(-0.0090983656, -0.3092644373, 0.4653633829, -0.0039776718, NA,
                                  -0.000808178), tolerance = 1e-6)
  expect_equal(rc$heterogeneity$Q[rc$heterogeneity$table == "impacts"], 0.3123068256, tolerance = 1e-6)
  # Welfare: two interactions adjusted for all, and a crude combination row.
  w <- read_five("impacts_welfare.csv")
  rw <- deconflate(w, n_draws = 0)
  expect_equal(rw$adjusted$disease, c(ids, "LAM;MAS", "MET;RP", "LAM;SCK"))
  expect_equal(rw$adjusted$adjusted, c(8.7442976935, 3.8319918047, 3.1977586251, 1.0855585373, 2.1515488417,
                                       2, 1.5, -0.6200833573), tolerance = 1e-7)
  expect_equal(rw$totals$adjusted_total, 4.1479565223, tolerance = 1e-8)
  expect_equal(rw$contributions$total, c(2.2400759141, 1.0837970506, 0.3357170535, 0.3433323825, 0.1450341215),
               tolerance = 1e-8)

  # Culling: event impacts with mixed measures, crude and stratified rows
  # and a combination (reference_v041.py).
  cu <- read_five("culling.csv")
  expect_equal(cu$impacts$measure, c("HR", "HR", "RR", "OR", "RD", "HR"))
  expect_equal(cu$impacts$adjusted_for[c(2, 5)], c("LAM;SCK", "all"))
  expect_error(deconflate(cu, n_draws = 0), class = "deconflate_unsupported")
  ev <- deconflate(cu, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_s3_class(ev, "cm_event_result")
  expect_equal(ev$adjusted$adjusted, c(1.060298968, 1.2342866774, 1.080381173, 1.1753668593, 1.509273874,
                                       2.2779452712), tolerance = 1e-7)
  expect_equal(ev$attributable$summary$attributable, 0.0541469811, tolerance = 1e-8)

  # Thresholds used in run_all_features.R (the ten pair rows).
  th <- cm_threshold(y10, "impacts:4", c(0, 2.5), conclusion = "sign")
  t1 <- th$thresholds[th$thresholds$status == "threshold", ]
  expect_equal(nrow(t1), 1)
  expect_equal(t1$threshold, 1.1935832731450715, tolerance = 1e-6)
  th <- cm_threshold(y10, "impacts:3", c(0, 15), conclusion = "total", target = 3)
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 9.730627246165994,
               tolerance = 1e-6)
  th <- cm_threshold(y10, "prob:SCK", c(0.1, 0.6), conclusion = "rank", diseases = c("MET", "SCK"))
  expect_equal(th$thresholds$threshold[th$thresholds$status == "threshold"], 0.33002218486278845,
               tolerance = 1e-6)
  skip_on_cran()
  th <- cm_threshold(w, "interaction:LAM;MAS", c(-5, 20), conclusion = "change", target = -0.1)
  expect_equal(sum(th$thresholds$status == "threshold"), 1)
  # Draws of every input: almost every draw is usable.
  dr <- deconflate(w, n_draws = 30, seed = 1)$draws
  expect_lte(dr$n_rejected, 3)
  expect_equal(colnames(dr$params), names(w$distributions))
})
