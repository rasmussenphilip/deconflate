#' Read model inputs from CSV files or data frames
#'
#' Builds a model ([cm_model()]) from three tables: the diseases, their
#' associations and one impact table. Each table can be a path to a CSV file
#' (with any name) or a data frame typed in R. Values with a distribution
#' (see Uncertainty) are kept in the model, and [deconflate()] uses them to
#' compute intervals. [cm_template()] writes an example set to start from.
#'
#' All tables are checked before anything is built, and every problem is
#' reported at once, with its table, row and column (row 1 is the first row
#' below the header). [cm_check_inputs()] runs the same checks without
#' stopping.
#'
#' For several outcomes (e.g. milk yield and calving interval), keep one
#' impact table per outcome and read and adjust each in turn: the disease and
#' association tables can be shared.
#'
#' @section Columns:
#' Column names are not case-sensitive. Optional columns can be left out or
#' left empty; unrecognised columns are ignored (with a note). Every table
#' can have a `note` column for free text (e.g. the source of a value),
#' which the package never uses. In every
#' table, `disease` (and `associated_with`) take a disease id or a
#' combination of ids separated by `;` (e.g. `MAS;SCK`), and `among` takes a
#' subgroup (`SCK`, `no SCK`, `SCK;no MET`, or `no others`).
#'
#' **diseases** (one row per disease; whole-population values only):
#' `disease`, `value` (required); `measure` (`prevalence` (default),
#' `probability` or `incidence_rate`, converted with `1 - exp(-value)`),
#' `note`. Ids must not contain `|`, `;` or `:` or start with `no `; `all`
#' and `others` are reserved. Prevalences, probabilities, associations,
#' impacts and the overall risk must all refer to the same period (e.g. a
#' lactation or a year) and population.
#'
#' **associations** (one row per estimate; any row about two or more
#' diseases): `disease`, `value` (required); `associated_with`, `among`,
#' `adjusted_for`, `measure` (`OR` (default), `RR`, `RD`, `phi`,
#' `prevalence` or `OR_ratio`), `study` (identifies the row in `$estimates`
#' and the diagnosis), `note`. See
#' [cm_associations()] for what each row means. Pairs that no row mentions
#' are unknown: the model fills in their association. An odds ratio of 1
#' states that two diseases are unrelated. Several rows about the same
#' association are pooled.
#'
#' **impacts** (one row per estimate): `disease`, `value` (required);
#' `among`, `adjusted_for`, `measure`, `study`, `note`. Every disease needs a row of
#' its own; a row for a combination (e.g. `SCK;SCM`) adds an interaction to
#' the model, and a combination adjusted for `all` is the interaction itself.
#' See [cm_impacts()]. Without a `measure` column the impacts are additive,
#' in any units (the same for every row; results come back in those units).
#' With a `measure` column they are event impacts: every row needs a measure
#' (`HR`, `rate_ratio`, `RR`, `OR` or `RD`), and they are adjusted with
#' `deconflate(..., event_model = TRUE)`. What the impacts are about and
#' their units are not columns: give them as `outcome` and `units` (see
#' Arguments).
#'
#' Files in the layout of earlier versions (columns `id` and `type`,
#' `disease1` and `disease2`, `estimand`; three-way and interaction files)
#' are reported with what to change. The columns `label` (of the diseases
#' and associations tables), `time_horizon`, `reference_population` and
#' `source` are no longer read: files with them still load, with a note.
#' The columns `outcome`, `label` and `units` of an impact table are no
#' longer read either (a note says which arguments to use instead), but a
#' table whose rows have different outcomes or units is an error.
#'
#' @section Uncertainty:
#' Every table can have the columns `dist` and `p1`-`p4`. A row with a
#' `dist` has an uncertain value: [deconflate()] draws it from its
#' distribution (`n_draws`), while `value` gives the central estimate. A row
#' with an empty `dist` is a fixed point value (and, when a disease or
#' association has several rows, an exact one), so point values and
#' uncertain values can be mixed freely. Values are drawn on the scale they
#' were entered on (e.g. an incidence rate, an odds ratio, a hazard ratio or
#' a percent impact). A point value outside the support of its distribution
#' is noted. Distributions and their parameters:
#'
#' | `dist` | `p1` | `p2` | `p3` | `p4` |
#' |---|---|---|---|---|
#' | `fixed` | value | | | |
#' | `normal` | mean | sd | lower bound (optional) | upper bound (optional) |
#' | `lognormal` | meanlog | sdlog | | |
#' | `lognormal_ci` | estimate | lower CI | upper CI | level (default 0.95) |
#' | `beta` | shape1 | shape2 | min (default 0) | max (default 1) |
#' | `pert` | min | mode | max | lambda (default 4) |
#' | `pert_mean` | min | mean | max | lambda (default 4) |
#' | `uniform` | min | max | | |
#'
#' @param diseases,associations,impacts Paths to CSV files (any names) or
#'   data frames. `diseases` and `impacts` are required; `associations` can
#'   be left out only for the sensitivity tools ([deconflate()] needs at
#'   least one association).
#' @param ... Not used (arguments of earlier versions give a message). The
#'   arguments after it must be named.
#' @param outcome,units Optional: what the impacts are about and their units,
#'   one text each (e.g. `outcome = "milk yield loss"`, `units = "% of
#'   yield"`). They are used in the printouts and the plot titles and axes;
#'   without them, these use general wording and the axes show no unit. See
#'   [cm_impacts()].
#' @return A [cm_model()], with the distributions of uncertain values in
#'   `$distributions` (keyed by table and row, e.g. `"associations:3"`) and
#'   any notes from the checks in attribute `"problems"`. When a check finds
#'   an error, the error (class `deconflate_input_problems`) holds every
#'   problem found in `$problems`, as [cm_check_inputs()] returns them. The notes are also
#'   shown as a message of class `deconflate_input_notes`, which holds them
#'   in `$problems`; `suppressMessages()` hides it.
#' @export
#' @examples
#' dir <- file.path(tempdir(), "deconflate-inputs")
#' cm_template(dir, overwrite = TRUE)
#' m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
#'                     associations = file.path(dir, "associations.csv"),
#'                     impacts = file.path(dir, "yield.csv"),
#'                     outcome = "milk yield loss", units = "% of yield")
#' m
#' deconflate(m, n_draws = 0)
#'
#' # Example files shipped with the package
#' ex <- system.file("extdata", "five_diseases", package = "deconflate")
#' five <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
#'                        associations = file.path(ex, "associations.csv"),
#'                        impacts = file.path(ex, "impacts_yield.csv"),
#'                        units = "% of yield")
#' five
#' err <- system.file("extdata", "example_with_errors", package = "deconflate")
#' cm_check_inputs(diseases = file.path(err, "diseases.csv"),
#'                 associations = file.path(err, "associations.csv"),
#'                 impacts = file.path(err, "impacts_yield.csv"))
#'
#' # The same kind of tables typed in R
#' m2 <- cm_read_inputs(
#'   diseases = data.frame(disease = c("d1", "d2"), value = c(0.10, 0.15)),
#'   associations = data.frame(disease = "d1", associated_with = "d2", value = 2),
#'   impacts = data.frame(disease = c("d1", "d2"), value = c(2.5, 5)),
#'   units = "% of yield"
#' )
#' deconflate(m2)
cm_read_inputs <- function(diseases, associations = NULL, impacts, ..., outcome = NULL, units = NULL) {
  old_args(list(...), c(three_way = three_way_message, interactions = interactions_message,
                        adjusted_associations = "`adjusted_associations` is no longer used: an association adjusted for diseases of the model is used exactly, and one adjusted only for other factors is used as crude (with a note)."))
  if (missing(diseases)) cm_abort("`diseases` is required: a path to a CSV file or a data frame.", class = "deconflate_invalid_input")
  if (missing(impacts)) cm_abort("`impacts` is required: a path to a CSV file or a data frame.", class = "deconflate_invalid_input")
  outcome <- meta_text(outcome, "outcome")
  units <- meta_text(units, "units")
  chk <- run_input_checks(diseases, associations, impacts, outcome, units)
  pr <- chk$problems
  if (any(pr$severity == "error")) {
    cm_abort(format_problems(pr[pr$severity == "error", , drop = FALSE]),
             class = c("deconflate_input_problems", "deconflate_invalid_input"), problems = pr)
  }
  if (nrow(pr)) cm_inform(format_problems(pr), class = "deconflate_input_notes", problems = pr)
  m <- chk$model
  attr(m, "problems") <- pr
  m
}

#' Check model input tables
#'
#' Runs the checks of [cm_read_inputs()] and returns every problem found,
#' without stopping.
#'
#' @inheritParams cm_read_inputs
#' @return A data frame (class `cm_problems`) with `table`, `row`, `column`,
#'   `severity` (`"error"` or `"note"`) and `problem`. No rows means the
#'   inputs can be read.
#' @export
#' @examples
#' cm_check_inputs(
#'   diseases = data.frame(disease = c("d1", "d2"), value = c(0.10, 1.5)),
#'   impacts = data.frame(disease = c("d1", "d3"), value = c(2, "x"))
#' )
cm_check_inputs <- function(diseases = NULL, associations = NULL, impacts = NULL, ..., outcome = NULL,
                            units = NULL) {
  old_args(list(...), c(three_way = three_way_message, interactions = interactions_message,
                        adjusted_associations = "`adjusted_associations` is no longer used."))
  outcome <- meta_text(outcome, "outcome")
  units <- meta_text(units, "units")
  run_input_checks(diseases, associations, impacts, outcome, units)$problems
}

#' Write a template set of input files
#'
#' Writes an example set of input files to a folder, to edit and read back
#' with [cm_read_inputs()]. The values are illustrative. Some values have a
#' distribution in their row (columns `dist` and `p1`-`p4`); the others are
#' point values. The files are:
#'
#' * `diseases.csv`: three dairy diseases;
#' * `associations.csv`: their associations (odds ratios);
#' * `yield.csv`: additive impacts on milk yield, in percent of yield (any
#'   units can be used; give them with `cm_read_inputs(units = )`);
#' * `culling.csv`: event impacts on culling (hazard ratios and a risk ratio),
#'   for `deconflate(..., event_model = TRUE, overall_risk = ...)`.
#'
#' The folder `system.file("extdata", "five_diseases", package =
#' "deconflate")` has a fuller example that uses every feature (several
#' rows per item, subgroups, adjusted rows, combinations), and the folder
#' `pigs` an example outside dairy cattle (four respiratory diseases of
#' growing pigs, with three outcomes; see `vignette("other-settings")`).
#'
#' @param dir Folder to write to (created if needed).
#' @param overwrite Overwrite existing files?
#' @return The file paths, invisibly.
#' @export
#' @examples
#' dir <- file.path(tempdir(), "my-inputs")
#' cm_template(dir, overwrite = TRUE)
#' yield <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
#'                         associations = file.path(dir, "associations.csv"),
#'                         impacts = file.path(dir, "yield.csv"),
#'                         outcome = "milk yield loss", units = "% of yield")
#' deconflate(yield, n_draws = 200, seed = 1)
#' culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
#'                           associations = file.path(dir, "associations.csv"),
#'                           impacts = file.path(dir, "culling.csv"), outcome = "culling")
#' deconflate(culling, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
cm_template <- function(dir, overwrite = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  tabs <- list(
    diseases = data.frame(
      disease = c("LAM", "SCK", "MET"), measure = c("prevalence", "incidence_rate", "prevalence"),
      value = c(0.25, 0.48, 0.10),
      dist = c("beta", NA, NA), p1 = c(78.29, NA, NA), p2 = c(227.42, NA, NA), p3 = NA, p4 = NA,
      note = c("lameness; illustrative; beta(shape1, shape2)", "subclinical ketosis; cases per lactation",
               "metritis"),
      stringsAsFactors = FALSE),
    associations = data.frame(
      disease = c("LAM", "MET", "MET"), associated_with = c("SCK", "SCK", "LAM"),
      among = NA, adjusted_for = NA, measure = "OR", value = c(2.01, 1.94, 6.10), study = NA,
      dist = c("normal", NA, NA), p1 = c(2.01, NA, NA), p2 = c(0.20, NA, NA), p3 = c(0, NA, NA), p4 = NA,
      note = c("illustrative; normal(mean, sd) truncated at 0", NA, NA),
      stringsAsFactors = FALSE),
    yield = data.frame(
      disease = c("LAM", "SCK", "MET"), among = NA, adjusted_for = NA, value = c(4.81, 8.40, 5.61),
      study = NA,
      dist = c("normal", "normal", NA), p1 = c(4.81, 8.40, NA), p2 = c(0.87, 1.19, NA), p3 = NA, p4 = NA,
      note = c("milk yield loss in % of yield; illustrative; normal(mean, sd)", "normal(mean, sd)", NA),
      stringsAsFactors = FALSE),
    culling = data.frame(
      disease = c("LAM", "SCK", "MET"), among = NA, adjusted_for = NA, measure = c("HR", "HR", "RR"),
      value = c(1.74, 1.92, 1.45), study = NA,
      dist = c("lognormal_ci", NA, NA), p1 = c(1.74, NA, NA), p2 = c(1.45, NA, NA), p3 = c(2.09, NA, NA),
      p4 = NA,
      note = c("illustrative; hazard ratio; lognormal_ci(estimate, lower, upper)", "hazard ratio",
               "risk ratio of culling within the lactation"),
      stringsAsFactors = FALSE)
  )
  paths <- file.path(dir, paste0(names(tabs), ".csv"))
  exists <- file.exists(paths)
  if (any(exists) && !overwrite) {
    cm_abort(sprintf("Files already exist (use overwrite = TRUE): %s",
                     paste(basename(paths[exists]), collapse = ", ")), class = "deconflate_invalid_input")
  }
  for (k in seq_along(tabs)) utils::write.csv(tabs[[k]], paths[k], row.names = FALSE, na = "")
  cm_inform(sprintf("Wrote %d files to %s", length(paths), dir))
  invisible(paths)
}

#' @export
print.cm_problems <- function(x, ...) {
  if (!nrow(x)) {
    cat("No problems found.\n")
  } else {
    d <- x
    class(d) <- "data.frame"
    cat(format_problems(d), "\n")
  }
  invisible(x)
}

# Internals -------------------------------------------------------------------

# Columns that describe the uncertainty of a row's value.
dist_columns <- c("dist", "p1", "p2", "p3", "p4")

input_table_specs <- list(
  diseases = list(required = c("disease", "value"),
                  optional = c("measure", "note", dist_columns)),
  associations = list(required = c("disease", "value"),
                      optional = c("associated_with", "among", "adjusted_for", "measure", "study",
                                   "note", dist_columns)),
  impacts = list(required = c("disease", "value"),
                 optional = c("among", "adjusted_for", "measure", "study",
                              "note", dist_columns))
)

# Read one table from a path or a data frame, with normalised column names.
read_input_table <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) cm_abort(sprintf("File not found: %s", x), class = "deconflate_invalid_input")
    # An empty file is an empty table. Columns are read as text and converted
    # by the checks, so that ids such as "T" or "01" are not type-converted.
    if (!isTRUE(file.size(x) > 0)) return(data.frame())
    x <- utils::read.csv(x, stringsAsFactors = FALSE, na.strings = c("", "NA"),
                         strip.white = TRUE, check.names = FALSE,
                         colClasses = "character", fileEncoding = "UTF-8-BOM")
  }
  if (!is.data.frame(x)) cm_abort("Each table must be a data frame or the path to a CSV file.", class = "deconflate_invalid_input")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  names(x) <- gsub("[ .]+", "_", tolower(trimws(names(x))))
  x
}

clean_chr <- function(v, default = NA_character_, n = NULL) {
  if (is.null(v)) return(rep(default, n %||% 0L))
  v <- trimws(as.character(v))
  v[!is.na(v) & !nzchar(v)] <- NA_character_
  v[is.na(v)] <- default
  v
}

# Build one distribution from a `dist` name and parameters p1-p4.
input_dist <- function(dist, p) {
  p <- as.numeric(p)
  length(p) <- 4L
  need <- function(k, what) {
    if (any(is.na(p[seq_len(k)]))) cm_abort(sprintf("'%s' needs %s.", dist, what), class = "deconflate_invalid_input")
  }
  opt <- function(i, default) if (is.na(p[i])) default else p[i]
  if (is.na(dist)) cm_abort("Missing distribution name.", class = "deconflate_invalid_input")
  switch(dist,
    fixed = {
      need(1, "p1 (value)")
      dist_fixed(p[1])
    },
    normal = {
      need(2, "p1 (mean) and p2 (sd)")
      dist_normal(p[1], p[2], lower = opt(3, -Inf), upper = opt(4, Inf))
    },
    lognormal = {
      need(2, "p1 (meanlog) and p2 (sdlog)")
      dist_lognormal(p[1], p[2])
    },
    lognormal_ci = {
      need(3, "p1 (estimate), p2 (lower) and p3 (upper)")
      dist_lognormal_ci(p[1], p[2], p[3], level = opt(4, 0.95))
    },
    beta = {
      need(2, "p1 (shape1) and p2 (shape2)")
      dist_beta(p[1], p[2], min = opt(3, 0), max = opt(4, 1))
    },
    pert = {
      need(3, "p1 (min), p2 (mode) and p3 (max)")
      dist_pert(p[1], p[2], p[3], lambda = opt(4, 4))
    },
    pert_mean = {
      need(3, "p1 (min), p2 (mean) and p3 (max)")
      dist_pert_mean(p[1], p[2], p[3], lambda = opt(4, 4))
    },
    uniform = {
      need(2, "p1 (min) and p2 (max)")
      dist_uniform(p[1], p[2])
    },
    cm_abort(sprintf("Unknown distribution '%s' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).", dist), class = "deconflate_invalid_input")
  )
}

format_problems <- function(pr) {
  where <- ifelse(is.na(pr$row), pr$table, sprintf("%s, row %d", pr$table, pr$row))
  where <- ifelse(is.na(pr$column), where, sprintf("%s, column '%s'", where, pr$column))
  n_err <- sum(pr$severity == "error")
  head <- if (n_err) {
    sprintf("Found %d problem(s) in the inputs:", n_err)
  } else {
    "Notes on the inputs:"
  }
  tag <- ifelse(pr$severity == "note", " (note)", "")
  paste(c(head, sprintf("  %s%s: %s", where, tag, pr$problem)), collapse = "\n")
}

# Problem collector shared by the table checks.
new_problems <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e$add <- function(table, row, column, problem, severity = "error") {
    e$rows[[length(e$rows) + 1L]] <- data.frame(
      table = table, row = as.integer(row), column = as.character(column),
      severity = severity, problem = problem, stringsAsFactors = FALSE)
  }
  e$table <- function() {
    pr <- if (length(e$rows)) do.call(rbind, e$rows) else {
      data.frame(table = character(0), row = integer(0), column = character(0),
                 severity = character(0), problem = character(0), stringsAsFactors = FALSE)
    }
    rownames(pr) <- NULL
    class(pr) <- c("cm_problems", "data.frame")
    pr
  }
  e$has_errors <- function() any(vapply(e$rows, function(r) r$severity == "error", logical(1)))
  e
}

# Read a table; problems reading it are recorded. Empty tables become NULL.
read_checked <- function(src, label, pc) {
  tab <- tryCatch(read_input_table(src), error = function(e) e)
  if (inherits(tab, "condition")) {
    pc$add(label, NA, NA, conditionMessage(tab))
    return(NULL)
  }
  if (!is.null(tab) && !nrow(tab)) return(NULL)
  tab
}

has_columns <- function(tab, spec, label, pc) {
  sp <- input_table_specs[[spec]]
  nm <- names(tab)
  renamed <- list(
    diseases = c(id = "disease", type = "measure"),
    associations = c(disease1 = "disease", disease2 = "associated_with"))[[spec]]
  if (spec == "associations" && all(c("disease3", "ratio") %in% nm)) {
    pc$add(label, NA, NA, paste("This looks like a three-way file, which is no longer used.", three_way_message))
    return(FALSE)
  }
  if (spec == "impacts" && all(c("disease1", "disease2") %in% nm) && !"disease" %in% nm) {
    pc$add(label, NA, NA, paste("This looks like an interactions file, which is no longer used.", interactions_message))
    return(FALSE)
  }
  for (old in intersect(names(renamed), nm)) {
    pc$add(label, NA, old, sprintf("Column `%s` is now called `%s`: rename it.", old, renamed[[old]]))
  }
  miss <- setdiff(sp$required, nm)
  miss <- setdiff(miss, unname(renamed[names(renamed) %in% nm]))
  for (col in miss) pc$add(label, NA, col, "Required column is missing.")
  extra <- setdiff(nm, c(sp$required, sp$optional, names(renamed)))
  # The outcome and units of an impact table are arguments now
  # (check_impact_table() reports these columns).
  if (spec == "impacts") extra <- setdiff(extra, meta_columns)
  if (spec == "impacts" && "estimand" %in% extra) {
    pc$add(label, NA, "estimand", paste(sub("`estimand` is", "The estimand column is", estimand_message), "Column ignored."),
           "note")
    extra <- setdiff(extra, "estimand")
  }
  if (spec == "impacts" && any(c("scale", "direction") %in% extra)) {
    pc$add(label, NA, NA,
           "Columns scale and direction are no longer used: impacts are adjusted in their own units; convert the adjusted impacts as needed afterwards.",
           "note")
    extra <- setdiff(extra, c("scale", "direction"))
  }
  if (spec == "associations" && "adjusted" %in% extra) {
    pc$add(label, NA, "adjusted",
           "The adjusted column is no longer used (column ignored): an association is adjusted when adjusted_for names diseases of the model; one adjusted only for other factors is used as crude.",
           "note")
    extra <- setdiff(extra, "adjusted")
  }
  if (spec == "associations" && any(c("n11", "n10", "n01", "n00") %in% extra)) {
    pc$add(label, NA, NA,
           "Contingency-table counts (n11, n10, n01, n00) are no longer read: compute each table's odds ratio, n11 * n00 / (n10 * n01), and enter it in `value` with measure OR.")
    extra <- setdiff(extra, c("n11", "n10", "n01", "n00"))
  }
  for (col in intersect(extra, unused_inputs)) {
    pc$add(label, NA, col, sprintf("Column `%s` is not used; put free text in `note`.", col), "note")
  }
  extra <- setdiff(extra, unused_inputs)
  for (col in extra) pc$add(label, NA, col, "Unrecognised column: ignored (free text goes in `note`).", "note")
  !length(miss) && !length(intersect(names(renamed), nm))
}

# Distributions given in a table's dist and p1-p4 columns: a named list
# (names = `keys`) for the rows with a `dist`. Rows without one are point
# values. `values` are the point values, checked against each distribution's
# support; `allowed` marks rows that may have a distribution.
row_dists <- function(tab, label, pc, keys, values, allowed = NULL, allowed_msg = NULL) {
  if (!any(c("dist", paste0("p", 1:4)) %in% names(tab))) return(list())
  dist <- tolower(chr_col(tab, "dist"))
  P <- vapply(paste0("p", 1:4), function(col) num_col(tab, col, label, pc, FALSE), numeric(nrow(tab)))
  P <- matrix(P, nrow = nrow(tab))
  # Parameters that are not numbers are reported by num_col(); do not report
  # the same row again as a distribution with missing parameters.
  given <- vapply(paste0("p", 1:4), function(col) !is.na(chr_col(tab, col)), logical(nrow(tab)))
  given <- matrix(given, nrow = nrow(tab))
  bad_p <- rowSums(given & is.na(P)) > 0
  out <- list()
  for (r in seq_len(nrow(tab))) {
    if (bad_p[r]) next
    if (is.na(dist[r])) {
      if (any(!is.na(P[r, ]))) {
        pc$add(label, r, "dist", "p1-p4 are filled in but `dist` is empty: name the distribution in `dist`, or leave p1-p4 empty.")
      }
      next
    }
    if (!is.null(allowed) && !isTRUE(allowed[r])) {
      pc$add(label, r, "dist", allowed_msg)
      next
    }
    dd <- tryCatch(input_dist(dist[r], P[r, ]), error = function(e) e)
    if (inherits(dd, "condition")) {
      pc$add(label, r, "dist", conditionMessage(dd))
      next
    }
    if (is.na(keys[r])) next
    v <- values[r]
    if (is.finite(v) && !any(dist_support(dd)[, 1] <= v & dist_support(dd)[, 2] >= v)) {
      pc$add(label, r, "dist",
             sprintf("The point value %g lies outside the distribution's support %s.", v, format_support(dd)),
             "note")
    }
    out[[keys[r]]] <- dd
  }
  out
}

num_col <- function(tab, col, label, pc, required = TRUE) {
  v <- tab[[col]]
  n <- nrow(tab)
  if (is.null(v)) return(rep(NA_real_, n))
  bad <- rep(FALSE, n)
  if (is.numeric(v)) {
    x <- as.numeric(v)
    bad <- is.nan(x)
    for (r in which(bad)) pc$add(label, r, col, "NaN is not a number.")
  } else {
    vv <- trimws(as.character(v))
    x <- suppressWarnings(as.numeric(vv))
    bad <- !is.na(vv) & nzchar(vv) & is.na(x)
    for (r in which(bad)) pc$add(label, r, col, sprintf("'%s' is not a number.", vv[r]))
  }
  if (required) for (r in which(is.na(x) & !bad)) pc$add(label, r, col, "Missing value (this column is required).")
  x
}

chr_col <- function(tab, col, default = NA_character_) clean_chr(tab[[col]], default, nrow(tab))

check_disease_table <- function(d, pc) {
  if (!has_columns(d, "diseases", "diseases", pc)) return(NULL)
  id <- chr_col(d, "disease")
  for (r in which(is.na(id))) pc$add("diseases", r, "disease", "Missing disease id.")
  for (r in which(duplicated(id) & !is.na(id))) {
    pc$add("diseases", r, "disease",
           sprintf("Duplicate disease '%s': the diseases table has one row per disease (rows about several diseases or subgroups belong in the associations table).", id[r]))
  }
  for (r in which(!is.na(id) & grepl("[|;:]", id))) {
    pc$add("diseases", r, "disease",
           sprintf("Disease id '%s' must not contain '|', ';' or ':' (a combination of diseases belongs in the associations table).", id[r]))
  }
  for (r in which(!is.na(id) & tolower(id) %in% c("all", "others"))) {
    pc$add("diseases", r, "disease", "'all' and 'others' are reserved and cannot be disease ids.")
  }
  for (r in which(!is.na(id) & grepl("^no\\s", id, ignore.case = TRUE))) {
    pc$add("diseases", r, "disease", "Disease ids must not start with 'no ' (it marks absent diseases in `among`).")
  }
  value <- num_col(d, "value", "diseases", pc)
  for (r in which(!is.na(value) & value < 0)) pc$add("diseases", r, "value", "The value must not be negative.")
  measure <- chr_col(d, "measure", "prevalence")
  bad <- !(measure %in% disease_measures)
  for (r in which(bad)) {
    pc$add("diseases", r, "measure", sprintf("Unknown measure '%s' (use prevalence, probability or incidence_rate).", measure[r]))
  }
  prob <- ifelse(measure == "incidence_rate", 1 - exp(-value), value)
  for (r in which(!bad & !is.na(prob) & value >= 0 & (prob <= 0 | prob >= 1))) {
    pc$add("diseases", r, "value", sprintf("Gives a probability of %g; it must be strictly between 0 and 1.", prob[r]))
  }
  list(id = id, value = value, measure = measure, note = chr_col(d, "note"),
       ids = unique(id[!is.na(id) & !grepl("[|;:]", id) & !(tolower(id) %in% c("all", "others"))]),
       dists = row_dists(d, "diseases", pc, keys = ifelse(is.na(id), NA, paste0("prob:", id)),
                         values = value))
}

# Parse each row with the model's own row checks, reporting failures by row.
check_rows <- function(tab, label, ids, pc, parse) {
  for (r in seq_len(nrow(tab))) {
    e <- tryCatch({
      parse(tab[r, , drop = FALSE])
      NULL
    }, deconflate_error = function(e) e)
    if (!is.null(e)) {
      msg <- sub("^[a-z]+, row [0-9?]+: ", "", conditionMessage(e))
      col <- regmatches(msg, regexpr("`(disease|associated_with|among|adjusted_for)`", msg))
      pc$add(label, tab$row[r], if (length(col)) gsub("`", "", col) else NA, msg)
    }
  }
}

adjusted_other_notes <- function(adj, label, ids, pc) {
  if (is.null(ids)) return(invisible(NULL))
  for (r in which(!is.na(adj))) {
    s <- split_ids(adj[r])
    if (length(s) == 1L && tolower(s) == "all") next
    other <- setdiff(s, ids)
    if (length(other)) {
      pc$add(label, r, "adjusted_for",
             sprintf("%s %s not %s of the model: treated as factors outside the model (ignored; %s).",
                     paste(sprintf("'%s'", other), collapse = ", "), if (length(other) > 1L) "are" else "is",
                     if (length(other) > 1L) "diseases" else "a disease",
                     if (length(intersect(s, ids))) "the adjustment for the model's diseases is used" else "the row is used as crude"),
             "note")
    }
  }
}

check_association_table <- function(a, ids, pc) {
  lab <- "associations"
  if (!has_columns(a, "associations", lab, pc)) return(NULL)
  dis <- chr_col(a, "disease")
  aw <- chr_col(a, "associated_with")
  for (r in which(is.na(dis))) pc$add(lab, r, "disease", "Missing disease.")
  measure <- chr_col(a, "measure", "OR")
  retired <- measure %in% names(retired_measures)
  for (r in which(retired)) pc$add(lab, r, "measure", retired_measures[[measure[r]]])
  bad_m <- !retired & !(measure %in% association_measures)
  for (r in which(bad_m)) {
    pc$add(lab, r, "measure", sprintf("Unknown measure '%s' (use %s).", measure[r],
                                      paste(association_measures, collapse = ", ")))
  }
  value <- num_col(a, "value", lab, pc, required = FALSE)
  for (r in which(!retired & is.na(chr_col(a, "value")))) {
    pc$add(lab, r, "value",
           "Missing value. Leave the pair out (no row) if its association is unknown; an odds ratio of 1 states that the diseases are unrelated.")
  }
  for (r in which(measure %in% c("OR", "RR", "OR_ratio") & !is.na(value) & value <= 0)) {
    pc$add(lab, r, "value", "Odds ratios, risk ratios and odds ratio ratios must be positive.")
  }
  for (r in which(measure == "prevalence" & !is.na(value) & (value <= 0 | value >= 1))) {
    pc$add(lab, r, "value", "A prevalence must lie strictly between 0 and 1.")
  }
  for (r in which(measure %in% c("RD", "phi") & !is.na(value) & (value < -1 | value > 1))) {
    pc$add(lab, r, "value", "Risk differences and phi must lie between -1 and 1.")
  }
  among <- chr_col(a, "among")
  adj <- chr_col(a, "adjusted_for")
  ok <- !is.na(dis) & !retired & !bad_m & !is.na(value)
  if (!is.null(ids) && any(ok)) {
    rows <- data.frame(row = seq_len(nrow(a)), disease = dis, associated_with = aw, among = among,
                       adjusted_for = adj, measure = measure, value = value, stringsAsFactors = FALSE)
    check_rows(rows[ok, , drop = FALSE], lab, ids, pc, function(x) association_spec(x, ids))
    adjusted_other_notes(adj, lab, ids, pc)
  }
  list(disease = dis, associated_with = aw, among = among, adjusted_for = adj, measure = measure,
       value = value, study = chr_col(a, "study"), note = chr_col(a, "note"),
       dists = row_dists(a, lab, pc, keys = paste0("associations:", seq_len(nrow(a))), values = value))
}

# Columns of an impact table that are arguments of cm_read_inputs() now.
meta_columns <- c("outcome", "label", "units")

# The outcome, label and units columns of an impact table: not read (a note
# with the argument to use), but different values in one table are an
# error, since a table holds one outcome in one unit.
meta_column_notes <- function(im, label, pc) {
  for (col in intersect(meta_columns, names(im))) {
    v <- unique(stats::na.omit(chr_col(im, col)))
    arg <- if (col == "units") "units" else "outcome"
    if (length(v) > 1L) {
      pc$add(label, NA, col, if (col == "units") {
        sprintf("The rows have different units (found: %s): every row of an impact table must be in the same units.",
                paste(v, collapse = ", "))
      } else {
        sprintf("This table holds several outcomes (found: %s): keep one impact table per outcome, and adjust each in turn.",
                paste(v, collapse = ", "))
      })
    } else {
      pc$add(label, NA, col, sprintf("Column `%s` is no longer read: give it as `cm_read_inputs(%s = %s)`.", col, arg,
                                     if (length(v)) sprintf("\"%s\"", v) else "..."), "note")
    }
  }
}

check_impact_table <- function(im, label, ids, pc) {
  if (!has_columns(im, "impacts", label, pc)) return(NULL)
  dis <- chr_col(im, "disease")
  value <- num_col(im, "value", label, pc)
  for (r in which(is.na(dis))) pc$add(label, r, "disease", "Missing disease.")
  # Event impacts have a measure column (with at least one entry).
  meas <- chr_col(im, "measure")
  event <- any(!is.na(meas))
  among <- chr_col(im, "among")
  adj <- chr_col(im, "adjusted_for")
  if (event) {
    for (r in which(is.na(meas))) {
      pc$add(label, r, "measure", sprintf("Missing measure: event impacts need a measure in every row (%s).",
                                          paste(event_measures, collapse = ", ")))
    }
    for (r in which(!is.na(meas) & !(meas %in% event_measures))) {
      pc$add(label, r, "measure", sprintf("Unknown measure '%s' (use %s).", meas[r],
                                          paste(event_measures, collapse = ", ")))
    }
    for (r in which(meas %in% setdiff(event_measures, "RD") & !is.na(value) & value <= 0)) {
      pc$add(label, r, "value", "Ratios must be positive.")
    }
    for (r in which(meas %in% "RD" & !is.na(value) & (value <= -1 | value >= 1))) {
      pc$add(label, r, "value", "A risk difference must lie between -1 and 1.")
    }
    no_effect <- "1 for ratios, 0 for risk differences"
  } else {
    no_effect <- "0"
  }
  if (!is.null(ids)) {
    ok <- !is.na(dis)
    rows <- data.frame(row = seq_len(nrow(im)), disease = dis, among = among, adjusted_for = adj,
                       measure = meas, value = value, stringsAsFactors = FALSE)
    if (any(ok)) check_rows(rows[ok, , drop = FALSE], label, ids, pc, function(x) impact_spec(x, ids))
    adjusted_other_notes(adj, label, ids, pc)
    own <- unique(trimws(dis[!is.na(dis) & !grepl(";", dis, fixed = TRUE)]))
    miss <- setdiff(ids, own)
    if (length(miss)) {
      pc$add(label, NA, NA, sprintf("No impact for: %s. Every disease needs a row of its own (for no effect: %s).",
                                    paste(miss, collapse = ", "), no_effect))
    }
  }
  meta_column_notes(im, label, pc)
  list(disease = dis, value = value, among = among, adjusted_for = adj,
       measure = if (event) meas else NULL, study = chr_col(im, "study"),
       note = chr_col(im, "note"),
       dists = row_dists(im, label, pc, keys = paste0("impacts:", seq_len(nrow(im))), values = value))
}

run_input_checks <- function(diseases, associations, impacts, outcome = NULL, units = NULL) {
  pc <- new_problems()
  tabs <- list(
    diseases = read_checked(diseases, "diseases", pc),
    associations = read_checked(associations, "associations", pc),
    impacts = read_checked(impacts, "impacts", pc)
  )
  dd <- if (is.null(tabs$diseases)) {
    pc$add("diseases", NA, NA, "A diseases table is required.")
    NULL
  } else check_disease_table(tabs$diseases, pc)
  ids <- dd$ids
  aa <- if (!is.null(tabs$associations)) check_association_table(tabs$associations, ids, pc)
  ii <- if (is.null(tabs$impacts)) {
    pc$add("impacts", NA, NA, "An impact table is required.")
    NULL
  } else check_impact_table(tabs$impacts, "impacts", ids, pc)
  dists <- c(dd$dists %||% list(), aa$dists %||% list(), ii$dists %||% list())

  out <- list(model = NULL, tables = tabs)
  if (is.null(dd) || is.null(ii) || pc$has_errors()) {
    out$problems <- pc$table()
    return(out)
  }
  built <- tryCatch({
    dis_obj <- cm_diseases(dd$id, dd$value, measure = dd$measure, note = dd$note)
    assoc_obj <- if (!is.null(aa)) {
      cm_associations(aa$disease, aa$associated_with, aa$value, measure = aa$measure, among = aa$among,
                      adjusted_for = aa$adjusted_for, study = aa$study, note = aa$note)
    }
    pop <- cm_population(dis_obj, assoc_obj)
    imp <- cm_impacts(ii$disease, ii$value, among = ii$among, adjusted_for = ii$adjusted_for,
                      measure = ii$measure, study = ii$study, note = ii$note,
                      outcome = outcome, units = units)
    m <- cm_model(pop, imp, distributions = if (length(dists)) dists else NULL)
    pair_tables(m)
    m
  }, error = function(e) e)
  if (inherits(built, "condition")) {
    pc$add("model", NA, NA, conditionMessage(built))
  } else {
    out$model <- built
  }
  out$problems <- pc$table()
  out
}
