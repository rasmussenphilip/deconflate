#' Read model inputs from CSV files or data frames
#'
#' Builds a model ([cm_model()]) from one table of each kind: the diseases,
#' their associations, one impact table and, optionally, interactions and
#' three-way terms. Each table can be a path to a CSV file (with any name) or
#' a data frame typed in R. Values with a distribution (see Uncertainty) are
#' kept in the model, and [deconflate()] uses them to compute intervals.
#' [cm_template()] writes an example set to start from.
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
#' can also have a free-text `note` column.
#'
#' **diseases**: `id`, `value` (required); `type` (`prevalence` (default),
#' `probability` or `incidence_rate`, converted with `1 - exp(-value)`),
#' `time_horizon`, `reference_population`, `source`. Ids must not contain
#' `|`, `;` or `:`, and `all` is reserved.
#'
#' **associations**: `disease1`, `disease2`, `value` (required); `measure`
#' (`OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 | disease2)) or
#' `phi`), `adjusted` (TRUE/FALSE: covariate-adjusted measures are rejected
#' unless `adjusted_associations = "use_as_marginal"`), `adjusted_for`,
#' `source`. Pairs without a row are unknown: the global model fills in their
#' association from the others. An odds ratio of 1 states that two diseases
#' are unrelated.
#'
#' **impacts**: `disease`, `value` (required); `estimand`, `adjusted_for`,
#' `measure`, `source`, `label`, `units` (one label and one unit per table).
#' Every disease needs a row. Without a `measure` column the impacts are
#' additive, in any units (the same for every row; results come back in
#' those units): `estimand` is `crude` (default) or `adjusted_linear` (with
#' `adjusted_for`: ids separated by `;`, or `all`), and a disease with no
#' impact has the value 0. With a `measure` column the impacts are event
#' impacts (see [cm_impacts()]): every row needs a measure (`HR`,
#' `rate_ratio`, `RR`, `OR` or `RD`) and an estimand (`snapshot_crude`, or
#' `snapshot_stratified` with `adjusted_for`), and they are adjusted with
#' `deconflate(..., event_model = TRUE)`.
#'
#' **interactions** (additive impacts only): `disease1`, `disease2`, `value`
#' (required); `source`.
#'
#' **three_way**: `disease1`, `disease2`, `disease3`, `ratio` (required);
#' `source`. See [cm_three_way()].
#'
#' @section Uncertainty:
#' Every table can have the columns `dist` and `p1`-`p4`. A row with a
#' `dist` has an uncertain value: [deconflate()] draws it from its
#' distribution (`n_draws`), while the point value (`value`, or `ratio` for
#' three-way terms) gives the central estimate. A row with an empty `dist` is
#' a fixed point value, so point values and uncertain values can be mixed
#' freely. Values are drawn on the scale they were entered on (e.g. an
#' incidence rate, an odds ratio, a hazard ratio or a percent impact). A
#' point value outside the support of its distribution is noted.
#' Distributions and their parameters:
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
#' @param diseases,associations,impacts,interactions,three_way Paths to CSV
#'   files (any names) or data frames. `diseases` and `impacts` are required;
#'   `associations` can be left out only for the sensitivity tools
#'   ([deconflate()] needs at least one association).
#' @param adjusted_associations Passed to [cm_population()].
#' @return A [cm_model()], with the distributions of uncertain values in
#'   `$distributions` and any notes from the checks in attribute `"problems"`.
#' @export
#' @examples
#' dir <- file.path(tempdir(), "deconflate-inputs")
#' cm_template(dir, overwrite = TRUE)
#' m <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
#'                     associations = file.path(dir, "associations.csv"),
#'                     impacts = file.path(dir, "yield.csv"))
#' m
#' deconflate(m, n_draws = 0)
#'
#' # Example files shipped with the package
#' ex <- system.file("extdata", "five_diseases", package = "deconflate")
#' five <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
#'                        associations = file.path(ex, "associations.csv"),
#'                        impacts = file.path(ex, "impacts_yield.csv"),
#'                        three_way = file.path(ex, "three_way.csv"))
#' five
#' err <- system.file("extdata", "example_with_errors", package = "deconflate")
#' cm_check_inputs(diseases = file.path(err, "diseases.csv"),
#'                 associations = file.path(err, "associations.csv"),
#'                 impacts = file.path(err, "impacts_yield.csv"))
#'
#' # The same kind of tables typed in R
#' m2 <- cm_read_inputs(
#'   diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
#'   associations = data.frame(disease1 = "d1", disease2 = "d2", value = 2),
#'   impacts = data.frame(disease = c("d1", "d2"), value = c(2.5, 5), units = "% of yield")
#' )
#' deconflate(m2)
cm_read_inputs <- function(diseases, associations = NULL, impacts, interactions = NULL,
                           three_way = NULL,
                           adjusted_associations = c("error", "use_as_marginal")) {
  if (missing(diseases)) cm_abort("`diseases` is required: a path to a CSV file or a data frame.")
  if (missing(impacts)) cm_abort("`impacts` is required: a path to a CSV file or a data frame.")
  adjusted_associations <- match.arg(adjusted_associations)
  chk <- run_input_checks(diseases, associations, three_way, impacts, interactions,
                          adjusted_associations)
  pr <- chk$problems
  if (any(pr$severity == "error")) {
    cm_abort(format_problems(pr[pr$severity == "error", , drop = FALSE]),
             class = "deconflate_input_problems")
  }
  if (nrow(pr)) message(format_problems(pr))
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
#'   diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 1.5)),
#'   impacts = data.frame(disease = c("d1", "d3"), value = c(2, "x"))
#' )
cm_check_inputs <- function(diseases = NULL, associations = NULL, impacts = NULL,
                            interactions = NULL, three_way = NULL,
                            adjusted_associations = c("error", "use_as_marginal")) {
  adjusted_associations <- match.arg(adjusted_associations)
  run_input_checks(diseases, associations, three_way, impacts, interactions,
                   adjusted_associations)$problems
}

#' Build distributions from a table
#'
#' Converts a table with columns `key`, `dist` and `p1`-`p4` (distributions
#' as in the Uncertainty section of [cm_read_inputs()]) into a named list of
#' `cm_dist` objects, e.g. for the `distributions` argument of [cm_model()].
#' Keys are kept as given.
#'
#' @param x A path to a CSV file or a data frame.
#' @return A named list of `cm_dist` objects (names = keys).
#' @export
#' @examples
#' cm_dist_table(data.frame(key = c("impact:d1", "assoc:d1:d2"),
#'                          dist = c("normal", "lognormal_ci"),
#'                          p1 = c(2.5, 2), p2 = c(0.5, 1.4), p3 = c(NA, 2.9)))
cm_dist_table <- function(x) {
  tab <- read_input_table(x)
  if (is.null(tab) || !nrow(tab)) cm_abort("The uncertainty table is empty.")
  if (!all(c("key", "dist") %in% names(tab))) cm_abort("The uncertainty table needs columns `key` and `dist`.")
  key <- clean_chr(tab$key)
  dist <- tolower(clean_chr(tab$dist))
  P <- dist_params(tab)
  errs <- character(0)
  out <- list()
  for (r in seq_len(nrow(tab))) {
    if (is.na(key[r])) {
      errs <- c(errs, sprintf("row %d: missing key.", r))
      next
    }
    if (!is.null(out[[key[r]]])) {
      errs <- c(errs, sprintf("row %d (%s): the key is given more than once.", r, key[r]))
      next
    }
    d <- tryCatch(input_dist(dist[r], P[r, ]), error = function(e) e)
    if (inherits(d, "condition")) {
      errs <- c(errs, sprintf("row %d (%s): %s", r, key[r], conditionMessage(d)))
    } else {
      out[[key[r]]] <- d
    }
  }
  if (length(errs)) cm_abort(paste(c("Problems in the uncertainty table:", errs), collapse = "\n  "))
  out
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
#'   units can be used);
#' * `yield_interactions.csv`: an empty interactions table for the yield
#'   impacts (pairwise interactions in the same units);
#' * `culling.csv`: event impacts on culling (hazard ratios and a risk ratio),
#'   for `deconflate(..., event_model = TRUE, overall_risk = ...)`;
#' * `three_way.csv`: an empty three-way table.
#'
#' The folder `system.file("extdata", "five_diseases", package =
#' "deconflate")` has a fuller example that uses every feature.
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
#'                         impacts = file.path(dir, "yield.csv"))
#' deconflate(yield, n_draws = 200, seed = 1)
#' culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
#'                           associations = file.path(dir, "associations.csv"),
#'                           impacts = file.path(dir, "culling.csv"))
#' deconflate(culling, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
cm_template <- function(dir, overwrite = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  src <- "Illustrative values"
  tabs <- list(
    diseases = data.frame(
      id = c("LAM", "SCK", "MET"), value = c(0.25, 0.48, 0.10),
      type = c("prevalence", "incidence_rate", "prevalence"),
      time_horizon = "lactation", source = src,
      dist = c("beta", NA, NA), p1 = c(78.29, NA, NA), p2 = c(227.42, NA, NA),
      p3 = NA, p4 = NA, note = c("beta(shape1, shape2)", NA, NA), stringsAsFactors = FALSE),
    associations = data.frame(
      disease1 = c("LAM", "MET", "MET"), disease2 = c("SCK", "SCK", "LAM"),
      value = c(2.01, 1.94, 6.10), measure = "OR", source = src,
      dist = c("normal", NA, NA), p1 = c(2.01, NA, NA), p2 = c(0.20, NA, NA),
      p3 = c(0, NA, NA), p4 = NA, note = c("normal(mean, sd) truncated at 0", NA, NA),
      stringsAsFactors = FALSE),
    yield = data.frame(
      disease = c("LAM", "SCK", "MET"), value = c(4.81, 8.40, 5.61), estimand = "crude",
      adjusted_for = NA, label = "milk yield loss", units = "% of yield", source = src,
      dist = c("normal", "normal", NA), p1 = c(4.81, 8.40, NA), p2 = c(0.87, 1.19, NA),
      p3 = NA, p4 = NA, note = c("normal(mean, sd)", "normal(mean, sd)", NA),
      stringsAsFactors = FALSE),
    yield_interactions = data.frame(
      disease1 = character(0), disease2 = character(0), value = numeric(0),
      source = character(0), dist = character(0), p1 = numeric(0), p2 = numeric(0),
      p3 = numeric(0), p4 = numeric(0), stringsAsFactors = FALSE),
    culling = data.frame(
      disease = c("LAM", "SCK", "MET"), value = c(1.74, 1.92, 1.45),
      measure = c("HR", "HR", "RR"), estimand = "snapshot_crude", adjusted_for = NA,
      label = "culling", units = NA, source = src,
      dist = c("lognormal_ci", NA, NA), p1 = c(1.74, NA, NA), p2 = c(1.45, NA, NA),
      p3 = c(2.09, NA, NA), p4 = NA,
      note = c("hazard ratio; lognormal_ci(estimate, lower, upper)", "hazard ratio",
               "risk ratio of culling within the lactation"),
      stringsAsFactors = FALSE),
    three_way = data.frame(
      disease1 = character(0), disease2 = character(0), disease3 = character(0),
      ratio = numeric(0), source = character(0), dist = character(0), p1 = numeric(0),
      p2 = numeric(0), p3 = numeric(0), p4 = numeric(0), stringsAsFactors = FALSE)
  )
  paths <- file.path(dir, paste0(names(tabs), ".csv"))
  exists <- file.exists(paths)
  if (any(exists) && !overwrite) {
    cm_abort(sprintf("Files already exist (use overwrite = TRUE): %s",
                     paste(basename(paths[exists]), collapse = ", ")))
  }
  for (k in seq_along(tabs)) utils::write.csv(tabs[[k]], paths[k], row.names = FALSE, na = "")
  message(sprintf("Wrote %d files to %s", length(paths), dir))
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
  diseases = list(required = c("id", "value"),
                  optional = c("type", "time_horizon", "reference_population", "source", "note",
                               dist_columns)),
  associations = list(required = c("disease1", "disease2", "value"),
                      optional = c("measure", "adjusted", "adjusted_for", "source", "note",
                                   dist_columns)),
  three_way = list(required = c("disease1", "disease2", "disease3", "ratio"),
                   optional = c("source", "note", dist_columns)),
  impacts = list(required = c("disease", "value"),
                 optional = c("estimand", "adjusted_for", "measure", "source", "label", "units",
                              "note", dist_columns)),
  interactions = list(required = c("disease1", "disease2", "value"),
                      optional = c("source", "note", dist_columns))
)

# Read one table from a path or a data frame, with normalised column names.
read_input_table <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) cm_abort(sprintf("File not found: %s", x))
    # An empty file is an empty table. Columns are read as text and converted
    # by the checks, so that ids such as "T" or "01" are not type-converted.
    if (!isTRUE(file.size(x) > 0)) return(data.frame())
    x <- utils::read.csv(x, stringsAsFactors = FALSE, na.strings = c("", "NA"),
                         strip.white = TRUE, check.names = FALSE,
                         colClasses = "character", fileEncoding = "UTF-8-BOM")
  }
  if (!is.data.frame(x)) cm_abort("Each table must be a data frame or the path to a CSV file.")
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

dist_params <- function(tab) {
  P <- vapply(paste0("p", 1:4), function(col) {
    v <- tab[[col]]
    if (is.null(v)) rep(NA_real_, nrow(tab)) else suppressWarnings(as.numeric(v))
  }, numeric(nrow(tab)))
  matrix(P, nrow = nrow(tab))
}

# Build one distribution from a `dist` name and parameters p1-p4.
input_dist <- function(dist, p) {
  p <- as.numeric(p)
  length(p) <- 4L
  need <- function(k, what) {
    if (any(is.na(p[seq_len(k)]))) cm_abort(sprintf("'%s' needs %s.", dist, what))
  }
  opt <- function(i, default) if (is.na(p[i])) default else p[i]
  if (is.na(dist)) cm_abort("Missing distribution name.")
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
    cm_abort(sprintf("Unknown distribution '%s' (use fixed, normal, lognormal, lognormal_ci, beta, pert, pert_mean or uniform).", dist))
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
  miss <- setdiff(sp$required, names(tab))
  for (col in miss) pc$add(label, NA, col, "Required column is missing.")
  extra <- setdiff(names(tab), c(sp$required, sp$optional))
  if (spec == "impacts" && "outcome" %in% extra) {
    pc$add(label, NA, "outcome",
           "An impact table holds one outcome: split this table into one table per outcome, drop the outcome column, and adjust each table in turn.")
    extra <- setdiff(extra, "outcome")
  }
  if (spec == "impacts" && any(c("scale", "direction") %in% extra)) {
    pc$add(label, NA, NA,
           "Columns scale and direction are no longer used: impacts are adjusted in their own units; convert the adjusted impacts as needed afterwards.",
           "note")
    extra <- setdiff(extra, c("scale", "direction"))
  }
  if (spec == "associations" && any(c("n11", "n10", "n01", "n00") %in% extra)) {
    pc$add(label, NA, NA,
           "Contingency-table counts (n11, n10, n01, n00) are no longer read: compute each table's odds ratio, n11 * n00 / (n10 * n01), and enter it in `value` with measure OR.")
    extra <- setdiff(extra, c("n11", "n10", "n01", "n00"))
  }
  for (col in extra) pc$add(label, NA, col, "Unrecognised column (ignored).", "note")
  !length(miss)
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
        pc$add(label, r, "dist", "Parameters p1-p4 are given without a distribution (dist).")
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
  if (required) for (r in which(is.na(x) & !bad)) pc$add(label, r, col, "Missing value.")
  x
}

chr_col <- function(tab, col, default = NA_character_) clean_chr(tab[[col]], default, nrow(tab))

check_id_col <- function(v, label, col, ids, pc) {
  for (r in which(is.na(v))) pc$add(label, r, col, "Missing disease id.")
  if (!is.null(ids)) {
    for (r in which(!is.na(v) & !(v %in% ids))) {
      pc$add(label, r, col, sprintf("Unknown disease '%s' (not in the diseases table).", v[r]))
    }
  }
}

check_adjusted_for <- function(adj, label, ids, pc) {
  for (r in which(!is.na(adj))) {
    s <- split_ids(adj[r])
    if (length(s) == 1L && tolower(s) == "all") next
    unk <- setdiff(s, ids)
    if (length(unk) && !is.null(ids)) {
      pc$add(label, r, "adjusted_for", sprintf("Unknown disease(s): %s.", paste(unk, collapse = ", ")))
    }
  }
}

check_disease_table <- function(d, pc) {
  if (!has_columns(d, "diseases", "diseases", pc)) return(NULL)
  id <- chr_col(d, "id")
  for (r in which(is.na(id))) pc$add("diseases", r, "id", "Missing disease id.")
  for (r in which(duplicated(id) & !is.na(id))) {
    pc$add("diseases", r, "id", sprintf("Duplicate disease id '%s'.", id[r]))
  }
  for (r in which(!is.na(id) & grepl("[|;:]", id))) {
    pc$add("diseases", r, "id", sprintf("Disease id '%s' must not contain '|', ';' or ':'.", id[r]))
  }
  for (r in which(!is.na(id) & tolower(id) == "all")) {
    pc$add("diseases", r, "id", "'all' is reserved and cannot be a disease id.")
  }
  value <- num_col(d, "value", "diseases", pc)
  for (r in which(!is.na(value) & value < 0)) pc$add("diseases", r, "value", "Must not be negative.")
  type <- chr_col(d, "type", "prevalence")
  bad_type <- !(type %in% c("prevalence", "probability", "incidence_rate"))
  for (r in which(bad_type)) {
    pc$add("diseases", r, "type", sprintf("Unknown type '%s' (use prevalence, probability or incidence_rate).", type[r]))
  }
  prob <- ifelse(type == "incidence_rate", 1 - exp(-value), value)
  for (r in which(!bad_type & !is.na(prob) & value >= 0 & (prob <= 0 | prob >= 1))) {
    pc$add("diseases", r, "value", sprintf("Gives a probability of %g; it must be strictly between 0 and 1.", prob[r]))
  }
  hz <- unique(stats::na.omit(chr_col(d, "time_horizon")))
  if (length(hz) > 1L) {
    pc$add("diseases", NA, "time_horizon",
           sprintf("Diseases refer to different time horizons (%s); probabilities, associations and impacts should share one period.",
                   paste(hz, collapse = ", ")), "note")
  }
  list(id = id, value = value, type = type, time_horizon = chr_col(d, "time_horizon"),
       reference_population = chr_col(d, "reference_population"), source = chr_col(d, "source"),
       ids = unique(id[!is.na(id)]),
       dists = row_dists(d, "diseases", pc, keys = ifelse(is.na(id), NA, paste0("prob:", id)),
                         values = value))
}

check_association_table <- function(a, ids, adjusted_associations, pc) {
  lab <- "associations"
  if (!has_columns(a, "associations", lab, pc)) return(NULL)
  a1 <- chr_col(a, "disease1")
  a2 <- chr_col(a, "disease2")
  check_id_col(a1, lab, "disease1", ids, pc)
  check_id_col(a2, lab, "disease2", ids, pc)
  for (r in which(!is.na(a1) & !is.na(a2) & a1 == a2)) {
    pc$add(lab, r, "disease2", "An association needs two different diseases.")
  }
  key <- ifelse(is.na(a1) | is.na(a2), NA, pair_key(a1, a2))
  for (r in which(duplicated(key) & !is.na(key))) {
    pc$add(lab, r, NA, sprintf("Duplicate pair %s:%s.", a1[r], a2[r]))
  }
  measure <- chr_col(a, "measure", "OR")
  retired <- measure %in% names(retired_measures)
  for (r in which(retired)) pc$add(lab, r, "measure", retired_measures[[measure[r]]])
  for (r in which(!retired & !(measure %in% association_measures))) {
    pc$add(lab, r, "measure", sprintf("Unknown measure '%s' (use %s).", measure[r],
                                      paste(association_measures, collapse = ", ")))
  }
  value <- num_col(a, "value", lab, pc, required = FALSE)
  # (A value that is not a number is reported by num_col().)
  for (r in which(!retired & is.na(chr_col(a, "value")))) {
    pc$add(lab, r, "value",
           "Missing value. Leave the pair out (no row) if its association is unknown; an odds ratio of 1 states that the diseases are unrelated.")
  }
  for (r in which(measure %in% c("OR", "RR") & !is.na(value) & value <= 0)) {
    pc$add(lab, r, "value", "Odds ratios and risk ratios must be positive.")
  }
  for (r in which(measure == "cond_prob" & !is.na(value) & (value < 0 | value > 1))) {
    pc$add(lab, r, "value", "A conditional probability must be between 0 and 1.")
  }
  for (r in which(measure %in% c("RD", "phi") & !is.na(value) & (value < -1 | value > 1))) {
    pc$add(lab, r, "value", "Must be between -1 and 1.")
  }
  adj <- chr_col(a, "adjusted", "FALSE")
  adj_l <- toupper(adj) %in% c("TRUE", "T", "YES", "1")
  for (r in which(!(toupper(adj) %in% c("TRUE", "T", "YES", "1", "FALSE", "F", "NO", "0")))) {
    pc$add(lab, r, "adjusted", sprintf("'%s' is not TRUE or FALSE.", adj[r]))
  }
  for (r in which(adj_l)) {
    if (adjusted_associations == "error") {
      pc$add(lab, r, "adjusted",
             "A covariate-adjusted measure is not a marginal 2x2 association. Use the crude measure, or set adjusted_associations = 'use_as_marginal' to use it as an approximation.")
    } else {
      pc$add(lab, r, "adjusted", "Covariate-adjusted measure used as if it were marginal (approximation).",
             "note")
    }
  }
  list(a1 = a1, a2 = a2, value = value, measure = measure, adjusted = adj_l,
       adjusted_for = chr_col(a, "adjusted_for"), source = chr_col(a, "source"),
       dists = row_dists(a, lab, pc, keys = ifelse(is.na(key), NA, paste0("assoc:", a1, ":", a2)),
                         values = value))
}

check_three_way_table <- function(t3, ids, pc) {
  lab <- "three_way"
  if (!has_columns(t3, "three_way", lab, pc)) return(NULL)
  d <- lapply(c("disease1", "disease2", "disease3"), function(col) {
    v <- chr_col(t3, col)
    check_id_col(v, lab, col, ids, pc)
    v
  })
  for (r in seq_len(nrow(t3))) {
    x <- c(d[[1]][r], d[[2]][r], d[[3]][r])
    if (!anyNA(x) && anyDuplicated(x)) pc$add(lab, r, NA, "A three-way term needs three different diseases.")
  }
  key <- vapply(seq_len(nrow(t3)), function(r) {
    x <- c(d[[1]][r], d[[2]][r], d[[3]][r])
    if (anyNA(x)) NA_character_ else paste(sort(x), collapse = "|")
  }, character(1))
  for (r in which(duplicated(key) & !is.na(key))) pc$add(lab, r, NA, "Duplicate three-way term.")
  ratio <- num_col(t3, "ratio", lab, pc)
  for (r in which(!is.na(ratio) & ratio <= 0)) pc$add(lab, r, "ratio", "Must be positive.")
  tkeys <- ifelse(is.na(d[[1]]) | is.na(d[[2]]) | is.na(d[[3]]), NA,
                  paste("three", d[[1]], d[[2]], d[[3]], sep = ":"))
  list(d1 = d[[1]], d2 = d[[2]], d3 = d[[3]], ratio = ratio, source = chr_col(t3, "source"),
       dists = row_dists(t3, lab, pc, keys = tkeys, values = ratio))
}

check_impact_table <- function(im, label, ids, pc) {
  if (!has_columns(im, "impacts", label, pc)) return(NULL)
  dis <- chr_col(im, "disease")
  value <- num_col(im, "value", label, pc)
  check_id_col(dis, label, "disease", ids, pc)
  for (r in which(duplicated(dis) & !is.na(dis))) {
    pc$add(label, r, "disease", sprintf("Disease '%s' has more than one impact.", dis[r]))
  }
  # Event impacts have a measure column (with at least one entry).
  meas <- chr_col(im, "measure")
  event <- any(!is.na(meas))
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
    est <- chr_col(im, "estimand")
    for (r in which(!(est %in% c("snapshot_crude", "snapshot_stratified")))) {
      pc$add(label, r, "estimand",
             if (is.na(est[r])) "Missing estimand: event impacts need snapshot_crude or snapshot_stratified (see ?cm_impacts)." else
               sprintf("Unknown estimand '%s' for event impacts (use snapshot_crude or snapshot_stratified; see ?cm_impacts).", est[r]))
    }
    for (r in which(est %in% "snapshot_crude" & !is.na(adj))) {
      pc$add(label, r, "adjusted_for", "adjusted_for is given for a snapshot_crude estimate; set estimand = snapshot_stratified.")
    }
    for (r in which(est %in% "snapshot_stratified" & is.na(adj))) {
      pc$add(label, r, "adjusted_for", "estimand = snapshot_stratified needs adjusted_for (disease ids or all).")
    }
    for (r in which(meas %in% setdiff(event_measures, "RD") & !is.na(value) & value <= 0)) {
      pc$add(label, r, "value", "Ratios must be positive.")
    }
    for (r in which(meas %in% "RD" & !is.na(value) & (value <= -1 | value >= 1))) {
      pc$add(label, r, "value", "A risk difference must lie between -1 and 1.")
    }
    no_effect <- "1 for ratios, 0 for risk differences"
  } else {
    est <- chr_col(im, "estimand", "crude")
    for (r in which(est %in% c("snapshot_crude", "snapshot_stratified"))) {
      pc$add(label, r, "estimand",
             sprintf("'%s' is an estimand of event impacts: add a measure column (e.g. HR) and use event_model = TRUE in deconflate().", est[r]))
    }
    for (r in which(!(est %in% c("crude", "adjusted_linear", "snapshot_crude", "snapshot_stratified")))) {
      pc$add(label, r, "estimand",
             sprintf("Unsupported estimand '%s' (use crude, or adjusted_linear for coefficients of an additive regression).", est[r]))
    }
    for (r in which(est == "crude" & !is.na(adj))) {
      pc$add(label, r, "adjusted_for",
             "adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.")
    }
    for (r in which(est == "adjusted_linear" & is.na(adj))) {
      pc$add(label, r, "adjusted_for", "estimand = adjusted_linear needs adjusted_for (disease ids or all).")
    }
    no_effect <- "0"
  }
  check_adjusted_for(adj, label, ids, pc)
  if (!is.null(ids)) {
    miss <- setdiff(ids, dis)
    if (length(miss)) {
      pc$add(label, NA, NA, sprintf("No impact for: %s. Add a row for each (for no effect: %s).",
                                    paste(miss, collapse = ", "), no_effect))
    }
  }
  one <- function(col) {
    v <- unique(stats::na.omit(chr_col(im, col)))
    if (length(v) > 1L) {
      pc$add(label, NA, col, sprintf("One %s per table is allowed (found: %s).", col,
                                     paste(v, collapse = ", ")))
    }
    if (length(v)) v[1] else NULL
  }
  list(disease = dis, value = value, estimand = est, adjusted_for = adj,
       measure = if (event) meas else NULL,
       source = chr_col(im, "source"), label = one("label"), units = one("units"),
       dists = row_dists(im, label, pc, keys = ifelse(is.na(dis), NA, paste0("impact:", dis)),
                         values = value))
}

check_interaction_table <- function(it, label, ids, pc) {
  if (!has_columns(it, "interactions", label, pc)) return(NULL)
  x1 <- chr_col(it, "disease1")
  x2 <- chr_col(it, "disease2")
  check_id_col(x1, label, "disease1", ids, pc)
  check_id_col(x2, label, "disease2", ids, pc)
  for (r in which(!is.na(x1) & !is.na(x2) & x1 == x2)) {
    pc$add(label, r, "disease2", "An interaction needs two different diseases.")
  }
  key <- ifelse(is.na(x1) | is.na(x2), NA, pair_key(x1, x2))
  for (r in which(duplicated(key) & !is.na(key))) {
    pc$add(label, r, NA, sprintf("Duplicate interaction %s:%s.", x1[r], x2[r]))
  }
  value <- num_col(it, "value", label, pc)
  list(d1 = x1, d2 = x2, value = value,
       source = chr_col(it, "source"),
       dists = row_dists(it, label, pc, keys = ifelse(is.na(key), NA, paste("inter", x1, x2, sep = ":")),
                         values = value))
}

run_input_checks <- function(diseases, associations, three_way, impacts, interactions,
                             adjusted_associations) {
  pc <- new_problems()
  tabs <- list(
    diseases = read_checked(diseases, "diseases", pc),
    associations = read_checked(associations, "associations", pc),
    impacts = read_checked(impacts, "impacts", pc),
    interactions = read_checked(interactions, "interactions", pc),
    three_way = read_checked(three_way, "three_way", pc)
  )
  dd <- if (is.null(tabs$diseases)) {
    pc$add("diseases", NA, NA, "A diseases table is required.")
    NULL
  } else check_disease_table(tabs$diseases, pc)
  ids <- dd$ids
  aa <- if (!is.null(tabs$associations)) {
    check_association_table(tabs$associations, ids, adjusted_associations, pc)
  }
  tt <- if (!is.null(tabs$three_way)) check_three_way_table(tabs$three_way, ids, pc)
  ii <- if (is.null(tabs$impacts)) {
    pc$add("impacts", NA, NA, "An impact table is required.")
    NULL
  } else check_impact_table(tabs$impacts, "impacts", ids, pc)
  xx <- if (!is.null(tabs$interactions)) {
    check_interaction_table(tabs$interactions, "interactions", ids, pc)
  }
  if (!is.null(xx) && !is.null(ii$measure)) {
    pc$add("interactions", NA, NA, "Interactions apply to additive impacts only, not to event impacts (measure column).")
  }
  dists <- c(dd$dists %||% list(), aa$dists %||% list(), tt$dists %||% list(),
             ii$dists %||% list(), xx$dists %||% list())

  out <- list(model = NULL, tables = tabs)
  if (is.null(dd) || is.null(ii) || pc$has_errors()) {
    out$problems <- pc$table()
    return(out)
  }
  built <- tryCatch({
    dis_obj <- suppressWarnings(cm_diseases(dd$id, dd$value, type = dd$type,
                                            time_horizon = dd$time_horizon,
                                            reference_population = dd$reference_population,
                                            source = dd$source))
    assoc_obj <- if (!is.null(aa)) {
      cm_associations(aa$a1, aa$a2, aa$value, measure = aa$measure,
                      adjusted = aa$adjusted, adjusted_for = aa$adjusted_for, source = aa$source)
    }
    tw_obj <- if (!is.null(tt)) cm_three_way(tt$d1, tt$d2, tt$d3, tt$ratio, source = tt$source)
    pop <- cm_population(dis_obj, assoc_obj, tw_obj, adjusted_associations = adjusted_associations)
    pair_tables(pop)
    imp <- cm_impacts(ii$disease, ii$value, estimand = ii$estimand, adjusted_for = ii$adjusted_for,
                      source = ii$source, label = ii$label, units = ii$units, measure = ii$measure)
    int <- if (!is.null(xx)) cm_interactions(xx$d1, xx$d2, xx$value, source = xx$source)
    cm_model(pop, imp, int, distributions = if (length(dists)) dists else NULL)
  }, error = function(e) e)
  if (inherits(built, "condition")) {
    pc$add("model", NA, NA, conditionMessage(built))
  } else {
    out$model <- built
  }
  out$problems <- pc$table()
  out
}
