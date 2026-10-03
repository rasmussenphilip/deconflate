#' Read model inputs from CSV files or data frames
#'
#' Builds a population ([cm_population()]), one analysis per impact table
#' ([cm_analyses()]), optionally a hazard-ratio model ([cm_hr_model()]) and,
#' if uncertainty tables are given, a Monte Carlo sampler. Each table can be
#' a path to a CSV file or a data frame typed in R; alternatively, `dir`
#' names a folder of CSV files (see Files). [cm_template()] writes an example
#' set to start from.
#'
#' All tables are checked before anything is built, and every problem is
#' reported at once, with its table, row and column (row 1 is the first row
#' below the header). [cm_check_inputs()] runs the same checks without
#' stopping.
#'
#' @section Files:
#' In a folder (`dir`), files are found by name:
#' * `diseases.csv` (required), `associations.csv`, `three_way.csv`;
#' * one impact table per analysis: `impacts_<analysis>.csv`, e.g.
#'   `impacts_yield.csv` and `impacts_fertility.csv`. A single `impacts.csv`
#'   is an analysis named `impacts`. All analyses share the diseases and
#'   associations;
#' * optional interactions per analysis: `interactions_<analysis>.csv` (or
#'   `interactions.csv` when there is one analysis);
#' * optional uncertainty: `uncertainty.csv` (shared: disease probabilities,
#'   associations, and impacts keyed by analysis) and/or
#'   `uncertainty_<analysis>.csv` (that analysis's impacts and interactions);
#' * optional `hazard_ratios.csv` (culling or mortality hazard ratios, for
#'   [deconflate_hr()]).
#'
#' In R, `impacts`, `interactions` and `uncertainty` can be one table or a
#' named list of tables, named after the analyses (`uncertainty` may also
#' have an element named `shared`).
#'
#' @section Columns:
#' Column names are not case-sensitive. Optional columns can be left out or
#' left empty; unrecognised columns are ignored (with a note).
#'
#' **diseases**: `id`, `value` (required); `type` (`prevalence` (default),
#' `probability` or `incidence_rate`, converted with `1 - exp(-value)`),
#' `time_horizon`, `reference_population`, `source`. Ids must not contain
#' `|`, `;` or `:`, and `all` is reserved.
#'
#' **associations**: `disease1`, `disease2` (required); `value`, `measure`
#' (`OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 | disease2)), `phi`,
#' `table`, `independent` or `unknown`), `n11`, `n10`, `n01`, `n00` (counts,
#' for `measure = table`; a zero cell gets 0.5 added to every cell, which is
#' reported), `adjusted` (TRUE/FALSE: covariate-adjusted measures are
#' rejected unless `adjusted_associations = "use_as_marginal"`),
#' `adjusted_for`, `source`. Pairs that are not listed are independent (or
#' unknown, see `missing_associations`).
#'
#' **three_way**: `disease1`, `disease2`, `disease3`, `ratio` (required);
#' `source`. See [cm_three_way()].
#'
#' **impacts** (one table per analysis): `disease`, `value` (required);
#' `estimand` (`crude` (default) or `adjusted_linear`), `adjusted_for` (ids
#' separated by `;`, or `all`; only with `adjusted_linear`), `source`,
#' `label`, `units` (one label and one unit per table). Every disease needs a
#' row (use 0 for no impact). Values are in any units, the same for every
#' row; results come back in those units.
#'
#' **interactions**: `disease1`, `disease2`, `value` (required); `source`.
#'
#' **hazard_ratios**: `disease`, `value` (required); `estimand` (`crude`
#' (default) or `adjusted`), `adjusted_for`, `source`. Every disease needs a
#' row (use 1 for no effect).
#'
#' **uncertainty**: `key`, `dist` (required); `p1`-`p4` (parameters),
#' `note`. Keys name the input: `prob:<disease>`, `assoc:<d1>:<d2>`,
#' `impact:<analysis>:<disease>` and `inter:<analysis>:<d1>:<d2>`. In an
#' analysis's own file (`uncertainty_<analysis>.csv`), or when there is only
#' one analysis, the analysis can be left out: `impact:<disease>`,
#' `inter:<d1>:<d2>`. Values are drawn on the scale the input was entered on
#' (e.g. an incidence rate, or a percent impact). Distributions and their
#' parameters:
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
#' @param diseases,associations,three_way,hazard_ratios Paths to CSV files or
#'   data frames.
#' @param impacts,interactions,uncertainty A path or data frame, or a named
#'   list of them (one per analysis; see Files).
#' @param dir Optional folder with CSV files named as in Files. Arguments
#'   given explicitly take precedence over files.
#' @param missing_associations,adjusted_associations Passed to
#'   [cm_population()].
#' @return A `cm_inputs` list with `population`, `analyses` (a
#'   [cm_analyses()] object, or `NULL` without impact tables), `model` (the
#'   [cm_model()] when there is exactly one analysis), `hr_model` (or
#'   `NULL`), `sampler` (a [cm_sampler()] for one analysis, a
#'   [cm_batch_sampler()] for several, or `NULL`), `tables` (as read; columns
#'   of CSV files are read as text and converted by the checks) and
#'   `problems` (notes only, since errors stop).
#' @export
#' @examples
#' dir <- file.path(tempdir(), "deconflate-inputs")
#' cm_template(dir, overwrite = TRUE)
#' inp <- cm_read_inputs(dir = dir)
#' inp
#' deconflate(inp$analyses)
#'
#' # Example files shipped with the package: the 2024 global dairy inputs,
#' # and a set with deliberate errors
#' gd <- cm_read_inputs(dir = system.file("extdata", "global_dairy_2024", package = "deconflate"))
#' gd$analyses
#' cm_check_inputs(dir = system.file("extdata", "example_with_errors", package = "deconflate"))
#'
#' # The same kind of tables typed in R
#' inp2 <- cm_read_inputs(
#'   diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
#'   associations = data.frame(disease1 = "d1", disease2 = "d2", value = 2),
#'   impacts = list(yield = data.frame(disease = c("d1", "d2"), value = c(2.5, 5),
#'                                     units = "% of yield"))
#' )
#' deconflate(inp2$model)
cm_read_inputs <- function(diseases = NULL, associations = NULL, three_way = NULL,
                           impacts = NULL, interactions = NULL, uncertainty = NULL,
                           hazard_ratios = NULL, dir = NULL,
                           missing_associations = c("independent", "unknown"),
                           adjusted_associations = c("error", "use_as_marginal")) {
  missing_associations <- match.arg(missing_associations)
  adjusted_associations <- match.arg(adjusted_associations)
  chk <- run_input_checks(diseases, associations, three_way, impacts, interactions,
                          uncertainty, hazard_ratios, dir, missing_associations,
                          adjusted_associations)
  pr <- chk$problems
  if (any(pr$severity == "error")) {
    cm_abort(format_problems(pr[pr$severity == "error", , drop = FALSE]),
             class = "deconflate_input_problems")
  }
  if (nrow(pr)) message(format_problems(pr))
  chk$problems <- pr
  structure(chk, class = "cm_inputs")
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
cm_check_inputs <- function(diseases = NULL, associations = NULL, three_way = NULL,
                            impacts = NULL, interactions = NULL, uncertainty = NULL,
                            hazard_ratios = NULL, dir = NULL,
                            missing_associations = c("independent", "unknown"),
                            adjusted_associations = c("error", "use_as_marginal")) {
  missing_associations <- match.arg(missing_associations)
  adjusted_associations <- match.arg(adjusted_associations)
  run_input_checks(diseases, associations, three_way, impacts, interactions, uncertainty,
                   hazard_ratios, dir, missing_associations, adjusted_associations)$problems
}

#' Build distributions from an uncertainty table
#'
#' Converts a table with columns `key`, `dist` and `p1`-`p4` (see the
#' Columns section of [cm_read_inputs()]) into a named list of `cm_dist`
#' objects. Keys are kept as given.
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
#' with [cm_read_inputs()]: three dairy diseases with their associations,
#' two analyses (`impacts_yield.csv` in percent of yield and
#' `impacts_calving_interval.csv` in days), culling hazard ratios, an empty
#' interactions file for the yield analysis, an empty three-way file and an
#' uncertainty file. The values are illustrative.
#'
#' @param dir Folder to write to (created if needed).
#' @param overwrite Overwrite existing files?
#' @return The file paths, invisibly.
#' @export
#' @examples
#' cm_template(file.path(tempdir(), "my-inputs"), overwrite = TRUE)
cm_template <- function(dir, overwrite = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  src <- "Illustrative values"
  tabs <- list(
    diseases = data.frame(
      id = c("LAM", "SCK", "MET"), value = c(0.25, 0.48, 0.10),
      type = c("prevalence", "incidence_rate", "prevalence"),
      time_horizon = "lactation", source = src, stringsAsFactors = FALSE),
    associations = data.frame(
      disease1 = c("LAM", "MET", "MET"), disease2 = c("SCK", "SCK", "LAM"),
      value = c(2.01, 1.94, 6.10), measure = "OR", source = src, stringsAsFactors = FALSE),
    three_way = data.frame(
      disease1 = character(0), disease2 = character(0), disease3 = character(0),
      ratio = numeric(0), source = character(0), stringsAsFactors = FALSE),
    impacts_yield = data.frame(
      disease = c("LAM", "SCK", "MET"), value = c(4.81, 8.40, 5.61), estimand = "crude",
      adjusted_for = NA, label = "milk yield loss", units = "% of yield", source = src,
      stringsAsFactors = FALSE),
    impacts_calving_interval = data.frame(
      disease = c("LAM", "SCK", "MET"), value = c(12, 4, 18), estimand = "crude",
      adjusted_for = NA, label = "calving interval increase", units = "days", source = src,
      stringsAsFactors = FALSE),
    interactions_yield = data.frame(
      disease1 = character(0), disease2 = character(0), value = numeric(0),
      source = character(0), stringsAsFactors = FALSE),
    hazard_ratios = data.frame(
      disease = c("LAM", "SCK", "MET"), value = c(1.74, 1.92, 1.50), estimand = "crude",
      adjusted_for = NA, source = src, stringsAsFactors = FALSE),
    uncertainty = data.frame(
      key = c("prob:LAM", "assoc:LAM:SCK", "impact:yield:LAM", "impact:yield:SCK",
              "impact:calving_interval:MET"),
      dist = c("beta", "normal", "normal", "normal", "pert"),
      p1 = c(78.29, 2.01, 4.81, 8.40, 6), p2 = c(227.42, 0.20, 0.87, 1.19, 18),
      p3 = c(NA, 0, NA, NA, 30), p4 = NA,
      note = c("beta(shape1, shape2)", "normal(mean, sd) truncated at 0",
               "normal(mean, sd), percent", "normal(mean, sd), percent",
               "pert(min, mode, max), days"),
      stringsAsFactors = FALSE)
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
print.cm_inputs <- function(x, ...) {
  cat("<cm_inputs>\n")
  print(x$population)
  if (!is.null(x$analyses)) {
    cat(sprintf("  Analyses: %s\n", paste(names(x$analyses$models), collapse = ", ")))
  }
  if (!is.null(x$hr_model)) cat("  Hazard ratios: yes (use $hr_model with deconflate_hr())\n")
  if (inherits(x$sampler, "cm_batch_sampler")) {
    cat(sprintf("  Uncertain inputs: batch sampler over %d analyses (use $sampler with cm_monte_carlo())\n",
                length(x$sampler$samplers)))
  } else if (!is.null(x$sampler)) {
    cat(sprintf("  Uncertain inputs: %d (use $sampler with cm_monte_carlo())\n",
                length(attr(x$sampler, "specs"))))
  }
  invisible(x)
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

input_table_specs <- list(
  diseases = list(required = c("id", "value"),
                  optional = c("type", "time_horizon", "reference_population", "source")),
  associations = list(required = c("disease1", "disease2"),
                      optional = c("value", "measure", "n11", "n10", "n01", "n00",
                                   "adjusted", "adjusted_for", "source")),
  three_way = list(required = c("disease1", "disease2", "disease3", "ratio"),
                   optional = "source"),
  impacts = list(required = c("disease", "value"),
                 optional = c("estimand", "adjusted_for", "source", "label", "units")),
  interactions = list(required = c("disease1", "disease2", "value"), optional = "source"),
  hazard_ratios = list(required = c("disease", "value"),
                       optional = c("estimand", "adjusted_for", "source")),
  uncertainty = list(required = c("key", "dist"),
                     optional = c("p1", "p2", "p3", "p4", "note"))
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

# Collect the sources of a table type that may have one table per analysis:
# returns a named list (name "" for the unsuffixed table) of paths or frames.
collect_sources <- function(arg, dir, stem) {
  if (!is.null(arg)) {
    if (is.data.frame(arg) || (is.character(arg) && length(arg) == 1L)) return(stats::setNames(list(arg), ""))
    if (is.list(arg)) {
      nms <- names(arg)
      if (is.null(nms) || any(!nzchar(nms))) {
        cm_abort(sprintf("`%s` must be a table or a list of tables named after the analyses.", stem))
      }
      return(arg)
    }
    cm_abort(sprintf("`%s` must be a data frame, a path, or a named list of them.", stem))
  }
  if (is.null(dir) || !dir.exists(dir)) return(list())
  files <- list.files(dir, pattern = sprintf("^%s(_.+)?\\.csv$", stem), ignore.case = TRUE)
  if (!length(files)) return(list())
  nm <- sub(sprintf("^%s_?", stem), "", sub("\\.csv$", "", files, ignore.case = TRUE),
            ignore.case = TRUE)
  out <- as.list(file.path(dir, files))
  names(out) <- nm
  out[order(nm)]
}

table_label <- function(stem, name) if (nzchar(name)) paste0(stem, "_", name) else stem

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
           "deconflate 0.2 uses one impact table per analysis: split this table into impacts_<outcome>.csv files (or a named list in R) and drop the outcome column.")
    extra <- setdiff(extra, "outcome")
  }
  if (spec == "impacts" && any(c("scale", "direction") %in% extra)) {
    pc$add(label, NA, NA,
           "Columns scale and direction are no longer used: impacts are adjusted in their own units; give direction and scale to productivity_gap() instead.",
           "note")
    extra <- setdiff(extra, c("scale", "direction"))
  }
  for (col in extra) pc$add(label, NA, col, "Unrecognised column (ignored).", "note")
  !length(miss)
}

num_col <- function(tab, col, label, pc, required = TRUE) {
  v <- tab[[col]]
  n <- nrow(tab)
  if (is.null(v)) return(rep(NA_real_, n))
  bad <- rep(FALSE, n)
  if (is.numeric(v)) {
    x <- as.numeric(v)
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
       ids = unique(id[!is.na(id)]))
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
  valid_m <- c("OR", "RR", "RD", "cond_prob", "phi", "table", "independent", "unknown")
  for (r in which(!(measure %in% valid_m))) {
    pc$add(lab, r, "measure", sprintf("Unknown measure '%s' (use %s).", measure[r],
                                      paste(valid_m, collapse = ", ")))
  }
  value <- num_col(a, "value", lab, pc, required = FALSE)
  needs <- measure %in% c("OR", "RR", "RD", "cond_prob", "phi")
  for (r in which(needs & is.na(value))) {
    pc$add(lab, r, "value", sprintf("Missing value for measure %s.", measure[r]))
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
  cnt <- do.call(cbind, lapply(c("n11", "n10", "n01", "n00"),
                               function(col) num_col(a, col, lab, pc, FALSE)))
  for (r in which(measure == "table")) {
    if (anyNA(cnt[r, ]) || any(cnt[r, ] < 0)) {
      pc$add(lab, r, NA, "measure = table needs non-negative counts n11, n10, n01 and n00.")
    } else if (sum(cnt[r, ]) == 0) {
      pc$add(lab, r, NA, "The contingency table is empty (all counts are zero).")
    } else if (any(cnt[r, ] == 0)) {
      pc$add(lab, r, NA, "The table has a zero cell; 0.5 was added to every cell (Haldane correction).",
             "note")
    }
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
  list(a1 = a1, a2 = a2, value = value, measure = measure, cnt = cnt, adjusted = adj_l,
       adjusted_for = chr_col(a, "adjusted_for"), source = chr_col(a, "source"),
       numeric_keys = stats::setNames(paste(a1, a2, sep = ":"), key)[needs & !is.na(key)])
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
    paste(sort(c(d[[1]][r], d[[2]][r], d[[3]][r])), collapse = "|")
  }, character(1))
  for (r in which(duplicated(key))) pc$add(lab, r, NA, "Duplicate three-way term.")
  ratio <- num_col(t3, "ratio", lab, pc)
  for (r in which(!is.na(ratio) & ratio <= 0)) pc$add(lab, r, "ratio", "Must be positive.")
  list(d1 = d[[1]], d2 = d[[2]], d3 = d[[3]], ratio = ratio, source = chr_col(t3, "source"))
}

check_impact_table <- function(im, label, ids, pc) {
  if (!has_columns(im, "impacts", label, pc)) return(NULL)
  dis <- chr_col(im, "disease")
  value <- num_col(im, "value", label, pc)
  check_id_col(dis, label, "disease", ids, pc)
  for (r in which(duplicated(dis) & !is.na(dis))) {
    pc$add(label, r, "disease", sprintf("Disease '%s' has more than one impact.", dis[r]))
  }
  est <- chr_col(im, "estimand", "crude")
  for (r in which(!(est %in% c("crude", "adjusted_linear")))) {
    pc$add(label, r, "estimand",
           sprintf("Unsupported estimand '%s' (use crude, or adjusted_linear for coefficients of an additive regression).", est[r]))
  }
  adj <- chr_col(im, "adjusted_for")
  for (r in which(est == "crude" & !is.na(adj))) {
    pc$add(label, r, "adjusted_for",
           "adjusted_for is given for a crude estimate. Set estimand = adjusted_linear if this is a coefficient from an additive regression adjusted for those diseases; other adjusted estimands are not supported.")
  }
  for (r in which(est == "adjusted_linear" & is.na(adj))) {
    pc$add(label, r, "adjusted_for", "estimand = adjusted_linear needs adjusted_for (disease ids or all).")
  }
  check_adjusted_for(adj, label, ids, pc)
  if (!is.null(ids)) {
    miss <- setdiff(ids, dis)
    if (length(miss)) {
      pc$add(label, NA, NA, sprintf("No impact for: %s. Add a row with value 0 for no impact.",
                                    paste(miss, collapse = ", ")))
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
       source = chr_col(im, "source"), label = one("label"), units = one("units"))
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
  list(d1 = x1, d2 = x2, value = num_col(it, "value", label, pc),
       source = chr_col(it, "source"),
       keys = stats::setNames(paste(x1, x2, sep = ":"), key)[!is.na(key)])
}

check_hr_table <- function(h, ids, pc) {
  lab <- "hazard_ratios"
  if (!has_columns(h, "hazard_ratios", lab, pc)) return(NULL)
  dis <- chr_col(h, "disease")
  check_id_col(dis, lab, "disease", ids, pc)
  for (r in which(duplicated(dis) & !is.na(dis))) {
    pc$add(lab, r, "disease", sprintf("Disease '%s' has more than one hazard ratio.", dis[r]))
  }
  value <- num_col(h, "value", lab, pc)
  for (r in which(!is.na(value) & value <= 0)) pc$add(lab, r, "value", "Hazard ratios must be positive.")
  est <- chr_col(h, "estimand", "crude")
  for (r in which(!(est %in% c("crude", "adjusted")))) {
    pc$add(lab, r, "estimand", sprintf("Unknown estimand '%s' (use crude or adjusted).", est[r]))
  }
  adj <- chr_col(h, "adjusted_for")
  for (r in which(est == "crude" & !is.na(adj))) {
    pc$add(lab, r, "adjusted_for", "adjusted_for is given for a crude hazard ratio; set estimand = adjusted.")
  }
  for (r in which(est == "adjusted" & is.na(adj))) {
    pc$add(lab, r, "adjusted_for", "estimand = adjusted needs adjusted_for (disease ids or all).")
  }
  check_adjusted_for(adj, lab, ids, pc)
  if (!is.null(ids)) {
    miss <- setdiff(ids, dis)
    if (length(miss)) {
      pc$add(lab, NA, NA, sprintf("No hazard ratio for: %s. Add a row with value 1 for no effect.",
                                  paste(miss, collapse = ", ")))
    }
  }
  list(disease = dis, value = value, estimand = est, adjusted_for = adj,
       source = chr_col(h, "source"))
}

# Check one uncertainty table. `own` is the analysis the table belongs to
# (or "" for the shared table). Returns the distributions with normalised
# keys: population keys as in cm_sampler(), and analysis keys as
# "<analysis>//impact:<disease>" or "<analysis>//inter:<d1>:<d2>".
check_uncertainty_table <- function(un, label, own, ids, assoc_keys, analyses, inter_keys, pc) {
  if (!has_columns(un, "uncertainty", label, pc)) return(list())
  key <- chr_col(un, "key")
  dist <- tolower(chr_col(un, "dist"))
  P <- vapply(paste0("p", 1:4), function(col) num_col(un, col, label, pc, FALSE), numeric(nrow(un)))
  P <- matrix(P, nrow = nrow(un))
  single <- if (nzchar(own)) own else if (length(analyses) == 1L) analyses else NA_character_
  out <- list()
  for (r in seq_len(nrow(un))) {
    k <- key[r]
    if (is.na(k)) {
      pc$add(label, r, "key", "Missing key.")
      next
    }
    parts <- strsplit(k, ":", fixed = TRUE)[[1]]
    norm <- NA_character_
    msg <- NULL
    type <- parts[1]
    if (type == "prob") {
      if (nzchar(own)) {
        msg <- "Disease probabilities are shared by all analyses; put prob: keys in uncertainty.csv."
      } else if (length(parts) == 2L && !is.null(ids) && parts[2] %in% ids) {
        norm <- k
      } else {
        msg <- sprintf("Key '%s' does not match a disease (expected prob:<disease>).", k)
      }
    } else if (type == "assoc") {
      if (nzchar(own)) {
        msg <- "Associations are shared by all analyses; put assoc: keys in uncertainty.csv."
      } else if (length(parts) == 3L && pair_key(parts[2], parts[3]) %in% names(assoc_keys)) {
        # Keys use the orientation of the association row, so that "a:b"
        # and "b:a" are recognised as the same input.
        norm <- paste0("assoc:", assoc_keys[[pair_key(parts[2], parts[3])]])
      } else {
        msg <- sprintf("Key '%s' does not match an association with a numeric measure (expected assoc:<d1>:<d2>).", k)
      }
    } else if (type == "impact") {
      an <- if (length(parts) == 3L) parts[2] else if (length(parts) == 2L) single else NA_character_
      dis <- parts[length(parts)]
      if (length(parts) > 3L) {
        msg <- sprintf("Key '%s' has too many parts (use impact:<analysis>:<disease> or impact:<disease>).", k)
      } else if (length(parts) == 3L && nzchar(own) && parts[2] != own) {
        msg <- sprintf("Key '%s' names analysis '%s', but this table belongs to '%s'.", k, parts[2], own)
      } else if (is.na(an)) {
        msg <- sprintf("Key '%s': name the analysis (impact:<analysis>:<disease>), since there are several.", k)
      } else if (!(an %in% analyses)) {
        msg <- sprintf("Key '%s': unknown analysis '%s'.", k, an)
      } else if (is.null(ids) || !(dis %in% ids)) {
        msg <- sprintf("Key '%s': unknown disease '%s'.", k, dis)
      } else {
        norm <- paste0(an, "//impact:", dis)
      }
    } else if (type == "inter") {
      an <- if (length(parts) == 4L) parts[2] else if (length(parts) == 3L) single else NA_character_
      pr <- parts[(length(parts) - 1L):length(parts)]
      if (length(parts) == 4L && nzchar(own) && parts[2] != own) {
        msg <- sprintf("Key '%s' names analysis '%s', but this table belongs to '%s'.", k, parts[2], own)
      } else if (is.na(an) || length(parts) < 3L) {
        msg <- sprintf("Key '%s': use inter:<analysis>:<d1>:<d2>.", k)
      } else if (!(an %in% analyses)) {
        msg <- sprintf("Key '%s': unknown analysis '%s'.", k, an)
      } else if (!(pair_key(pr[1], pr[2]) %in% names(inter_keys[[an]]))) {
        msg <- sprintf("Key '%s' does not match an interaction of analysis '%s'.", k, an)
      } else {
        norm <- paste0(an, "//inter:", inter_keys[[an]][[pair_key(pr[1], pr[2])]])
      }
    } else {
      msg <- sprintf("Key '%s' must start with prob:, assoc:, impact: or inter:.", k)
    }
    if (!is.null(msg)) pc$add(label, r, "key", msg)
    if (!is.na(norm) && !is.null(out[[norm]])) {
      pc$add(label, r, "key", sprintf("Key '%s': the same input is given more than once.", k))
      norm <- NA_character_
    }
    dd <- tryCatch(input_dist(dist[r], P[r, ]), error = function(e) e)
    if (inherits(dd, "condition")) {
      pc$add(label, r, "dist", conditionMessage(dd))
    } else if (!is.na(norm)) {
      attr(dd, "where") <- c(label = label, row = r)
      out[[norm]] <- dd
    }
  }
  out
}

run_input_checks <- function(diseases, associations, three_way, impacts, interactions,
                             uncertainty, hazard_ratios, dir, missing_associations,
                             adjusted_associations) {
  pc <- new_problems()
  if (!is.null(dir) && !dir.exists(dir)) pc$add("inputs", NA, NA, sprintf("Folder not found: %s", dir))
  single_src <- function(arg, stem) {
    if (!is.null(arg)) return(arg)
    if (!is.null(dir)) {
      f <- file.path(dir, paste0(stem, ".csv"))
      if (file.exists(f)) return(f)
    }
    NULL
  }
  multi <- function(arg, stem) {
    tryCatch(collect_sources(arg, dir, stem), error = function(e) {
      pc$add(stem, NA, NA, conditionMessage(e))
      list()
    })
  }
  tabs <- list(
    diseases = read_checked(single_src(diseases, "diseases"), "diseases", pc),
    associations = read_checked(single_src(associations, "associations"), "associations", pc),
    three_way = read_checked(single_src(three_way, "three_way"), "three_way", pc),
    hazard_ratios = read_checked(single_src(hazard_ratios, "hazard_ratios"), "hazard_ratios", pc)
  )
  imp_src <- multi(impacts, "impacts")
  int_src <- multi(interactions, "interactions")
  unc_src <- multi(uncertainty, "uncertainty")

  # Analysis names: the suffix of impacts_<name>, or "impacts".
  an_names <- ifelse(nzchar(names(imp_src)), names(imp_src), "impacts")
  for (nm in an_names[grepl("[:|; ]", an_names)]) {
    pc$add("impacts", NA, NA, sprintf("Analysis name '%s' must not contain ':', '|', ';' or spaces.", nm))
  }
  for (nm in unique(an_names[duplicated(an_names)])) {
    pc$add("impacts", NA, NA, sprintf("Analysis '%s' is given twice.", nm))
  }
  for (nm in intersect(an_names, c("population", "interactions"))) {
    pc$add("impacts", NA, NA, sprintf("'%s' cannot be used as an analysis name.", nm))
  }
  imp_labels <- stats::setNames(vapply(names(imp_src), function(nm) table_label("impacts", nm),
                                       character(1)), an_names)
  imp_tabs <- lapply(seq_along(imp_src), function(j) read_checked(imp_src[[j]], imp_labels[[j]], pc))
  names(imp_tabs) <- an_names
  imp_tabs <- imp_tabs[!vapply(imp_tabs, is.null, logical(1))]
  analyses <- names(imp_tabs)

  # Interactions: interactions_<analysis>, or an unsuffixed table for a
  # single analysis (or for the analysis named "impacts").
  int_tabs <- list()
  for (j in seq_along(int_src)) {
    nm <- names(int_src)[j]
    lab <- table_label("interactions", nm)
    tab <- read_checked(int_src[[j]], lab, pc)
    if (is.null(tab)) next
    target <- if (nzchar(nm)) nm else if (length(analyses) == 1L) analyses else if ("impacts" %in% analyses) "impacts" else NA
    if (is.na(target)) {
      pc$add(lab, NA, NA, "With several analyses, name interaction tables after their analysis (interactions_<analysis>.csv).")
    } else if (!(target %in% analyses)) {
      pc$add(lab, NA, NA, sprintf("There is no impact table for analysis '%s'.", target))
    } else if (!is.null(int_tabs[[target]])) {
      pc$add(lab, NA, NA, sprintf("Analysis '%s' has more than one interaction table.", target))
    } else {
      attr(tab, "label") <- lab
      int_tabs[[target]] <- tab
    }
  }

  # Diseases and population tables.
  dd <- if (is.null(tabs$diseases)) {
    pc$add("diseases", NA, NA, "A diseases table is required.")
    NULL
  } else check_disease_table(tabs$diseases, pc)
  ids <- dd$ids
  aa <- if (!is.null(tabs$associations)) {
    check_association_table(tabs$associations, ids, adjusted_associations, pc)
  }
  tt <- if (!is.null(tabs$three_way)) check_three_way_table(tabs$three_way, ids, pc)
  hh <- if (!is.null(tabs$hazard_ratios)) check_hr_table(tabs$hazard_ratios, ids, pc)
  ii <- lapply(analyses, function(nm) {
    check_impact_table(imp_tabs[[nm]], imp_labels[[nm]], ids, pc)
  })
  names(ii) <- analyses
  xx <- lapply(names(int_tabs), function(nm) {
    check_interaction_table(int_tabs[[nm]], attr(int_tabs[[nm]], "label"), ids, pc)
  })
  names(xx) <- names(int_tabs)
  inter_keys <- lapply(analyses, function(nm) xx[[nm]]$keys %||% stats::setNames(character(0), character(0)))
  names(inter_keys) <- analyses

  # Uncertainty tables.
  dists <- list()
  for (j in seq_along(unc_src)) {
    nm <- names(unc_src)[j]
    own <- if (nm %in% c("", "shared")) "" else nm
    lab <- table_label("uncertainty", if (nm == "shared") "" else nm)
    if (nzchar(own) && !(own %in% analyses)) {
      pc$add(lab, NA, NA, sprintf("There is no impact table for analysis '%s'.", own))
      next
    }
    tab <- read_checked(unc_src[[j]], lab, pc)
    if (is.null(tab)) next
    d <- check_uncertainty_table(tab, lab, own, ids,
                                 aa$numeric_keys %||% stats::setNames(character(0), character(0)),
                                 analyses, inter_keys, pc)
    for (k in names(d)) {
      if (!is.null(dists[[k]])) {
        w <- attr(d[[k]], "where")
        pc$add(w[["label"]], as.integer(w[["row"]]), "key", "The same input is given more than once.")
      } else {
        dists[[k]] <- d[[k]]
      }
    }
  }
  if (length(dists) && !length(analyses)) {
    pc$add("uncertainty", NA, NA, "Uncertainty needs at least one impact table (a Monte Carlo run adjusts impacts).", "note")
  }

  # Build ---------------------------------------------------------------------
  out <- list(population = NULL, analyses = NULL, model = NULL, hr_model = NULL,
              sampler = NULL,
              tables = c(tabs, list(impacts = imp_tabs, interactions = int_tabs)))
  if (is.null(dd) || pc$has_errors()) {
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
                      n11 = aa$cnt[, 1], n10 = aa$cnt[, 2], n01 = aa$cnt[, 3], n00 = aa$cnt[, 4],
                      adjusted = aa$adjusted, adjusted_for = aa$adjusted_for, source = aa$source)
    }
    tw_obj <- if (!is.null(tt)) cm_three_way(tt$d1, tt$d2, tt$d3, tt$ratio, source = tt$source)
    pop <- cm_population(dis_obj, assoc_obj, tw_obj, missing_associations = missing_associations,
                         adjusted_associations = adjusted_associations)
    pair_tables(pop)
    imps <- lapply(analyses, function(nm) {
      x <- ii[[nm]]
      cm_impacts(x$disease, x$value, estimand = x$estimand, adjusted_for = x$adjusted_for,
                 source = x$source, label = x$label %||% nm, units = x$units)
    })
    names(imps) <- analyses
    ints <- lapply(names(xx), function(nm) {
      x <- xx[[nm]]
      cm_interactions(x$d1, x$d2, x$value, source = x$source)
    })
    names(ints) <- names(xx)
    an <- if (length(imps)) do.call(cm_analyses, c(list(population = pop), imps, list(interactions = ints)))
    hr <- if (!is.null(hh)) {
      cm_hr_model(pop, cm_hazard_ratios(hh$disease, hh$value, estimand = hh$estimand,
                                        adjusted_for = hh$adjusted_for, source = hh$source))
    }
    list(pop = pop, an = an, hr = hr)
  }, error = function(e) e)
  if (inherits(built, "condition")) {
    pc$add("model", NA, NA, conditionMessage(built))
    out$problems <- pc$table()
    return(out)
  }
  out$population <- built$pop
  out$analyses <- built$an
  out$hr_model <- built$hr
  if (!is.null(built$an) && length(built$an$models) == 1L) out$model <- built$an$models[[1]]

  if (length(dists) && !is.null(built$an)) {
    pick <- function(prefix, an = NULL) {
      sel <- if (is.null(an)) startsWith(names(dists), prefix) else startsWith(names(dists), paste0(an, "//", prefix))
      x <- dists[sel]
      nms <- sub("^.*//", "", names(x))
      stats::setNames(x, substring(nms, nchar(prefix) + 1L))
    }
    smp <- tryCatch({
      if (length(analyses) == 1L) {
        cm_sampler(out$model, diseases = pick("prob:"), associations = pick("assoc:"),
                   impacts = pick("impact:", analyses), interactions = pick("inter:", analyses))
      } else {
        imp_d <- lapply(analyses, function(nm) pick("impact:", nm))
        int_d <- lapply(analyses, function(nm) pick("inter:", nm))
        names(imp_d) <- names(int_d) <- analyses
        cm_batch_sampler(built$an, diseases = pick("prob:"), associations = pick("assoc:"),
                         impacts = imp_d, interactions = int_d)
      }
    }, error = function(e) e)
    if (inherits(smp, "condition")) {
      pc$add("uncertainty", NA, NA, conditionMessage(smp))
    } else {
      out$sampler <- smp
    }
  }
  out$problems <- pc$table()
  out
}
