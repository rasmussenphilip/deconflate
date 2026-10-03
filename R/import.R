#' Read model inputs from CSV files or data frames
#'
#' Builds a [cm_model()] (and, if an uncertainty table is given, a
#' [cm_sampler()]) from up to five tables. Each table can be a path to a CSV
#' file or a data frame typed in R; alternatively, `dir` names a folder
#' containing `diseases.csv`, `associations.csv`, `impacts.csv`,
#' `interactions.csv` and `uncertainty.csv` (only `diseases.csv` is
#' required). [cm_template()] writes an example set to start from.
#'
#' All tables are checked before anything is built, and every problem is
#' reported at once, with its table, row and column (row 1 is the first row
#' below the header). [cm_check_inputs()] runs the same checks without
#' stopping.
#'
#' @section Columns:
#' Column names are not case-sensitive. Optional columns can be left out or
#' left empty; unrecognised columns are ignored (with a note).
#'
#' **diseases**: `id`, `value` (required); `type` (`prevalence` (default),
#' `probability` or `incidence_rate`, converted with `1 - exp(-value)`),
#' `time_horizon`, `reference_population`, `source`.
#'
#' **associations**: `disease1`, `disease2` (required); `value`, `measure`
#' (`OR` (default), `RR`, `RD`, `cond_prob` (P(disease1 | disease2)), `phi`,
#' `table`, `independent` or `unknown`), `n11`, `n10`, `n01`, `n00` (counts,
#' for `measure = table`), `adjusted` (TRUE/FALSE), `adjusted_for`, `source`.
#' Pairs that are not listed are independent (or unknown, see
#' `missing_associations`).
#'
#' **impacts**: `disease`, `outcome`, `value` (required); `scale`
#' (`proportion` (default), `percent`, `absolute` or `hazard_ratio`),
#' `units` (required for `absolute`), `direction` (`decrease` (default) or
#' `increase`; ignored for hazard ratios), `adjusted_for` (ids separated by
#' `;`), `source`. Every disease needs one row per outcome (use 0 for no
#' impact).
#'
#' **interactions**: `disease1`, `disease2`, `value`, `outcome` (required);
#' `scale` (`proportion` (default) or `percent`), `source`.
#'
#' **uncertainty**: `key`, `dist` (required); `p1`-`p4` (parameters),
#' `note`. Keys name the input: `prob:<disease>`, `assoc:<d1>:<d2>`,
#' `impact:<outcome>:<disease>` or `inter:<outcome>:<d1>:<d2>`. Values are
#' drawn on the scale the input was entered on (e.g. an incidence rate, or a
#' percent impact). Distributions and their parameters:
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
#' @param diseases,associations,impacts,interactions,uncertainty Paths to CSV
#'   files or data frames. `NULL` tables are taken from `dir` (if given) or
#'   left out.
#' @param dir Optional folder with the CSV files named as above.
#' @param missing_associations Passed to [cm_model()].
#' @param outcome_correlation Passed to [cm_sampler()].
#' @return A `cm_inputs` list with `model`, `sampler` (or `NULL`), `tables`
#'   (the tables as read) and `problems` (notes only, since errors stop).
#' @export
#' @examples
#' dir <- file.path(tempdir(), "deconflate-inputs")
#' cm_template(dir, overwrite = TRUE)
#' inp <- cm_read_inputs(dir = dir)
#' inp
#' deconflate(inp$model, method = "global")
#'
#' # The same tables typed in R
#' inp2 <- cm_read_inputs(
#'   diseases = data.frame(id = c("d1", "d2"), value = c(0.10, 0.15)),
#'   associations = data.frame(disease1 = "d1", disease2 = "d2", value = 2),
#'   impacts = data.frame(disease = c("d1", "d2"), outcome = "yield",
#'                        value = c(2.5, 5), scale = "percent")
#' )
#' deconflate(inp2$model)
cm_read_inputs <- function(diseases = NULL, associations = NULL, impacts = NULL,
                           interactions = NULL, uncertainty = NULL, dir = NULL,
                           missing_associations = c("independent", "unknown"),
                           outcome_correlation = NULL) {
  missing_associations <- match.arg(missing_associations)
  chk <- run_input_checks(diseases, associations, impacts, interactions, uncertainty, dir,
                          missing_associations, outcome_correlation)
  pr <- chk$problems
  if (any(pr$severity == "error")) {
    cm_abort(format_problems(pr[pr$severity == "error", , drop = FALSE]),
             class = "deconflate_input_problems")
  }
  if (nrow(pr)) message(format_problems(pr))
  structure(list(model = chk$model, sampler = chk$sampler, tables = chk$tables,
                 problems = pr),
            class = "cm_inputs")
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
#'   impacts = data.frame(disease = c("d1", "d3"), outcome = "yield", value = c(2, "x"))
#' )
cm_check_inputs <- function(diseases = NULL, associations = NULL, impacts = NULL,
                            interactions = NULL, uncertainty = NULL, dir = NULL,
                            missing_associations = c("independent", "unknown"),
                            outcome_correlation = NULL) {
  missing_associations <- match.arg(missing_associations)
  run_input_checks(diseases, associations, impacts, interactions, uncertainty, dir,
                   missing_associations, outcome_correlation)$problems
}

#' Build distributions from an uncertainty table
#'
#' Converts a table with columns `key`, `dist` and `p1`-`p4` (see the
#' Columns section of [cm_read_inputs()]) into a named list of `cm_dist`
#' objects.
#'
#' @param x A path to a CSV file or a data frame.
#' @return A named list of `cm_dist` objects (names = keys).
#' @export
#' @examples
#' cm_dist_table(data.frame(key = c("impact:yield:d1", "assoc:d1:d2"),
#'                          dist = c("normal", "lognormal_ci"),
#'                          p1 = c(2.5, 2), p2 = c(0.5, 1.4), p3 = c(NA, 2.9)))
cm_dist_table <- function(x) {
  tab <- read_input_table(x)
  if (is.null(tab) || !nrow(tab)) cm_abort("The uncertainty table is empty.")
  if (!all(c("key", "dist") %in% names(tab))) cm_abort("The uncertainty table needs columns `key` and `dist`.")
  key <- clean_chr(tab$key)
  dist <- tolower(clean_chr(tab$dist))
  P <- vapply(paste0("p", 1:4), function(col) {
    v <- tab[[col]]
    if (is.null(v)) rep(NA_real_, nrow(tab)) else suppressWarnings(as.numeric(v))
  }, numeric(nrow(tab)))
  P <- matrix(P, nrow = nrow(tab))
  errs <- character(0)
  out <- list()
  for (r in seq_len(nrow(tab))) {
    if (is.na(key[r])) {
      errs <- c(errs, sprintf("row %d: missing key.", r))
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
#' Writes example `diseases.csv`, `associations.csv`, `impacts.csv`,
#' `interactions.csv` (header only) and `uncertainty.csv` files to a folder,
#' to edit and read back with [cm_read_inputs()]. The example has three dairy
#' diseases, with yield impacts and culling hazard ratios; the values are
#' illustrative.
#'
#' @param dir Folder to write to (created if needed).
#' @param overwrite Overwrite existing files?
#' @return The file paths, invisibly.
#' @export
#' @examples
#' cm_template(file.path(tempdir(), "my-inputs"), overwrite = TRUE)
cm_template <- function(dir, overwrite = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  tabs <- list(
    diseases = data.frame(
      id = c("LAM", "SCK", "MET"), value = c(0.25, 0.48, 0.10),
      type = c("prevalence", "incidence_rate", "prevalence"),
      time_horizon = "lactation", source = "Illustrative values",
      stringsAsFactors = FALSE),
    associations = data.frame(
      disease1 = c("LAM", "MET", "MET"), disease2 = c("SCK", "SCK", "LAM"),
      value = c(2.01, 1.94, 6.10), measure = "OR", source = "Illustrative values",
      stringsAsFactors = FALSE),
    impacts = data.frame(
      disease = rep(c("LAM", "SCK", "MET"), 2),
      outcome = rep(c("yield", "culling"), each = 3),
      value = c(4.81, 8.40, 5.61, 1.74, 1.92, 1.12),
      scale = rep(c("percent", "hazard_ratio"), each = 3),
      direction = rep(c("decrease", NA), each = 3),
      units = NA, adjusted_for = NA, source = "Illustrative values",
      stringsAsFactors = FALSE),
    interactions = data.frame(
      disease1 = character(0), disease2 = character(0), value = numeric(0),
      outcome = character(0), scale = character(0), source = character(0),
      stringsAsFactors = FALSE),
    uncertainty = data.frame(
      key = c("prob:LAM", "assoc:LAM:SCK", "impact:yield:LAM", "impact:yield:SCK",
              "impact:culling:LAM"),
      dist = c("beta", "normal", "normal", "normal", "normal"),
      p1 = c(78.29, 2.01, 4.81, 8.40, 1.74), p2 = c(227.42, 0.20, 0.87, 1.19, 0.17),
      p3 = c(NA, 0, NA, NA, 0), p4 = NA,
      note = c("beta(shape1, shape2)", "normal(mean, sd) truncated at 0",
               "normal(mean, sd), percent", "normal(mean, sd), percent",
               "normal(mean, sd) truncated at 0"),
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
  print(x$model)
  if (!is.null(x$sampler)) {
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
  impacts = list(required = c("disease", "outcome", "value"),
                 optional = c("scale", "units", "direction", "adjusted_for", "source")),
  interactions = list(required = c("disease1", "disease2", "value", "outcome"),
                      optional = c("scale", "source")),
  uncertainty = list(required = c("key", "dist"),
                     optional = c("p1", "p2", "p3", "p4", "note"))
)

# Read one table from a path or a data frame, with normalised column names.
read_input_table <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) cm_abort(sprintf("File not found: %s", x))
    x <- utils::read.csv(x, stringsAsFactors = FALSE, na.strings = c("", "NA"),
                         strip.white = TRUE, check.names = FALSE,
                         fileEncoding = "UTF-8-BOM")
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

run_input_checks <- function(diseases, associations, impacts, interactions, uncertainty, dir,
                             missing_associations, outcome_correlation) {
  probs <- list()
  add <- function(table, row, column, problem, severity = "error") {
    probs[[length(probs) + 1L]] <<- data.frame(
      table = table, row = as.integer(row), column = as.character(column),
      severity = severity, problem = problem, stringsAsFactors = FALSE)
  }
  args <- list(diseases = diseases, associations = associations, impacts = impacts,
               interactions = interactions, uncertainty = uncertainty)
  tabs <- list()
  for (nm in names(args)) {
    src <- args[[nm]]
    if (is.null(src) && !is.null(dir)) {
      f <- file.path(dir, paste0(nm, ".csv"))
      if (file.exists(f)) src <- f
    }
    tab <- tryCatch(read_input_table(src), error = function(e) e)
    if (inherits(tab, "condition")) {
      add(nm, NA, NA, conditionMessage(tab))
      tab <- NULL
    }
    if (!is.null(tab) && !nrow(tab)) tab <- NULL
    tabs[nm] <- list(tab)
  }
  if (!is.null(dir) && !dir.exists(dir)) add("inputs", NA, NA, sprintf("Folder not found: %s", dir))

  has_cols <- function(tab, nm) {
    sp <- input_table_specs[[nm]]
    miss <- setdiff(sp$required, names(tab))
    for (col in miss) add(nm, NA, col, "Required column is missing.")
    extra <- setdiff(names(tab), c(sp$required, sp$optional))
    for (col in extra) add(nm, NA, col, "Unrecognised column (ignored).", "note")
    !length(miss)
  }
  num <- function(tab, col, nm, required = TRUE) {
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
      for (r in which(bad)) add(nm, r, col, sprintf("'%s' is not a number.", vv[r]))
    }
    if (required) for (r in which(is.na(x) & !bad)) add(nm, r, col, "Missing value.")
    x
  }
  chr <- function(tab, col, default = NA_character_) clean_chr(tab[[col]], default, nrow(tab))
  check_ids <- function(v, nm, col, ids) {
    for (r in which(is.na(v))) add(nm, r, col, "Missing disease id.")
    if (!is.null(ids)) {
      for (r in which(!is.na(v) & !(v %in% ids))) {
        add(nm, r, col, sprintf("Unknown disease '%s' (not in the diseases table).", v[r]))
      }
    }
  }

  # Diseases ---------------------------------------------------------------
  ids <- NULL
  d <- tabs$diseases
  ok_d <- FALSE
  if (is.null(d)) {
    add("diseases", NA, NA, "A diseases table is required.")
  } else if (has_cols(d, "diseases")) {
    id <- chr(d, "id")
    for (r in which(is.na(id))) add("diseases", r, "id", "Missing disease id.")
    for (r in which(duplicated(id) & !is.na(id))) {
      add("diseases", r, "id", sprintf("Duplicate disease id '%s'.", id[r]))
    }
    for (r in which(!is.na(id) & grepl("[|;:]", id))) {
      add("diseases", r, "id", sprintf("Disease id '%s' must not contain '|', ';' or ':'.", id[r]))
    }
    d_value <- num(d, "value", "diseases")
    for (r in which(!is.na(d_value) & d_value < 0)) add("diseases", r, "value", "Must not be negative.")
    d_type <- chr(d, "type", "prevalence")
    bad_type <- !(d_type %in% c("prevalence", "probability", "incidence_rate"))
    for (r in which(bad_type)) {
      add("diseases", r, "type", sprintf("Unknown type '%s' (use prevalence, probability or incidence_rate).", d_type[r]))
    }
    prob <- ifelse(d_type == "incidence_rate", 1 - exp(-d_value), d_value)
    for (r in which(!bad_type & !is.na(prob) & d_value >= 0 & (prob <= 0 | prob >= 1))) {
      add("diseases", r, "value", sprintf("Gives a probability of %g; it must be strictly between 0 and 1.", prob[r]))
    }
    ids <- unique(id[!is.na(id)])
    ok_d <- TRUE
  }

  # Associations -----------------------------------------------------------
  a <- tabs$associations
  ok_a <- TRUE
  if (!is.null(a)) {
    ok_a <- has_cols(a, "associations")
    if (ok_a) {
      a1 <- chr(a, "disease1")
      a2 <- chr(a, "disease2")
      check_ids(a1, "associations", "disease1", ids)
      check_ids(a2, "associations", "disease2", ids)
      for (r in which(!is.na(a1) & !is.na(a2) & a1 == a2)) {
        add("associations", r, "disease2", "An association needs two different diseases.")
      }
      key <- ifelse(is.na(a1) | is.na(a2), NA, pair_key(a1, a2))
      for (r in which(duplicated(key) & !is.na(key))) {
        add("associations", r, NA, sprintf("Duplicate pair %s:%s.", a1[r], a2[r]))
      }
      a_measure <- chr(a, "measure", "OR")
      valid_m <- c("OR", "RR", "RD", "cond_prob", "phi", "table", "independent", "unknown")
      for (r in which(!(a_measure %in% valid_m))) {
        add("associations", r, "measure", sprintf("Unknown measure '%s' (use %s).", a_measure[r],
                                                  paste(valid_m, collapse = ", ")))
      }
      a_value <- num(a, "value", "associations", required = FALSE)
      needs <- a_measure %in% c("OR", "RR", "RD", "cond_prob", "phi")
      for (r in which(needs & is.na(a_value))) {
        add("associations", r, "value", sprintf("Missing value for measure %s.", a_measure[r]))
      }
      rng <- function(sel, lo, hi, txt) {
        for (r in which(sel & !is.na(a_value) & (a_value < lo | a_value > hi))) {
          add("associations", r, "value", txt)
        }
      }
      for (r in which(a_measure %in% c("OR", "RR") & !is.na(a_value) & a_value <= 0)) {
        add("associations", r, "value", "Odds ratios and risk ratios must be positive.")
      }
      rng(a_measure == "cond_prob", 0, 1, "A conditional probability must be between 0 and 1.")
      rng(a_measure %in% c("RD", "phi"), -1, 1, "Must be between -1 and 1.")
      cnt <- lapply(c("n11", "n10", "n01", "n00"), function(col) num(a, col, "associations", FALSE))
      cnt <- do.call(cbind, cnt)
      for (r in which(a_measure == "table")) {
        if (anyNA(cnt[r, ]) || any(cnt[r, ] < 0)) {
          add("associations", r, NA, "measure = table needs non-negative counts n11, n10, n01 and n00.")
        }
      }
      a_adj <- chr(a, "adjusted", "FALSE")
      a_adj_l <- toupper(a_adj) %in% c("TRUE", "T", "YES", "1")
      for (r in which(!(toupper(a_adj) %in% c("TRUE", "T", "YES", "1", "FALSE", "F", "NO", "0")))) {
        add("associations", r, "adjusted", sprintf("'%s' is not TRUE or FALSE.", a_adj[r]))
      }
    }
  }

  # Impacts ----------------------------------------------------------------
  im <- tabs$impacts
  ok_i <- TRUE
  outcomes <- character(0)
  out_scale <- character(0)
  if (!is.null(im)) {
    ok_i <- has_cols(im, "impacts")
    if (ok_i) {
      i_dis <- chr(im, "disease")
      i_out <- chr(im, "outcome")
      i_value <- num(im, "value", "impacts")
      check_ids(i_dis, "impacts", "disease", ids)
      for (r in which(is.na(i_out))) add("impacts", r, "outcome", "Missing outcome.")
      i_scale <- chr(im, "scale", "proportion")
      valid_s <- c("proportion", "percent", "absolute", "hazard_ratio")
      for (r in which(!(i_scale %in% valid_s))) {
        add("impacts", r, "scale", sprintf("Unknown scale '%s' (use %s).", i_scale[r],
                                           paste(valid_s, collapse = ", ")))
      }
      i_dir <- chr(im, "direction", "decrease")
      for (r in which(!(i_dir %in% c("decrease", "increase")) & i_scale != "hazard_ratio")) {
        add("impacts", r, "direction", sprintf("Unknown direction '%s' (use decrease or increase).", i_dir[r]))
      }
      i_dir[i_scale == "hazard_ratio"] <- "increase"
      i_units <- chr(im, "units")
      for (r in which(i_scale == "absolute" & is.na(i_units))) {
        add("impacts", r, "units", "Units are required on the absolute scale.")
      }
      for (r in which(i_scale == "hazard_ratio" & !is.na(i_value) & i_value <= 0)) {
        add("impacts", r, "value", "Hazard ratios must be positive.")
      }
      ikey <- paste(i_out, i_dis, sep = ":")
      for (r in which(duplicated(ikey) & !is.na(i_out) & !is.na(i_dis))) {
        add("impacts", r, NA, sprintf("Disease '%s' has more than one impact on outcome '%s'.", i_dis[r], i_out[r]))
      }
      i_adj <- chr(im, "adjusted_for")
      for (r in which(!is.na(i_adj))) {
        unk <- setdiff(split_ids(i_adj[r]), ids)
        if (length(unk) && !is.null(ids)) {
          add("impacts", r, "adjusted_for", sprintf("Unknown disease(s): %s.", paste(unk, collapse = ", ")))
        }
      }
      outcomes <- unique(i_out[!is.na(i_out)])
      for (o in outcomes) {
        sel <- !is.na(i_out) & i_out == o
        sc <- unique(ifelse(i_scale[sel] == "percent", "proportion", i_scale[sel]))
        out_scale[o] <- sc[1]
        if (length(sc) > 1L) {
          add("impacts", NA, "scale", sprintf("Outcome '%s' mixes scales: %s.", o, paste(sc, collapse = ", ")))
        }
        if (!identical(sc[1], "hazard_ratio") && length(unique(i_dir[sel])) > 1L) {
          add("impacts", NA, "direction", sprintf("Outcome '%s' mixes directions.", o))
        }
        un <- unique(i_units[sel])
        if (length(un) > 1L) {
          add("impacts", NA, "units", sprintf("Outcome '%s' mixes units: %s.", o, paste(un, collapse = ", ")))
        }
        if (!is.null(ids)) {
          miss <- setdiff(ids, i_dis[sel])
          if (length(miss)) {
            add("impacts", NA, NA, sprintf("Outcome '%s' has no impact for: %s. Add a row with value 0 for no impact.",
                                           o, paste(miss, collapse = ", ")))
          }
        }
      }
    }
  }

  # Interactions -----------------------------------------------------------
  it <- tabs$interactions
  ok_x <- TRUE
  if (!is.null(it)) {
    ok_x <- has_cols(it, "interactions")
    if (ok_x) {
      x1 <- chr(it, "disease1")
      x2 <- chr(it, "disease2")
      check_ids(x1, "interactions", "disease1", ids)
      check_ids(x2, "interactions", "disease2", ids)
      for (r in which(!is.na(x1) & !is.na(x2) & x1 == x2)) {
        add("interactions", r, "disease2", "An interaction needs two different diseases.")
      }
      x_value <- num(it, "value", "interactions")
      x_out <- chr(it, "outcome")
      for (r in which(is.na(x_out))) add("interactions", r, "outcome", "Missing outcome.")
      for (r in which(!is.na(x_out) & !(x_out %in% outcomes))) {
        add("interactions", r, "outcome", sprintf("Outcome '%s' is not in the impacts table.", x_out[r]))
      }
      for (r in which(!is.na(x_out) & x_out %in% outcomes & out_scale[x_out] != "proportion")) {
        add("interactions", r, "outcome", sprintf("Interactions need an outcome on the proportion or percent scale ('%s' is not).", x_out[r]))
      }
      x_scale <- chr(it, "scale", "proportion")
      for (r in which(!(x_scale %in% c("proportion", "percent")))) {
        add("interactions", r, "scale", sprintf("Unknown scale '%s' (use proportion or percent).", x_scale[r]))
      }
    }
  }

  # Uncertainty ------------------------------------------------------------
  un <- tabs$uncertainty
  dists <- list()
  if (!is.null(un) && has_cols(un, "uncertainty")) {
    u_key <- chr(un, "key")
    u_dist <- tolower(chr(un, "dist"))
    P <- vapply(paste0("p", 1:4), function(col) num(un, col, "uncertainty", FALSE), numeric(nrow(un)))
    P <- matrix(P, nrow = nrow(un))
    assoc_keys <- if (!is.null(a) && ok_a) {
      ifelse(is.na(a1) | is.na(a2), NA, pair_key(a1, a2))[a_measure %in% c("OR", "RR", "RD", "cond_prob", "phi")]
    } else character(0)
    imp_keys <- if (!is.null(im) && ok_i) paste(i_out, i_dis, sep = ":") else character(0)
    int_keys <- if (!is.null(it) && ok_x) paste(x_out, pair_key(x1, x2), sep = ":") else character(0)
    for (r in which(is.na(u_key))) add("uncertainty", r, "key", "Missing key.")
    norm_key <- vapply(u_key, function(k) {
      if (is.na(k)) return(NA_character_)
      parts <- strsplit(k, ":", fixed = TRUE)[[1]]
      if (parts[1] == "assoc" && length(parts) == 3L) return(paste0("assoc:", pair_key(parts[2], parts[3])))
      if (parts[1] == "inter" && length(parts) == 4L) {
        return(paste0("inter:", parts[2], ":", pair_key(parts[3], parts[4])))
      }
      k
    }, character(1), USE.NAMES = FALSE)
    for (r in which(duplicated(norm_key) & !is.na(norm_key))) {
      add("uncertainty", r, "key", sprintf("Duplicate key '%s' (the same input appears twice).", u_key[r]))
    }
    for (r in which(!is.na(u_key))) {
      parts <- strsplit(u_key[r], ":", fixed = TRUE)[[1]]
      found <- switch(parts[1],
        prob = length(parts) == 2L && !is.null(ids) && parts[2] %in% ids,
        assoc = length(parts) == 3L && pair_key(parts[2], parts[3]) %in% assoc_keys,
        impact = length(parts) == 3L && paste(parts[2], parts[3], sep = ":") %in% imp_keys,
        inter = length(parts) == 4L &&
          paste(parts[2], pair_key(parts[3], parts[4]), sep = ":") %in% int_keys,
        NA)
      if (is.na(found)) {
        add("uncertainty", r, "key", sprintf("Key '%s' must start with prob:, assoc:, impact: or inter:.", u_key[r]))
      } else if (!found) {
        add("uncertainty", r, "key", sprintf("Key '%s' does not match an input (expected prob:<disease>, assoc:<d1>:<d2> with a numeric measure, impact:<outcome>:<disease> or inter:<outcome>:<d1>:<d2>).", u_key[r]))
      }
      dd <- tryCatch(input_dist(u_dist[r], P[r, ]), error = function(e) e)
      if (inherits(dd, "condition")) {
        add("uncertainty", r, "dist", conditionMessage(dd))
      } else if (isTRUE(found)) {
        dists[[u_key[r]]] <- dd
      }
    }
  }

  # Build ------------------------------------------------------------------
  pr <- if (length(probs)) do.call(rbind, probs) else {
    data.frame(table = character(0), row = integer(0), column = character(0),
               severity = character(0), problem = character(0), stringsAsFactors = FALSE)
  }
  model <- NULL
  sampler <- NULL
  if (ok_d && !any(pr$severity == "error")) {
    built <- tryCatch({
      dis_obj <- cm_diseases(id, d_value, type = d_type,
                             time_horizon = chr(d, "time_horizon"),
                             reference_population = chr(d, "reference_population"),
                             source = chr(d, "source"))
      assoc_obj <- if (!is.null(a)) {
        cm_associations(a1, a2, a_value, measure = a_measure,
                        n11 = cnt[, 1], n10 = cnt[, 2], n01 = cnt[, 3], n00 = cnt[, 4],
                        adjusted = a_adj_l, adjusted_for = chr(a, "adjusted_for"),
                        source = chr(a, "source"))
      } else NULL
      imp_obj <- if (!is.null(im)) {
        cm_impacts(i_dis, i_value, outcome = i_out, scale = i_scale, units = i_units,
                   direction = i_dir, adjusted_for = i_adj, source = chr(im, "source"))
      } else NULL
      int_obj <- if (!is.null(it)) {
        cm_interactions(x1, x2, x_value, outcome = x_out, scale = x_scale,
                        source = chr(it, "source"))
      } else NULL
      m <- cm_model(dis_obj, assoc_obj, imp_obj, int_obj,
                    missing_associations = missing_associations)
      pair_tables(m)
      m
    }, error = function(e) e)
    if (inherits(built, "condition")) {
      add("model", NA, NA, conditionMessage(built))
    } else {
      model <- built
      if (length(dists)) {
        pick <- function(prefix) {
          x <- dists[startsWith(names(dists), prefix)]
          stats::setNames(x, substring(names(x), nchar(prefix) + 1L))
        }
        sampler <- tryCatch(
          cm_sampler(model, diseases = pick("prob:"), associations = pick("assoc:"),
                     impacts = pick("impact:"), interactions = pick("inter:"),
                     outcome_correlation = outcome_correlation),
          error = function(e) e)
        if (inherits(sampler, "condition")) {
          add("uncertainty", NA, NA, conditionMessage(sampler))
          sampler <- NULL
        }
      }
    }
    pr <- do.call(rbind, probs)
  }
  if (is.null(pr)) {
    pr <- data.frame(table = character(0), row = integer(0), column = character(0),
                     severity = character(0), problem = character(0), stringsAsFactors = FALSE)
  }
  rownames(pr) <- NULL
  class(pr) <- c("cm_problems", "data.frame")
  list(model = model, sampler = sampler, tables = tabs, problems = pr)
}
