# Internal helpers -----------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Raise a classed error. All package errors inherit from "deconflate_error".
# Classes used for failures that callers (Monte Carlo, comparisons, screens)
# handle separately:
#   deconflate_infeasible      inputs that admit no valid probabilities
#   deconflate_singular        a singular or non-identifiable system
#   deconflate_nonconvergence  a numerical procedure did not converge
#   deconflate_unsupported     an estimand/method combination that is not supported
#   deconflate_nonfinite       non-finite results
cm_abort <- function(msg, class = NULL) {
  stop(errorCondition(msg, class = c(class, "deconflate_error"), call = NULL))
}

cm_warn <- function(msg, class = NULL) {
  warning(warningCondition(msg, class = c(class, "deconflate_warning"),
                           call = NULL))
}

# Short label for the class of a deconflate condition.
condition_type <- function(e) {
  cls <- c("deconflate_infeasible", "deconflate_singular", "deconflate_nonconvergence",
           "deconflate_unsupported", "deconflate_nonfinite")
  hit <- cls[cls %in% class(e)]
  if (length(hit)) sub("deconflate_", "", hit[1]) else "error"
}

check_numeric <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    cm_abort(sprintf("`%s` must be numeric, finite and without missing values.",
                     name))
  }
  invisible(x)
}

check_scalar_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || !(x %in% choices)) {
    cm_abort(sprintf("`%s` must be one of: %s.", name,
                     paste(sprintf("'%s'", choices), collapse = ", ")))
  }
  invisible(x)
}

check_choices <- function(x, choices, name) {
  bad <- setdiff(unique(x), choices)
  if (length(bad)) {
    cm_abort(sprintf("Invalid `%s`: %s. Allowed: %s.", name,
                     paste(sprintf("'%s'", bad), collapse = ", "),
                     paste(sprintf("'%s'", choices), collapse = ", ")))
  }
  invisible(x)
}

# Disease ids are used in keys such as "assoc:d1:d2" and in "a; b" lists.
check_ids <- function(id) {
  if (anyNA(id) || any(!nzchar(id))) cm_abort("Disease ids must be non-empty.")
  if (anyDuplicated(id)) cm_abort("Disease ids must be unique.")
  bad <- grepl("[|;:]", id) | id != trimws(id)
  if (any(bad)) {
    cm_abort(sprintf("Disease ids must not contain '|', ';' or ':' or surrounding spaces (check: %s).",
                     paste(id[bad], collapse = ", ")))
  }
  if (any(tolower(id) == "all")) cm_abort("'all' is reserved and cannot be a disease id.")
  invisible(id)
}

# Is an adjustment result usable as an estimate? The published approximation
# can return infinite or undefined values (a zero denominator); such a result
# is kept for inspection but never counts as a successful estimate. Used by
# Monte Carlo, method comparisons, sensitivity screens and threshold searches.
result_is_finite <- function(r) {
  if (inherits(r, "cm_hr_result")) {
    a <- r$adjusted$adjusted
    return(length(a) > 0 && all(is.finite(a)) && all(a > 0))
  }
  inherits(r, "cm_result") && all(is.finite(r$adjusted$adjusted)) &&
    all(is.finite(r$contributions$total)) && is.finite(r$totals$adjusted_total)
}

# Names of the diseases with non-finite adjusted values.
nonfinite_diseases <- function(r) {
  a <- r$adjusted
  paste(a$disease[!is.finite(a$adjusted)], collapse = ", ")
}

# Recycle an optional per-row argument to length n.
recycle_arg <- function(x, n, name) {
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) != n) {
    cm_abort(sprintf("`%s` must have length 1 or %d.", name, n))
  }
  x
}

# Unordered pair key, e.g. pair_key("b", "a") == "a|b".
pair_key <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = "|"), paste(b, a, sep = "|"))
}

# Split "a; b; c" into c("a", "b", "c"); NA or "" gives character(0).
split_ids <- function(x) {
  if (length(x) != 1L) {
    if (length(x) == 0L) return(character(0))
    cm_abort("split_ids() takes a single string.")
  }
  if (is.na(x) || !nzchar(trimws(x))) return(character(0))
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

# Resolve an adjustment set for disease `i`: "all" means every other disease.
resolve_adjusted_for <- function(x, i, ids) {
  s <- split_ids(x)
  if (length(s) == 1L && tolower(s) == "all") return(setdiff(ids, i))
  setdiff(s, i)
}

# All 2^n combinations of n binary disease indicators (rows = combinations).
disease_cells <- function(n) {
  idx <- seq_len(2^n) - 1
  m <- vapply(seq_len(n), function(j) (idx %/% 2^(j - 1)) %% 2,
              numeric(length(idx)))
  matrix(m, ncol = n)
}

# Crude difference E[f | D_i = 1] - E[f | D_i = 0] for every disease i, given
# the joint distribution over cells.
crude_difference <- function(cells, prob, f) {
  vapply(seq_len(ncol(cells)), function(i) {
    w1 <- prob * cells[, i]
    w0 <- prob * (1 - cells[, i])
    sum(w1 * f) / sum(w1) - sum(w0 * f) / sum(w0)
  }, numeric(1))
}

weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]
  w <- w[o] / sum(w)
  cw <- cumsum(w)
  vapply(probs, function(p) x[which(cw >= p - 1e-12)[1]], numeric(1))
}
