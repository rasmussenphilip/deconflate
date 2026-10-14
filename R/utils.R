# Internal helpers -----------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Raise a classed error. All package errors inherit from "deconflate_error".
# Classes used for failures that callers (draws, comparisons, screens)
# handle separately:
#   deconflate_infeasible      inputs that admit no valid probabilities
#   deconflate_singular        a singular or non-identifiable system
#   deconflate_nonconvergence  a numerical procedure did not converge
#   deconflate_unsupported     an estimand/method combination that is not supported
#   deconflate_nonfinite       non-finite results
#   deconflate_invalid_input   an argument, table, row or input value that is not valid
# Two more mark particular errors:
#   deconflate_unknown_pairs   (with deconflate_unsupported) a calculation needs
#                              a crude pair row for every pair
#   deconflate_input_problems  cm_read_inputs() found errors in the tables
# Warnings inherit from "deconflate_warning"; their classes are
# deconflate_sign_change, deconflate_nonfinite_warning,
# deconflate_nonconvergence and deconflate_incomplete_allocation. Users find
# the classes in ?deconflate-package (Conventions); keep both lists in step.
# `...`: named fields of the condition (e.g. `rows`, `terms`), so that callers
# can use what the message describes without parsing it.
cm_abort <- function(msg, class = NULL, ...) {
  stop(errorCondition(msg, ..., class = c(class, "deconflate_error"), call = NULL))
}

cm_warn <- function(msg, class = NULL, ...) {
  warning(warningCondition(msg, ..., class = c(class, "deconflate_warning"),
                           call = NULL))
}

# Short labels of the condition classes above.
condition_types <- c("infeasible", "singular", "nonconvergence", "unsupported", "nonfinite", "invalid_input")

# Short label for the class of a deconflate condition.
condition_type <- function(e) {
  hit <- condition_types[paste0("deconflate_", condition_types) %in% class(e)]
  if (length(hit)) hit[1] else "error"
}

# The class of a short label (NULL for "error" and other labels), to raise a
# condition of the same kind as the one that caused it.
type_class <- function(type) {
  if (length(type) == 1L && type %in% condition_types) paste0("deconflate_", type)
}

check_numeric <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    cm_abort(sprintf("`%s` must be numeric, finite and without missing values.",
                     name), class = "deconflate_invalid_input")
  }
  invisible(x)
}

check_scalar_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) != 1L || !(x %in% choices)) {
    cm_abort(sprintf("`%s` must be one of: %s.", name,
                     paste(sprintf("'%s'", choices), collapse = ", ")), class = "deconflate_invalid_input")
  }
  invisible(x)
}

check_choices <- function(x, choices, name) {
  bad <- setdiff(unique(x), choices)
  if (length(bad)) {
    cm_abort(sprintf("Invalid `%s`: %s. Allowed: %s.", name,
                     paste(sprintf("'%s'", bad), collapse = ", "),
                     paste(sprintf("'%s'", choices), collapse = ", ")), class = "deconflate_invalid_input")
  }
  invisible(x)
}

# match.arg() for the arguments of exported functions, with an invalid choice
# as a classed error that names the argument. Without `choices`, they are the
# argument's default in the calling function, as for match.arg().
cm_match_arg <- function(arg, choices, several_ok = FALSE, name = deparse(substitute(arg))) {
  if (missing(choices)) choices <- eval(formals(sys.function(sys.parent()))[[name]], envir = parent.frame())
  tryCatch(match.arg(arg, choices, several.ok = several_ok), error = function(e) {
    cm_abort(sprintf("`%s` must be %s: %s.", name, if (several_ok) "one or more of" else "one of",
                     paste(sprintf("'%s'", choices), collapse = ", ")), class = "deconflate_invalid_input")
  })
}

# Disease ids are used in keys such as "assoc:d1:d2" and in "a; b" lists.
check_ids <- function(id) {
  if (anyNA(id) || any(!nzchar(id))) cm_abort("Disease ids must be non-empty.", class = "deconflate_invalid_input")
  if (anyDuplicated(id)) cm_abort("Disease ids must be unique.", class = "deconflate_invalid_input")
  bad <- grepl("[|;:]", id) | id != trimws(id)
  if (any(bad)) {
    cm_abort(sprintf("Disease ids must not contain '|', ';' or ':' or surrounding spaces (check: %s).",
                     paste(id[bad], collapse = ", ")), class = "deconflate_invalid_input")
  }
  if (any(tolower(id) %in% c("all", "others"))) cm_abort("'all' and 'others' are reserved and cannot be disease ids.", class = "deconflate_invalid_input")
  if (any(grepl("^no\\s", id, ignore.case = TRUE))) {
    cm_abort("Disease ids must not start with 'no ' (it marks absent diseases in `among`).", class = "deconflate_invalid_input")
  }
  invisible(id)
}

# Is an adjustment result usable as an estimate? The published approximation
# can return infinite or undefined values (a zero denominator); such a result
# is kept for inspection but never counts as a successful estimate. Used by
# Monte Carlo, method comparisons, sensitivity screens and threshold searches.
result_is_finite <- function(r) {
  if (inherits(r, "cm_event_result")) {
    a <- r$adjusted$adjusted
    ok <- length(a) > 0 && all(is.finite(a)) && all(a > 0)
    if (ok && !is.null(r$attributable)) {
      ok <- is.finite(r$attributable$summary$attributable) &&
        (is.null(r$attributable$by_disease) || all(is.finite(r$attributable$by_disease$attributable)))
    }
    return(ok)
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
    cm_abort(sprintf("`%s` must have length 1 or %d.", name, n), class = "deconflate_invalid_input")
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

# Method argument of the exported adjustment functions: the published
# approximation is no longer one of their methods.
public_method <- function(method, choices, fun) {
  if (is.character(method) && "published" %in% method && !identical(method, choices)) {
    cm_abort(sprintf("The published approximation is not a method of %s. It is kept for comparison and reproduction: use compare_methods() or the reproduce_*() functions.",
                     fun), class = "deconflate_unsupported")
  }
  match.arg(method, choices)
}

# Check the overall risk of an event model: a proportion in (0, 1), or a
# distribution (cm_dist) whose support lies in (0, 1). Returns the central
# value (the distribution's mean) and the distribution (or NULL).
check_overall_risk <- function(overall_risk, event_model) {
  if (!isTRUE(event_model)) {
    if (!is.null(overall_risk)) {
      cm_abort("`overall_risk` is used only with event_model = TRUE.", class = "deconflate_unsupported")
    }
    return(list(value = NULL, dist = NULL))
  }
  if (is.null(overall_risk)) {
    cm_abort("event_model = TRUE needs `overall_risk`: the overall risk of the event in the population over the period of the estimates (a proportion, e.g. 0.25, or a distribution such as dist_beta(250, 750)).", class = "deconflate_invalid_input")
  }
  if (inherits(overall_risk, "cm_dist")) {
    sp <- dist_support(overall_risk)
    if (min(sp[, 1]) < 0 || max(sp[, 2]) > 1) {
      cm_abort(sprintf("The distribution of `overall_risk` must lie between 0 and 1 (its support is %s).",
                       format_support(overall_risk)), class = "deconflate_invalid_input")
    }
    v <- overall_risk$mean
    if (!is.finite(v) || v <= 0 || v >= 1) cm_abort("The mean of the `overall_risk` distribution must lie strictly between 0 and 1.", class = "deconflate_invalid_input")
    return(list(value = v, dist = overall_risk))
  }
  if (!is.numeric(overall_risk) || length(overall_risk) != 1L || !is.finite(overall_risk) ||
      overall_risk <= 0 || overall_risk >= 1) {
    cm_abort("`overall_risk` must be a single proportion strictly between 0 and 1, or a distribution.", class = "deconflate_invalid_input")
  }
  list(value = overall_risk, dist = NULL)
}

# Show a message of class "deconflate_message" (and `class`), so that callers
# can muffle it by class; `...`: named fields of the condition. Built by hand,
# as messageCondition() needs R 4.6; message() adds no newline to a condition.
cm_inform <- function(msg, class = NULL, ...) {
  message(structure(class = c(class, "deconflate_message", "message", "condition"),
                    list(message = paste0(msg, "\n"), call = NULL, ...)))
}
