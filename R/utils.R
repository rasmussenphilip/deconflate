# Internal helpers -----------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Raise a classed error. All package errors inherit from "deconflate_error";
# infeasible inputs additionally carry "deconflate_infeasible" so that
# cm_monte_carlo() can count them as rejected draws.
cm_abort <- function(msg, class = NULL) {
  stop(errorCondition(msg, class = c(class, "deconflate_error"), call = NULL))
}

cm_warn <- function(msg, class = NULL) {
  warning(warningCondition(msg, class = c(class, "deconflate_warning"),
                           call = NULL))
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
  if (length(x) == 0L || is.na(x) || !nzchar(trimws(x))) return(character(0))
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
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
