#' Compare adjustment methods
#'
#' Runs several adjustment methods on the same inputs and tabulates the
#' results side by side:
#' * `"published"`: the simple proportional approximation used in Rasmussen
#'   et al. (2022, 2024) (eq. 16; for hazard ratios, `HR - 1` as in 2024);
#' * `"simultaneous"`: the exact solution of the system of equations;
#' * `"global"`: the iterative (maximum-entropy) model of disease
#'   combinations, which can also include interactions.
#'
#' For a model, every method is run on the model's inputs. For a Monte Carlo
#' run made with several methods (`cm_monte_carlo(..., method = c(...))`),
#' the methods were applied to identical draws, and the table compares their
#' summaries.
#'
#' @param x A [cm_model()], or a `cm_mc` object from [cm_monte_carlo()] run
#'   with more than one method.
#' @param methods Methods to compare (models only).
#' @param economics Optional economics list (see [value_losses()]; for
#'   hazard-ratio outcomes, [attributable_risk()]). Adds gaps and values per
#'   method to `totals`.
#' @param stat For Monte Carlo runs: the statistic to compare, `"mean"`,
#'   `"median"` or `"trimmed_mean"`.
#' @param ... Passed to [deconflate()] (models) or [summary.cm_mc()] (Monte
#'   Carlo runs).
#' @return A `cm_comparison` object with:
#'   * `impacts`: one row per outcome and disease, with the raw impact and
#'     one column of adjusted impacts per method;
#'   * `change`: the same layout with the relative change from the raw
#'     impact;
#'   * `long`: the same information in long format, with sign-change flags
#'     (Monte Carlo: the summary statistics and stability flags per method);
#'   * `totals`: per outcome and method, the expected loss before and after
#'     adjustment and, with `economics`, the gap and its value;
#'   * `total_value`: total value per method (with `economics`);
#'   * `diagnostics`: [deconflate()] diagnostics per method (models);
#'   * `failed`: methods that could not be run, with the reason.
#' @export
#' @examples
#' cmp <- compare_methods(example_uk_dairy_2022(), economics = uk_dairy_2022_economics())
#' cmp
#' cmp$totals
compare_methods <- function(x, ...) UseMethod("compare_methods")

#' @rdname compare_methods
#' @export
compare_methods.cm_model <- function(x, methods = c("published", "simultaneous", "global"),
                                     economics = NULL, ...) {
  methods <- unique(match.arg(methods, c("published", "simultaneous", "global"),
                              several.ok = TRUE))
  res <- list()
  failed <- character(0)
  for (m in methods) {
    r <- tryCatch(
      withCallingHandlers(deconflate(x, method = m, warn = FALSE, ...),
                          deconflate_nonconvergence = function(w) invokeRestart("muffleWarning")),
      deconflate_error = function(e) e
    )
    if (inherits(r, "condition")) failed[m] <- conditionMessage(r) else res[[m]] <- r
  }
  if (!length(res)) {
    cm_abort(sprintf("All methods failed: %s", paste(names(failed), failed, sep = ": ", collapse = "; ")))
  }
  ok <- names(res)
  base <- res[[1]]$adjusted
  impacts <- base[, c("outcome", "disease", "scale", "raw")]
  change <- impacts
  for (m in ok) {
    impacts[[m]] <- res[[m]]$adjusted$adjusted
    change[[m]] <- res[[m]]$adjusted$change
  }
  long <- do.call(rbind, lapply(ok, function(m) {
    a <- res[[m]]$adjusted
    hr <- a$scale == "hazard_ratio"
    ex_raw <- ifelse(hr, a$raw - 1, a$raw)
    ex_adj <- ifelse(hr, a$adjusted - 1, a$adjusted)
    data.frame(outcome = a$outcome, disease = a$disease, method = m, raw = a$raw,
               adjusted = a$adjusted, change = a$change,
               sign_change = is.finite(ex_adj) & abs(ex_raw) > 1e-12 & abs(ex_adj) > 1e-12 &
                 sign(ex_adj) != sign(ex_raw),
               stringsAsFactors = FALSE)
  }))
  P <- stats::setNames(x$diseases$prob, x$diseases$id)
  total_value <- NULL
  totals <- do.call(rbind, lapply(ok, function(m) {
    r <- res[[m]]
    b <- attribute_burden(r)
    add <- r$adjusted[r$adjusted$scale != "hazard_ratio", , drop = FALSE]
    tt <- if (nrow(add)) {
      do.call(rbind, lapply(split(add, add$outcome), function(a) {
        data.frame(outcome = a$outcome[1], method = m,
                   raw_loss = sum(a$raw * P[a$disease]),
                   adjusted_loss = sum(b$total[b$outcome == a$outcome[1]]),
                   stringsAsFactors = FALSE)
      }))
    } else {
      data.frame(outcome = character(0), method = character(0), raw_loss = numeric(0),
                 adjusted_loss = numeric(0), stringsAsFactors = FALSE)
    }
    if (!is.null(economics)) {
      ev <- evaluate_economics(r, economics, allocate = FALSE)
      ev_tab <- merge(ev$summary, ev$by_outcome, by = "outcome", sort = FALSE)
      ev_tab$method <- m
      tt <- merge(tt, ev_tab, by = c("outcome", "method"), all = TRUE, sort = FALSE)
      total_value <<- c(total_value, stats::setNames(ev$total, m))
    }
    tt
  }))
  rownames(totals) <- NULL
  diagnostics <- do.call(rbind, lapply(ok, function(m) res[[m]]$diagnostics))
  rownames(diagnostics) <- NULL
  structure(list(impacts = impacts, change = change, long = long, totals = totals,
                 total_value = total_value, diagnostics = diagnostics, failed = failed,
                 methods = ok, source = "model", results = res),
            class = "cm_comparison")
}

#' @rdname compare_methods
#' @export
compare_methods.cm_mc <- function(x, stat = c("mean", "median", "trimmed_mean"), ...) {
  stat <- match.arg(stat)
  if (length(x$method) < 2L) {
    cm_abort("This Monte Carlo run used one method; run cm_monte_carlo(..., method = c(\"published\", \"simultaneous\")) to compare methods.")
  }
  s <- summary(x, diagnose = FALSE, ...)
  col <- switch(stat, mean = "mean", median = "q0.5", trimmed_mean = "trimmed_mean")
  if (is.null(s[[col]])) cm_abort("The median needs probs to include 0.5.")
  d <- x$draws
  if (is.null(d$scale)) d$scale <- NA_character_
  key <- unique(d[, c("outcome", "disease", "scale")])
  rownames(key) <- NULL
  raw <- vapply(seq_len(nrow(key)), function(i) {
    sel <- d$outcome == key$outcome[i] & d$disease == key$disease[i] & d$method == x$method[1]
    w <- x$weights[match(d$draw[sel], x$params$draw)]
    sum(w * d$raw[sel]) / sum(w)
  }, numeric(1))
  impacts <- data.frame(key, raw_mean = raw, stringsAsFactors = FALSE)
  stability <- key[, c("outcome", "disease")]
  for (m in x$method) {
    sm <- s[s$method == m, , drop = FALSE]
    idx <- match(paste(key$outcome, key$disease), paste(sm$outcome, sm$disease))
    impacts[[m]] <- sm[[col]][idx]
    stability[[m]] <- sm$stability[idx]
  }
  rownames(impacts) <- NULL
  totals <- NULL
  total_value <- NULL
  if (!is.null(x$losses)) {
    tt <- summary(x, what = "total", diagnose = FALSE)
    totals <- tt[, c("outcome", "method", "mean", "q0.5", "trimmed_mean", "stability")]
    tv <- tt[tt$outcome == "total", , drop = FALSE]
    total_value <- stats::setNames(tv[[col]], tv$method)
  }
  structure(list(impacts = impacts, stability = stability, long = s, totals = totals,
                 total_value = total_value, failed = character(0), methods = x$method,
                 source = "mc", stat = stat, n_draws = x$n_draws, n_rejected = x$n_rejected),
            class = "cm_comparison")
}

#' @export
print.cm_comparison <- function(x, digits = 3, ...) {
  cat(sprintf("<cm_comparison> methods: %s\n", paste(x$methods, collapse = ", ")))
  if (identical(x$source, "mc")) {
    cat(sprintf("Monte Carlo: %d draws (%d rejected); statistic: %s\n",
                x$n_draws, x$n_rejected, x$stat))
  }
  tab <- x$impacts
  num <- vapply(tab, is.numeric, logical(1))
  if (!is.null(tab$scale)) {
    pct <- !is.na(tab$scale) & tab$scale == "proportion"
    tab[pct, num] <- tab[pct, num] * 100
    tab$scale <- ifelse(pct, "%", ifelse(!is.na(tab$scale) & tab$scale == "hazard_ratio", "HR", tab$scale))
    names(tab)[names(tab) == "scale"] <- "unit"
  }
  tab[num] <- lapply(tab[num], signif, digits = digits)
  cat("\nAdjusted impacts:\n")
  print(tab, row.names = FALSE)
  if (identical(x$source, "model")) {
    sc <- x$long[x$long$sign_change, , drop = FALSE]
    if (nrow(sc)) {
      cat("\nSign changes (adjusted impact on the other side of zero, or of 1 for hazard ratios):\n")
      for (m in unique(sc$method)) {
        s <- sc[sc$method == m, , drop = FALSE]
        cat(sprintf("  %s: %s\n", m, paste(sprintf("%s %s", s$outcome, s$disease), collapse = ", ")))
      }
    }
  } else {
    st <- x$stability
    flagged <- vapply(x$methods, function(m) sum(st[[m]] != "ok", na.rm = TRUE), numeric(1))
    if (any(flagged > 0)) {
      cat(sprintf("\nUnstable estimates per method: %s (see cm_diagnose()).\n",
                  paste(sprintf("%s %d", x$methods, flagged), collapse = ", ")))
    }
  }
  if (!is.null(x$totals) && nrow(x$totals)) {
    cat("\nTotals:\n")
    tt <- x$totals
    nn <- vapply(tt, is.numeric, logical(1))
    tt[nn] <- lapply(tt[nn], signif, digits = 4)
    print(tt, row.names = FALSE)
  }
  if (!is.null(x$total_value)) {
    cat(sprintf("\nTotal value: %s\n",
                paste(sprintf("%s %.2f", names(x$total_value), x$total_value), collapse = "; ")))
  }
  if (length(x$failed)) {
    cat("\nFailed:\n")
    for (m in names(x$failed)) cat(sprintf("  %s: %s\n", m, x$failed[[m]]))
  }
  invisible(x)
}
