#' Compare adjustment methods
#'
#' Runs several adjustment methods on the same inputs and tabulates the
#' results side by side:
#' * `"published"`: the simple proportional approximation used in Rasmussen
#'   et al. (2022, 2024) (eq. 16);
#' * `"simultaneous"`: the exact solution of the system of equations;
#' * `"global"`: the iterative (maximum-entropy) model of disease
#'   combinations, which can also include interactions.
#'
#' For a model, every method is run on the model's inputs; methods that
#' cannot be run (e.g. the published approximation for adjusted estimands)
#' are reported with the reason. For [cm_analyses()], each analysis is
#' compared and the tables are stacked. For a hazard-ratio model, the
#' methods of [deconflate_hr()] are compared. For a Monte Carlo run made with
#' several methods, the methods were applied to identical draws, and their
#' summaries are compared.
#'
#' @param x A [cm_model()], [cm_analyses()], [cm_hr_model()], or a `cm_mc`
#'   object from [cm_monte_carlo()] run with more than one method.
#' @param methods Methods to compare.
#' @param valuation Optional valuation list for models (see
#'   [contribution_table()]), or a named list of them (one per analysis) for
#'   [cm_analyses()] (analyses without one are not valued); adds the gap and
#'   value per method to `totals`.
#' @param overall_risk For hazard-ratio models: optional overall risk, adding
#'   the attributable risk per method (see [attributable_risk()]).
#' @param stat For Monte Carlo runs: `"mean"`, `"median"` or
#'   `"trimmed_mean"`.
#' @param ... Passed to [deconflate()] or [deconflate_hr()] (for Monte Carlo
#'   runs, to [summary.cm_mc()]).
#' @return A `cm_comparison` object with `impacts` (raw and adjusted values,
#'   one column per method), `change` (relative change from raw), `long`
#'   (long format with sign-change flags, or the Monte Carlo summaries),
#'   `totals` (aggregate per method, with gap and value if requested),
#'   `diagnostics`, `failed` (methods that could not be run, or whose
#'   valuation or attributable risk could not be computed, with reasons; their
#'   `totals` are `NA`) and `methods`.
#' @export
#' @examples
#' compare_methods(example_supplement())
compare_methods <- function(x, ...) UseMethod("compare_methods")

#' @rdname compare_methods
#' @export
compare_methods.cm_model <- function(x, methods = c("published", "simultaneous", "global"),
                                     valuation = NULL, ...) {
  methods <- unique(match.arg(methods, c("published", "simultaneous", "global"), several.ok = TRUE))
  res <- list()
  failed <- character(0)
  for (m in methods) {
    r <- tryCatch(deconflate(x, method = m, warn = FALSE, ...), deconflate_error = function(e) e)
    if (inherits(r, "condition")) failed[m] <- conditionMessage(r) else res[[m]] <- r
  }
  if (!length(res)) {
    cm_abort(sprintf("All methods failed: %s", paste(names(failed), failed, sep = ": ", collapse = "; ")))
  }
  # A valuation that cannot be evaluated for a method (e.g. an aggregate loss
  # of 100% or more) is reported, not fatal.
  vals <- list()
  if (!is.null(valuation)) {
    check_valuation(valuation)
    for (m in names(res)) {
      ev <- tryCatch(evaluate_valuation(res[[m]], valuation), deconflate_error = function(e) e)
      if (inherits(ev, "condition")) {
        failed[m] <- paste("valuation:", conditionMessage(ev))
      } else {
        vals[[m]] <- ev
      }
    }
  }
  ok <- names(res)
  base <- res[[1]]$adjusted
  impacts <- base[, c("disease", "raw")]
  change <- impacts
  for (m in ok) {
    impacts[[m]] <- res[[m]]$adjusted$adjusted
    change[[m]] <- res[[m]]$adjusted$change
  }
  long <- do.call(rbind, lapply(ok, function(m) {
    a <- res[[m]]$adjusted
    data.frame(disease = a$disease, method = m, raw = a$raw, adjusted = a$adjusted,
               change = a$change, contribution = res[[m]]$contributions$total,
               sign_change = is.finite(a$adjusted) & abs(a$raw) > 1e-12 & abs(a$adjusted) > 1e-12 &
                 sign(a$adjusted) != sign(a$raw),
               stringsAsFactors = FALSE)
  }))
  totals <- do.call(rbind, lapply(ok, function(m) {
    r <- res[[m]]
    tt <- data.frame(method = m, raw_sum = r$totals$raw_sum, adjusted_total = r$totals$adjusted_total,
                     stringsAsFactors = FALSE)
    if (!is.null(valuation)) {
      ev <- vals[[m]]
      tt$gap <- if (is.null(ev)) NA_real_ else ev$gap$summary$gap
      if (!is.null(valuation$unit_value)) tt$value <- if (is.null(ev)) NA_real_ else ev$value$value
    }
    tt
  }))
  diagnostics <- do.call(rbind, lapply(ok, function(m) res[[m]]$diagnostics))
  rownames(totals) <- NULL
  rownames(diagnostics) <- NULL
  structure(list(impacts = impacts, change = change, long = long, totals = totals,
                 diagnostics = diagnostics, failed = failed, methods = ok, source = "model",
                 units = res[[1]]$units, label = res[[1]]$label, results = res),
            class = "cm_comparison")
}

#' @rdname compare_methods
#' @export
compare_methods.cm_analyses <- function(x, methods = c("published", "simultaneous", "global"),
                                        valuation = NULL, ...) {
  nms <- names(x$models)
  # One valuation per analysis: a single valuation list would otherwise be
  # ignored silently (its elements are not named after analyses).
  if (!is.null(valuation) && (!is.list(valuation) || is.null(names(valuation)) ||
                              length(setdiff(names(valuation), nms)))) {
    cm_abort(sprintf("For cm_analyses, `valuation` must be a list named after the analyses (%s), each element a valuation list.",
                     paste(nms, collapse = ", ")))
  }
  # The global method needs one joint distribution for all analyses.
  dots <- list(...)
  if ("global" %in% methods && is.null(dots$joint)) {
    fj <- dots[intersect(names(dots), c("tol", "max_iter", "max_diseases"))]
    j <- tryCatch(withCallingHandlers(do.call(fit_joint, c(list(x$population), fj)),
                                      deconflate_nonconvergence = function(w) invokeRestart("muffleWarning")),
                  deconflate_error = function(e) NULL)
    if (!is.null(j)) dots$joint <- j
  }
  cmps <- lapply(nms, function(nm) {
    do.call(compare_methods, c(list(x$models[[nm]], methods = methods, valuation = valuation[[nm]]),
                               dots))
  })
  names(cmps) <- nms
  # Analyses can have different columns (e.g. a method that failed in one
  # analysis only, or a valuation for some analyses): fill with NA.
  stack <- function(el) {
    parts <- lapply(nms, function(nm) {
      d <- cmps[[nm]][[el]]
      if (is.null(d) || !nrow(d)) return(NULL)
      cbind(analysis = nm, d, stringsAsFactors = FALSE)
    })
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (!length(parts)) return(NULL)
    cols <- unique(unlist(lapply(parts, names)))
    parts <- lapply(parts, function(d) {
      for (cl in setdiff(cols, names(d))) d[[cl]] <- NA
      d[cols]
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    out
  }
  failed <- unlist(lapply(nms, function(nm) {
    f <- cmps[[nm]]$failed
    if (length(f)) stats::setNames(f, paste(nm, names(f), sep = ": ")) else NULL
  }))
  structure(list(impacts = stack("impacts"), change = stack("change"), long = stack("long"),
                 totals = stack("totals"), diagnostics = stack("diagnostics"),
                 failed = failed %||% character(0),
                 methods = unique(unlist(lapply(cmps, `[[`, "methods"))),
                 source = "analyses", comparisons = cmps),
            class = "cm_comparison")
}

#' @rdname compare_methods
#' @export
compare_methods.cm_hr_model <- function(x, methods = c("published", "first_order", "snapshot"),
                                        overall_risk = NULL, ...) {
  methods <- unique(match.arg(methods, c("published", "first_order", "snapshot"), several.ok = TRUE))
  if (!is.null(overall_risk)) {
    check_numeric(overall_risk, "overall_risk")
    if (length(overall_risk) != 1L || overall_risk <= 0 || overall_risk >= 1) {
      cm_abort("`overall_risk` must be a single proportion between 0 and 1.")
    }
  }
  res <- list()
  failed <- character(0)
  for (m in methods) {
    r <- tryCatch(deconflate_hr(x, method = m, warn = FALSE, ...), deconflate_error = function(e) e)
    if (inherits(r, "condition")) failed[m] <- conditionMessage(r) else res[[m]] <- r
  }
  if (!length(res)) cm_abort("All methods failed.")
  ok <- names(res)
  impacts <- res[[1]]$adjusted[, c("disease", "raw")]
  change <- impacts
  for (m in ok) {
    impacts[[m]] <- res[[m]]$adjusted$adjusted
    change[[m]] <- res[[m]]$adjusted$change
  }
  totals <- NULL
  if (!is.null(overall_risk)) {
    joint <- NULL
    totals <- do.call(rbind, lapply(ok, function(m) {
      r <- res[[m]]
      joint <<- joint %||% r$joint %||% res$snapshot$joint %||% fit_joint(x$population)
      ar <- tryCatch(attributable_risk(r, overall_risk, joint = joint, allocate = FALSE),
                     deconflate_error = function(e) e)
      if (inherits(ar, "condition")) {
        failed[m] <<- paste("attributable risk:", conditionMessage(ar))
        return(data.frame(method = m, overall_risk = overall_risk, disease_free_risk = NA_real_,
                          attributable = NA_real_, stringsAsFactors = FALSE))
      }
      data.frame(method = m, overall_risk = overall_risk,
                 disease_free_risk = ar$summary$disease_free_risk,
                 attributable = ar$summary$attributable, stringsAsFactors = FALSE)
    }))
  }
  diagnostics <- do.call(rbind, lapply(ok, function(m) res[[m]]$diagnostics))
  structure(list(impacts = impacts, change = change, long = NULL, totals = totals,
                 diagnostics = diagnostics, failed = failed, methods = ok, source = "hr",
                 units = "hazard ratio", results = res),
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
  ids <- unique(d$disease)
  raw <- vapply(ids, function(id) {
    sel <- d$disease == id & d$method == x$method[1]
    w <- x$weights[match(d$draw[sel], x$params$draw)]
    sum(w * d$raw[sel]) / sum(w)
  }, numeric(1))
  impacts <- data.frame(disease = ids, raw_mean = unname(raw), stringsAsFactors = FALSE)
  stability <- data.frame(disease = ids, stringsAsFactors = FALSE)
  for (m in x$method) {
    sm <- s[s$method == m, , drop = FALSE]
    idx <- match(ids, sm$disease)
    impacts[[m]] <- sm[[col]][idx]
    stability[[m]] <- sm$stability[idx]
  }
  tt <- summary(x, what = "total", diagnose = FALSE)
  structure(list(impacts = impacts, stability = stability, long = s,
                 totals = tt[, c("quantity", "method", "mean", "q0.5", "trimmed_mean", "mcse", "stability")],
                 failed = character(0), methods = x$method, source = "mc", stat = stat,
                 n_draws = x$n_draws, n_rejected = x$n_rejected, units = x$units),
            class = "cm_comparison")
}

#' @export
print.cm_comparison <- function(x, digits = 3, ...) {
  cat(sprintf("<cm_comparison> methods: %s\n", paste(x$methods, collapse = ", ")))
  if (!is.null(x$units)) cat(sprintf("Units: %s\n", x$units))
  if (identical(x$source, "mc")) {
    cat(sprintf("Monte Carlo: %d draws (%d rejected); statistic: %s\n",
                x$n_draws, x$n_rejected, x$stat))
  }
  tab <- x$impacts
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = digits)
  cat("\nAdjusted values:\n")
  print(tab, row.names = FALSE)
  if (!is.null(x$long) && !is.null(x$long$sign_change)) {
    sc <- x$long[x$long$sign_change, , drop = FALSE]
    if (nrow(sc)) {
      cat("\nSign changes (adjusted impact on the other side of zero from the raw impact):\n")
      for (m in unique(sc$method)) {
        s <- sc[sc$method == m, , drop = FALSE]
        lab <- if (!is.null(s$analysis)) paste(s$analysis, s$disease) else s$disease
        cat(sprintf("  %s: %s\n", m, paste(lab, collapse = ", ")))
      }
    }
  }
  if (!is.null(x$stability)) {
    flagged <- vapply(x$methods, function(m) sum(x$stability[[m]] != "ok", na.rm = TRUE), numeric(1))
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
  if (length(x$failed)) {
    cat("\nNot run (or not valued):\n")
    for (m in names(x$failed)) cat(sprintf("  %s: %s\n", m, x$failed[[m]]))
  }
  invisible(x)
}
