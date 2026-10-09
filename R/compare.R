#' Compare adjustment methods
#'
#' Runs several adjustment methods on the same inputs and tabulates the
#' results side by side. For additive impacts:
#' * `"simultaneous"`: the exact solution from the pairwise tables;
#' * `"global"`: the exact solution from the maximum-entropy distribution of
#'   disease combinations, which also handles interactions, three-way terms
#'   and unknown pairs;
#' * `"published"`: the proportional approximation used in Rasmussen et al.
#'   (2022, 2024) (eq. 16 of the 2022 paper). It is not a method of
#'   [deconflate()]; it is kept here (and in the `reproduce_*()` functions)
#'   to compare earlier results with the exact solution. It applies to crude
#'   estimates only, can mask incompatible inputs, and is undefined when its
#'   denominator is zero.
#'
#' For event impacts (`event_model = TRUE`): `"snapshot"` (the model of
#' [deconflate()]), `"first_order"` (its log-linear approximation, from the
#' pairwise tables) and `"published"` (the approach of Rasmussen et al.
#' 2024: `HR - 1` adjusted with eq. 16), each with the attributable risk at
#' `overall_risk`. The first-order and published methods need hazard (or
#' rate) ratios.
#'
#' Methods run exactly as asked (there is no automatic switch): methods that
#' cannot be run (e.g. a pairwise method with unknown pairs, or the published
#' approximation for adjusted estimands) are reported with the reason, and
#' results with non-finite values are kept apart in `undefined`.
#'
#' With `n_draws`, the methods are also applied to the same draws of the
#' uncertain inputs (a draw is rejected if any method fails on it), and the
#' means and 95% intervals are reported in `draws`.
#'
#' @param model A [cm_model()].
#' @param methods Methods to compare (default: all for the kind of impacts).
#' @param event_model,overall_risk As in [deconflate()].
#' @param n_draws Number of draws (default 0: point estimates only).
#' @param seed Optional random seed for the draws.
#' @param ... Passed to the adjustment (e.g. `joint`, `feasibility`, or
#'   arguments of [fit_joint()]).
#' @return A `cm_comparison` object with `impacts` (raw and adjusted values,
#'   one column per method), `change` (relative change from raw), `long`
#'   (long format with sign-change flags), `totals` (naive and adjusted
#'   aggregate per method; for event impacts, the attributable risk),
#'   `diagnostics`, `failed` (methods that could not be run, or gave an
#'   undefined result, with reasons), `undefined`, `methods`, `results` and,
#'   with draws, `draws` (summary of each quantity by method, rejections).
#' @export
#' @examples
#' compare_methods(example_supplement())
compare_methods <- function(model, methods = NULL, event_model = FALSE, overall_risk = NULL,
                            n_draws = 0, seed = NULL, ...) {
  check_model(model)
  check_event_model(model, event_model)
  risk <- check_overall_risk(overall_risk, event_model)
  all_m <- if (event_model) c("published", "first_order", "snapshot") else c("published", "simultaneous", "global")
  methods <- unique(match.arg(methods %||% all_m, all_m, several.ok = TRUE))
  check_draw_args(n_draws, seed)
  # As in deconflate(): beyond 20 diseases the joint distribution is fitted
  # with the sampled backend (unless a backend or a joint is given).
  auto_args <- if (nrow(model$diseases) > 20L && is.null(list(...)$backend) && is.null(list(...)$joint)) {
    list(backend = "sampled", seed = seed)
  } else list()
  run_m <- function(mod, m, r, ...) {
    args <- utils::modifyList(auto_args, list(...))
    if (event_model) do.call(adjust_event, c(list(mod, method = m, overall_risk = r, warn = FALSE), args)) else
      do.call(adjust_impacts, c(list(mod, method = m, warn = FALSE), args))
  }
  res <- list()
  undefined <- list()
  failed <- character(0)
  for (m in methods) {
    r <- tryCatch(withCallingHandlers(run_m(model, m, risk$value, ...),
                                      deconflate_nonconvergence = function(w) {
                                        if (inherits(w, "warning")) invokeRestart("muffleWarning")
                                      }),
                  deconflate_error = function(e) e)
    if (inherits(r, "condition")) {
      failed[m] <- conditionMessage(r)
    } else if (!result_is_finite(r)) {
      # Kept for inspection in `undefined`, but not presented as an estimate.
      failed[m] <- sprintf("undefined: non-finite adjusted values (%s)", nonfinite_diseases(r))
      undefined[[m]] <- r
    } else {
      res[[m]] <- r
    }
  }
  if (!length(res)) {
    cm_abort(sprintf("All methods failed: %s", paste(names(failed), failed, sep = ": ", collapse = "; ")))
  }
  ok <- names(res)
  base <- res[[1]]$adjusted
  impacts <- base[, c("disease", "raw")]
  if (event_model) impacts <- base[, c("disease", "measure", "raw")]
  change <- impacts
  for (m in ok) {
    impacts[[m]] <- res[[m]]$adjusted$adjusted
    change[[m]] <- res[[m]]$adjusted$change
  }
  long <- do.call(rbind, lapply(ok, function(m) {
    a <- res[[m]]$adjusted
    contrib <- if (event_model) res[[m]]$attributable$by_disease$attributable else res[[m]]$contributions$total
    raw_dir <- if (event_model) sign(event_scale(a$raw, a$measure)) else sign(a$raw)
    adj_dir <- if (event_model) sign(log(pmax(a$adjusted, 1e-300))) else sign(a$adjusted)
    data.frame(disease = a$disease, method = m, raw = a$raw, adjusted = a$adjusted,
               change = a$change, contribution = contrib %||% NA_real_,
               sign_change = is.finite(a$adjusted) & raw_dir != 0 & adj_dir != 0 & adj_dir != raw_dir,
               stringsAsFactors = FALSE)
  }))
  totals <- do.call(rbind, lapply(ok, function(m) {
    r <- res[[m]]
    if (event_model) {
      s <- r$attributable$summary
      data.frame(method = m, overall_risk = s$overall_risk, disease_free_risk = s$disease_free_risk,
                 attributable = s$attributable, stringsAsFactors = FALSE)
    } else {
      data.frame(method = m, raw_sum = r$totals$raw_sum, adjusted_total = r$totals$adjusted_total,
                 stringsAsFactors = FALSE)
    }
  }))
  diagnostics <- do.call(rbind, lapply(ok, function(m) res[[m]]$diagnostics))
  rownames(totals) <- NULL
  rownames(diagnostics) <- NULL
  out <- structure(list(impacts = impacts, change = change, long = long, totals = totals,
                        diagnostics = diagnostics, failed = failed, methods = ok,
                        event_model = event_model, units = res[[1]]$units, label = res[[1]]$label,
                        results = res, undefined = undefined, draws = NULL),
                   class = "cm_comparison")
  specs <- draw_specs(model, risk$dist)
  if (n_draws > 0 && length(specs)) {
    dots <- list(...)
    dots$joint <- NULL   # each draw has its own joint distribution
    evaluate <- function(m2, r) {
      q <- lapply(ok, function(m) {
        rr <- do.call(run_m, c(list(m2, m, r %||% risk$value), dots))
        if (!result_is_finite(rr)) {
          cm_abort(sprintf("Non-finite adjusted values (method %s).", m), class = "deconflate_nonfinite")
        }
        result_quantities(rr, prefix = paste0(m, "|"))
      })
      unlist(q)
    }
    dr <- run_draws(model, specs, as.integer(n_draws), seed, "random", 10L, evaluate)
    dr$summary <- summarise_draws(dr)
    if (!is.null(dr$summary)) {
      parts <- strsplit(dr$summary$quantity, "|", fixed = TRUE)
      dr$summary$method <- vapply(parts, `[`, character(1), 1)
      dr$summary$quantity <- vapply(parts, `[`, character(1), 2)
      dr$summary <- dr$summary[, c("method", setdiff(names(dr$summary), "method"))]
    }
    out$draws <- dr
    out$notes <- c(out$notes, draw_notes(dr))
  } else if (n_draws > 0) {
    out$notes <- "No input has a distribution, so no draws were run."
  }
  out
}

#' @export
print.cm_comparison <- function(x, digits = 3, ...) {
  cat(sprintf("<cm_comparison> methods: %s\n", paste(x$methods, collapse = ", ")))
  if (!is.null(x$label)) cat(sprintf("Impacts: %s\n", x$label))
  if (!is.null(x$units) && !is.na(x$units)) cat(sprintf("Units: %s\n", x$units))
  tab <- x$impacts
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], signif, digits = digits)
  cat(if (isTRUE(x$event_model)) "\nAdjusted hazard ratios:\n" else "\nAdjusted values:\n")
  print(tab, row.names = FALSE)
  sc <- x$long[x$long$sign_change, , drop = FALSE]
  if (nrow(sc)) {
    cat("\nSign changes (adjusted value on the other side of no effect from the raw estimate):\n")
    for (m in unique(sc$method)) {
      cat(sprintf("  %s: %s\n", m, paste(sc$disease[sc$method == m], collapse = ", ")))
    }
  }
  if (!is.null(x$totals) && nrow(x$totals)) {
    cat("\nTotals:\n")
    tt <- x$totals
    nn <- vapply(tt, is.numeric, logical(1))
    tt[nn] <- lapply(tt[nn], signif, digits = 4)
    print(tt, row.names = FALSE)
  }
  dr <- x$draws
  if (!is.null(dr) && !is.null(dr$summary)) {
    cat(sprintf("\nUncertainty (%d draws, %d rejected): mean and 95%% interval of the total\n",
                dr$n_draws, dr$n_rejected))
    s <- dr$summary[dr$summary$quantity == "total", c("method", "mean", "lower", "upper", "stability")]
    nn <- vapply(s, is.numeric, logical(1))
    s[nn] <- lapply(s[nn], signif, digits = 4)
    print(s, row.names = FALSE)
  }
  if (length(x$failed)) {
    cat("\nNot run:\n")
    for (m in names(x$failed)) cat(sprintf("  %s: %s\n", m, x$failed[[m]]))
  }
  if (length(x$notes)) cat(sprintf("\nNote: %s\n", x$notes), sep = "")
  invisible(x)
}
