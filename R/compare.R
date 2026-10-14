#' Compare the model with the published approximation
#'
#' Runs the package's model and the earlier approximations on the same inputs
#' and tabulates the results side by side. For additive impacts:
#' * `"deconflate"`: the model of [deconflate()];
#' * `"published"`: the proportional approximation used in Rasmussen et al.
#'   (2022, 2024) (eq. 16 of the 2022 paper). It is not part of the model;
#'   it is kept here (and in the `reproduce_*()` functions) to compare
#'   earlier results with the exact solution. It needs one crude row per
#'   disease and a crude pair row for every pair, can mask incompatible
#'   inputs, and is undefined when its denominator is zero.
#'
#' For event impacts (`event_model = TRUE`): `"deconflate"` (the snapshot
#' model of [deconflate()]), `"first_order"` (its log-linear approximation,
#' from the pairwise tables) and `"published"` (the approach of Rasmussen et
#' al. 2024: `HR - 1` adjusted with eq. 16), each with the attributable risk
#' at `overall_risk`. The first-order and published methods need one hazard
#' (or rate) ratio per disease and every pair.
#'
#' Methods that cannot be run on the inputs are reported with the reason,
#' and results with non-finite values are kept apart in `undefined`.
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
#' @param ... Passed to [deconflate()]'s model (e.g. `joint`, `pool`, or
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
  retired_args(list(...))
  check_event_model(model, event_model)
  risk <- check_overall_risk(overall_risk, event_model)
  all_m <- if (event_model) c("published", "first_order", "deconflate") else c("published", "deconflate")
  methods <- unique(cm_match_arg(methods %||% all_m, all_m, several_ok = TRUE, name = "methods"))
  check_draw_args(n_draws, seed)
  dots <- list(...)
  pool <- dots$pool %||% "random"
  dots$pool <- NULL
  plan <- plan_model(model, event_model, dots, has_joint = !is.null(dots$joint), pool = pool)
  if (!is.null(seed) && identical(plan$fit_args$backend, "sampled")) plan$fit_args$seed <- seed
  run_m <- function(mod, m, r, central = TRUE, joint = dots$joint, frozen_a = NULL) {
    switch(m,
      deconflate = run_point(mod, plan, risk = r, joint = joint, warn = FALSE, dots = dots, central = central),
      published = if (event_model) adjust_event_pairwise(mod, "published", r, pool = pool, frozen = frozen_a) else
        adjust_published(mod, pool = pool, frozen = frozen_a),
      first_order = adjust_event_pairwise(mod, "first_order", r, pool = pool, frozen = frozen_a))
  }
  res <- list()
  undefined <- list()
  failed <- failed_type <- character(0)
  for (m in methods) {
    r <- tryCatch(withCallingHandlers(run_m(model, m, risk$value),
                                      deconflate_nonconvergence = function(w) {
                                        if (inherits(w, "warning")) invokeRestart("muffleWarning")
                                      }),
                  deconflate_error = function(e) e)
    if (inherits(r, "condition")) {
      failed[m] <- conditionMessage(r)
      failed_type[m] <- condition_type(r)
    } else if (!result_is_finite(r)) {
      # Kept for inspection in `undefined`, but not presented as an estimate.
      failed[m] <- sprintf("undefined: non-finite adjusted values (%s)", nonfinite_diseases(r))
      failed_type[m] <- "nonfinite"
      undefined[[m]] <- r
    } else {
      res[[m]] <- r
    }
  }
  if (!length(res)) {
    # (Of the kind of the first failure, usually the model's.)
    cm_abort(sprintf("All methods failed: %s", paste(names(failed), failed, sep = ": ", collapse = "; ")),
             class = type_class(failed_type[[1]]))
  }
  ok <- names(res)
  # One row per term of the model's result (diseases first).
  base <- (res$deconflate %||% res[[1]])$adjusted
  impacts <- base[, c("disease", "raw")]
  if (event_model) impacts <- base[, c("disease", "measure", "raw")]
  change <- impacts
  for (m in ok) {
    a <- res[[m]]$adjusted
    k <- match(impacts$disease, a$disease)
    impacts[[m]] <- a$adjusted[k]
    change[[m]] <- a$change[k]
  }
  long <- do.call(rbind, lapply(ok, function(m) {
    a <- res[[m]]$adjusted
    by <- if (event_model) res[[m]]$attributable$by_disease else res[[m]]$contributions
    contrib <- if (is.null(by)) NA_real_ else by[[if (event_model) "attributable" else "total"]][match(a$disease, by$disease)]
    raw_dir <- if (event_model) sign(event_scale(a$raw, ifelse(is.na(a$measure), "HR", a$measure))) else sign(a$raw)
    adj_dir <- if (event_model) sign(log(pmax(a$adjusted, 1e-300))) else sign(a$adjusted)
    data.frame(disease = a$disease, method = m, raw = a$raw, adjusted = a$adjusted,
               change = a$change, contribution = contrib,
               sign_change = is.finite(a$adjusted) & !is.na(raw_dir) & raw_dir != 0 & adj_dir != 0 & adj_dir != raw_dir,
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
  diagnostics <- do.call(rbind, lapply(ok, function(m) {
    d <- res[[m]]$diagnostics
    d$rank <- d$rank %||% NA_real_
    d[, c("method", "max_reconstruction_residual", "condition_number", "n_sign_changes", "sign_changes", "feasibility")]
  }))
  diagnostics$method <- ok
  rownames(totals) <- NULL
  rownames(diagnostics) <- NULL
  out <- structure(list(impacts = impacts, change = change, long = long, totals = totals,
                        diagnostics = diagnostics, failed = failed, methods = ok,
                        event_model = event_model, units = res[[1]]$units, label = res[[1]]$label,
                        results = res, undefined = undefined, draws = NULL),
                   class = "cm_comparison")
  frozen_pool <- res$deconflate$pooling$frozen
  specs <- draw_specs(model, risk$dist, frozen_pool)
  if (n_draws > 0 && length(specs)) {
    central <- res$deconflate
    frozen <- if (!is.null(central)) list(associations = central$pooling$frozen, impacts = central$row_se)
    evaluate <- function(m2, r) {
      q <- lapply(ok, function(m) {
        rr <- if (m == "deconflate") {
          run_point(m2, plan, risk = r %||% risk$value, joint = NULL, warn = FALSE, dots = dots,
                    joint_start = if (identical(central$joint$backend, "exact")) central$joint, frozen = frozen)
        } else run_m(m2, m, r %||% risk$value, central = FALSE, joint = NULL, frozen_a = frozen$associations)
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
