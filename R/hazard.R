#' Convert culling (or mortality) hazard ratios to period risks
#'
#' Converts a hazard ratio for disease `i` into the period risk of the event
#' (e.g. culling within a year) among animals with and without the disease,
#' given the disease's prevalence and the overall period risk.
#'
#' @param hr Hazard ratio(s).
#' @param prevalence Disease probability (same length as `hr`, or recycled).
#' @param overall_risk Overall period risk of the event in the population.
#' @param method
#'   * `"proportional_hazards"` (default): assumes proportional hazards over
#'     the period, so `risk_exposed = 1 - (1 - risk_unexposed)^hr`, with
#'     `risk_unexposed` solved so that the population risk equals
#'     `overall_risk`.
#'   * `"or_approx"`: treats the hazard ratio as an odds ratio in the 2x2
#'     table of disease by event (as in Rasmussen et al. 2022, section
#'     2.4.4 and Table 6). Provided to reproduce the published values.
#'   * `"overall_odds"`: treats the hazard ratio as an odds ratio relative to
#'     the overall risk, `risk_exposed = hr * r / (hr * r + 1 - r)`, with the
#'     overall risk `r` as the reference (`risk_unexposed = r`), as in the
#'     loss calculations of Rasmussen et al. (2024). `prevalence` is not used.
#'
#' @return A data frame with `risk_exposed`, `risk_unexposed` (the reference
#'   risk) and `excess` (their difference).
#' @export
#' @examples
#' # Displaced abomasum, Rasmussen et al. (2022) Table 6: HR 3.83, prevalence 0.03,
#' # culling rate 0.27; published excess probability 0.31.
#' hr_to_risk(3.83, 0.03, 0.27, method = "or_approx")
#' hr_to_risk(3.83, 0.03, 0.27)
#' # Rasmussen et al. (2024): adjusted HR relative to the overall culling risk
#' hr_to_risk(2.75, NA, 0.27, method = "overall_odds")
hr_to_risk <- function(hr, prevalence, overall_risk,
                       method = c("proportional_hazards", "or_approx", "overall_odds")) {
  method <- match.arg(method)
  n <- max(length(hr), length(prevalence), length(overall_risk))
  hr <- recycle_arg(hr, n, "hr")
  prevalence <- recycle_arg(prevalence, n, "prevalence")
  overall_risk <- recycle_arg(overall_risk, n, "overall_risk")
  check_numeric(hr, "hr")
  if (any(hr <= 0)) cm_abort("Hazard ratios must be positive.")
  res <- t(vapply(seq_len(n), function(j) {
    h <- hr[j]
    P <- prevalence[j]
    r <- overall_risk[j]
    if (method == "overall_odds") {
      c(h * r / (h * r + 1 - r), r)
    } else if (method == "or_approx") {
      p11 <- or_to_joint(h, r, P)
      c(p11 / P, (r - p11) / (1 - P))
    } else {
      f <- function(p0) (1 - P) * p0 + P * (1 - (1 - p0)^h) - r
      p0 <- stats::uniroot(f, c(0, 1), tol = 1e-14)$root
      c(1 - (1 - p0)^h, p0)
    }
  }, numeric(2)))
  # A hazard ratio of exactly 1 means no excess risk (avoid rounding noise).
  res[hr == 1, 1] <- overall_risk[hr == 1]
  res[hr == 1, 2] <- overall_risk[hr == 1]
  data.frame(hr = hr, prevalence = prevalence, overall_risk = overall_risk,
             risk_exposed = res[, 1], risk_unexposed = res[, 2],
             excess = res[, 1] - res[, 2])
}

#' Convert an adjusted excess risk back to a hazard ratio
#'
#' Inverse of the proportional-hazards conversion in [hr_to_risk()], for
#' reporting adjusted culling impacts as hazard ratios.
#'
#' @param excess Excess period risk among diseased animals.
#' @param risk_unexposed Period risk among animals without the disease.
#' @return Hazard ratio(s): `log(1 - (risk_unexposed + excess)) / log(1 - risk_unexposed)`.
#' @export
excess_to_hr <- function(excess, risk_unexposed) {
  log(1 - (risk_unexposed + excess)) / log(1 - risk_unexposed)
}

#' Convert culling hazard ratios into excess-risk impacts
#'
#' Converts disease-specific culling (or mortality) hazard ratios into the
#' excess period risk among diseased animals, using [hr_to_risk()] with each
#' disease's probability from the model. Diseases without a hazard ratio get
#' `hr = 1` (no excess). Use [as_impacts()] to turn the result into impacts
#' for [deconflate()], and [adjusted_hr()] to convert adjusted impacts back
#' to hazard ratios.
#'
#' @param diseases A [cm_diseases()] or [cm_model()] object.
#' @param hr Named numeric vector of hazard ratios (names = disease ids).
#' @param overall_risk Overall period risk of the event in the population
#'   (e.g. the annual culling rate as a proportion).
#' @param method `"proportional_hazards"` (default), `"or_approx"` (the
#'   published approach of Rasmussen et al. 2022) or `"overall_odds"`; see
#'   [hr_to_risk()].
#' @return A `cm_hr` data frame with one row per disease.
#' @export
#' @examples
#' conv <- hr_conversion(example_supplement(), c(d1 = 2, d2 = 1.5), overall_risk = 0.27)
#' conv
#' as_impacts(conv)
hr_conversion <- function(diseases, hr, overall_risk,
                          method = c("proportional_hazards", "or_approx", "overall_odds")) {
  method <- match.arg(method)
  if (inherits(diseases, "cm_model")) diseases <- diseases$diseases
  if (!inherits(diseases, "cm_diseases")) cm_abort("`diseases` must be a cm_diseases or cm_model object.")
  if (is.null(names(hr))) cm_abort("`hr` must be named by disease id.")
  unk <- setdiff(names(hr), diseases$id)
  if (length(unk)) cm_abort(sprintf("Unknown diseases in `hr`: %s.", paste(unk, collapse = ", ")))
  h <- stats::setNames(rep(1, nrow(diseases)), diseases$id)
  h[names(hr)] <- hr
  r <- hr_to_risk(unname(h), diseases$prob, overall_risk, method = method)
  out <- data.frame(disease = diseases$id, r, method = method, stringsAsFactors = FALSE)
  class(out) <- c("cm_hr", "data.frame")
  out
}

#' Turn a hazard-ratio conversion into impacts
#'
#' @param x A [hr_conversion()] result.
#' @param outcome Outcome label.
#' @param scale How the excess risk enters the productivity gap:
#'   * `"absolute"` (default): excess risk in probability units, so the
#'     disease-free risk is `observed - sum(excess * P)`. Use with an
#'     observed risk given as a proportion (e.g. `0.27`).
#'   * `"proportion"`: excess risk treated as a proportional increase in the
#'     observed rate, so the disease-free rate is `observed / (1 + L)`, as in
#'     Rasmussen et al. (2022). Use this to reproduce the published results.
#' @param ... Passed to [cm_impacts()] (e.g. `source`).
#' @return A [cm_impacts()] object with direction `"increase"`.
#' @export
as_impacts <- function(x, outcome = "culling", scale = c("absolute", "proportion"), ...) {
  if (!inherits(x, "cm_hr")) cm_abort("`x` must come from hr_conversion().")
  scale <- match.arg(scale)
  cm_impacts(x$disease, x$excess, outcome = outcome, scale = scale,
             units = if (scale == "absolute") "probability" else NA_character_,
             direction = "increase", ...)
}

#' Convert adjusted culling impacts back to hazard ratios
#'
#' @param result A [deconflate()] result.
#' @param conversion The [hr_conversion()] used to build the impacts. Not
#'   needed for `method = "excess_hr"`.
#' @param outcome Outcome label of the culling impacts.
#' @param method
#'   * `"proportional_hazards"` (default): inverts the proportional-hazards
#'     conversion with each disease's unexposed risk ([excess_to_hr()]).
#'   * `"published"`: rescales the raw hazard ratio by the ratio of adjusted
#'     to raw excess risk (Rasmussen et al. 2022, eq. 23).
#'   * `"excess_hr"`: the impacts are hazard ratios minus 1, adjusted
#'     directly (Rasmussen et al. 2024; see [example_global_dairy()]), so the
#'     adjusted hazard ratio is the adjusted impact plus 1.
#'
#'   Outcomes on the `"hazard_ratio"` scale are already hazard ratios; they
#'   are returned as they are, whatever the method.
#' @return A data frame with raw and adjusted impacts (`excess`,
#'   `excess_adjusted`: excess risks, or HR - 1 for `"excess_hr"`) and hazard
#'   ratios (`hr`, `hr_adjusted`).
#' @export
adjusted_hr <- function(result, conversion = NULL, outcome = "culling",
                        method = c("proportional_hazards", "published", "excess_hr")) {
  method <- match.arg(method)
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  a <- result$adjusted[result$adjusted$outcome == outcome, , drop = FALSE]
  if (!nrow(a)) cm_abort(sprintf("Outcome '%s' is not in the result.", outcome))
  if (a$scale[1] == "hazard_ratio" || method == "excess_hr") {
    if (a$scale[1] == "hazard_ratio") {
      return(data.frame(disease = a$disease, hr = a$raw, excess = a$raw - 1,
                        excess_adjusted = a$adjusted - 1, hr_adjusted = a$adjusted,
                        stringsAsFactors = FALSE))
    }
    return(data.frame(disease = a$disease, hr = a$raw + 1, excess = a$raw,
                      excess_adjusted = a$adjusted, hr_adjusted = a$adjusted + 1,
                      stringsAsFactors = FALSE))
  }
  if (!inherits(conversion, "cm_hr")) cm_abort("`conversion` must come from hr_conversion().")
  if (identical(conversion$method[1], "overall_odds")) {
    cm_abort("Conversions with method 'overall_odds' cannot be inverted; use them for losses only.")
  }
  cv <- conversion[match(a$disease, conversion$disease), , drop = FALSE]
  hr_adj <- if (method == "published") {
    ifelse(cv$excess != 0, a$adjusted * cv$hr / cv$excess, cv$hr)
  } else {
    excess_to_hr(a$adjusted, cv$risk_unexposed)
  }
  data.frame(disease = a$disease, hr = cv$hr, excess = cv$excess,
             excess_adjusted = a$adjusted, hr_adjusted = hr_adj,
             stringsAsFactors = FALSE)
}

#' Culling (or mortality) attributable to disease
#'
#' Converts adjusted hazard ratios into the part of an event's overall risk
#' (e.g. annual culling) that is attributable to disease, without counting
#' an animal with several diseases more than once, and allocates it to
#' diseases.
#'
#' Within the period, an animal with disease combination `x` has a constant
#' hazard `h0 * exp(sum_i beta[i] * x[i])`, with `beta = log(adjusted HR)`.
#' The baseline hazard `h0` is chosen so that the population risk, averaged
#' over the distribution of disease combinations ([fit_joint()]), equals
#' `overall_risk`:
#'
#' `sum_x P(x) * (1 - exp(-h0 * exp(beta . x))) = overall_risk`.
#'
#' The disease-free risk is `1 - exp(-h0)`, and the attributable risk is
#' `overall_risk - (1 - exp(-h0))`. An animal's risk cannot exceed 1, so the
#' attributable risk is smaller than the sum of per-disease excess risks
#' when diseases co-occur.
#'
#' The attributable risk is allocated to diseases by Shapley values over
#' disease combinations ([shapley_by_cell()]), with the loss
#' `1 - exp(-h0 * exp(beta . x)) - (1 - exp(-h0))`.
#'
#' The model is consistent with `method = "global"` in [deconflate()]: there,
#' the adjusted hazard ratios reproduce the raw ones exactly over the same
#' joint distribution. Hazard ratios from the other methods are used as they
#' are.
#'
#' @param result A [deconflate()] result with a hazard-ratio outcome (see
#'   [cm_impacts()]).
#' @param overall_risk Overall period risk of the event, as a proportion
#'   (e.g. `0.25` for an annual culling rate of 25%).
#' @param outcome Outcome label.
#' @param unit_value Optional value per animal removed (e.g. replacement
#'   cost less salvage value); adds `value` columns.
#' @param joint Optional [fit_joint()] result; by default the result's own
#'   joint distribution (global method) or a new fit.
#' @param allocate Logical: allocate the attributable risk to diseases?
#'   This is the slow step for many co-occurring diseases.
#' @param max_present Passed to [shapley_by_cell()].
#' @return A `cm_attributable` list with `summary` (overall, disease-free and
#'   attributable risk, attributable fraction and value), `by_disease`
#'   (adjusted hazard ratio, attributable risk, share and value),
#'   `skipped_mass` and `baseline_hazard`.
#' @export
#' @examples
#' m <- example_supplement()
#' m$impacts <- combine_impacts(m$impacts,
#'   cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), outcome = "culling",
#'              scale = "hazard_ratio"))
#' res <- deconflate(m, method = "global")
#' attributable_risk(res, overall_risk = 0.25, unit_value = 1300)
attributable_risk <- function(result, overall_risk, outcome = "culling", unit_value = NULL,
                              joint = NULL, allocate = TRUE, max_present = 10L) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  a <- result$adjusted[result$adjusted$outcome == outcome, , drop = FALSE]
  if (!nrow(a)) cm_abort(sprintf("Outcome '%s' is not in the result.", outcome))
  if (a$scale[1] != "hazard_ratio") {
    cm_abort(sprintf("Outcome '%s' is not on the hazard-ratio scale.", outcome))
  }
  check_numeric(overall_risk, "overall_risk")
  if (length(overall_risk) != 1L || overall_risk <= 0 || overall_risk >= 1) {
    cm_abort("`overall_risk` must be a single proportion between 0 and 1.")
  }
  bad <- !is.finite(a$adjusted) | a$adjusted <= 0
  if (any(bad)) {
    cm_abort(sprintf("Adjusted hazard ratios must be positive and finite (check: %s).",
                     paste(a$disease[bad], collapse = ", ")))
  }
  joint <- joint %||% result$joint %||% fit_joint(result$model)
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().")
  ids <- joint$diseases
  beta <- stats::setNames(log(a$adjusted), a$disease)[ids]
  rel <- exp(as.vector(joint$cells %*% beta))
  risk_at <- function(h) sum(joint$prob * (1 - exp(-h * rel)))
  h0 <- stats::uniroot(function(h) risk_at(h) - overall_risk,
                       c(0, -log(1 - overall_risk)), extendInt = "upX", tol = 1e-14)$root
  r0 <- 1 - exp(-h0)
  summ <- data.frame(outcome = outcome, overall_risk = overall_risk, disease_free_risk = r0,
                     attributable = overall_risk - r0,
                     attributable_fraction = (overall_risk - r0) / overall_risk,
                     stringsAsFactors = FALSE)
  if (!is.null(unit_value)) summ$value <- summ$attributable * unit_value
  by <- NULL
  skipped <- 0
  if (allocate) {
    loss <- function(x) 1 - exp(-h0 * exp(sum(beta * x[ids]))) - r0
    sh <- shapley_by_cell(joint, loss, max_present = max_present)
    by <- data.frame(disease = ids, hr_adjusted = unname(exp(beta)),
                     attributable = sh$shapley, share = sh$share, stringsAsFactors = FALSE)
    if (!is.null(unit_value)) by$value <- by$attributable * unit_value
    skipped <- attr(sh, "skipped_mass")
  }
  structure(list(summary = summ, by_disease = by, skipped_mass = skipped,
                 baseline_hazard = h0),
            class = "cm_attributable")
}

#' @export
print.cm_attributable <- function(x, ...) {
  s <- x$summary
  cat(sprintf("<cm_attributable> outcome: %s\n", s$outcome))
  cat(sprintf("  Overall risk %.4g; disease-free risk %.4g; attributable %.4g (%.1f%% of the overall risk)\n",
              s$overall_risk, s$disease_free_risk, s$attributable, 100 * s$attributable_fraction))
  if (!is.null(s$value)) cat(sprintf("  Value: %.4g\n", s$value))
  if (!is.null(x$by_disease)) {
    cat("\n")
    print(x$by_disease, row.names = FALSE, digits = 4)
    if (x$skipped_mass > 0) {
      cat(sprintf("\nCombinations skipped by `max_present`: probability %.2e\n", x$skipped_mass))
    }
  }
  invisible(x)
}
