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
#'
#' @return A data frame with `risk_exposed`, `risk_unexposed` and `excess`
#'   (their difference).
#' @export
#' @examples
#' # Displaced abomasum, Rasmussen et al. (2022) Table 6: HR 3.83, prevalence 0.03,
#' # culling rate 0.27; published excess probability 0.31.
#' hr_to_risk(3.83, 0.03, 0.27, method = "or_approx")
#' hr_to_risk(3.83, 0.03, 0.27)
hr_to_risk <- function(hr, prevalence, overall_risk,
                       method = c("proportional_hazards", "or_approx")) {
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
    if (method == "or_approx") {
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
#' @param method `"proportional_hazards"` (default) or `"or_approx"` (the
#'   published approach of Rasmussen et al. 2022); see [hr_to_risk()].
#' @return A `cm_hr` data frame with one row per disease.
#' @export
#' @examples
#' conv <- hr_conversion(example_supplement(), c(d1 = 2, d2 = 1.5), overall_risk = 0.27)
#' conv
#' as_impacts(conv)
hr_conversion <- function(diseases, hr, overall_risk,
                          method = c("proportional_hazards", "or_approx")) {
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

#' Convert adjusted excess risks back to hazard ratios
#'
#' @param result A [deconflate()] result.
#' @param conversion The [hr_conversion()] used to build the impacts.
#' @param outcome Outcome label of the culling impacts.
#' @param method `"proportional_hazards"` (default): inverts the
#'   proportional-hazards conversion with each disease's unexposed risk
#'   ([excess_to_hr()]). `"published"`: rescales the raw hazard ratio by the
#'   ratio of adjusted to raw excess risk (Rasmussen et al. 2022, eq. 23).
#' @return A data frame with raw and adjusted excess risks and hazard ratios.
#' @export
adjusted_hr <- function(result, conversion, outcome = "culling",
                        method = c("proportional_hazards", "published")) {
  method <- match.arg(method)
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  if (!inherits(conversion, "cm_hr")) cm_abort("`conversion` must come from hr_conversion().")
  a <- result$adjusted[result$adjusted$outcome == outcome, , drop = FALSE]
  if (!nrow(a)) cm_abort(sprintf("Outcome '%s' is not in the result.", outcome))
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
