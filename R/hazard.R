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
