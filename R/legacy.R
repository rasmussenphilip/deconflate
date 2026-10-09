# Legacy conversions used only to reproduce published results ----------------
#
# They are not supported for new analyses (see deconflate_hr() and
# attributable_risk()) and are not exported. The gap and valuation helpers
# below reproduce the economic tables of 2022; the package itself does not
# value impacts. reproduce_rasmussen_2022() uses
# the hazard ratio as an odds ratio (uk_dairy_2022_analyses(culling = TRUE))
# and eq. 23; reproduce_rasmussen_2024() adjusts HR - 1 additively
# (global_dairy_analyses(inputs, culling = TRUE)). The overall-odds excess
# risk of the 2024 losses is kept for reference and tested, but not used by
# the reproduction functions.

# Rasmussen et al. (2022), section 2.4.4: the hazard ratio is used as the odds
# ratio of a 2x2 table of disease (prevalence P) by culling (overall risk r).
legacy_hr_as_or_excess <- function(hr, prevalence, overall_risk) {
  vapply(seq_along(hr), function(j) {
    if (abs(hr[j] - 1) < 1e-12) return(0)
    p11 <- or_to_joint(hr[j], overall_risk, prevalence[j])
    p11 / prevalence[j] - (overall_risk - p11) / (1 - prevalence[j])
  }, numeric(1))
}

# Rasmussen et al. (2022), eq. 23: adjusted hazard ratio by rescaling.
legacy_eq23 <- function(hr, excess, excess_adjusted) {
  ifelse(excess != 0, excess_adjusted * hr / excess, hr)
}

# Rasmussen et al. (2024) losses: excess culling risk relative to the overall
# rate, treating the hazard ratio as an odds ratio.
legacy_overall_odds_excess <- function(hr, overall_risk) {
  hr * overall_risk / (hr * overall_risk + 1 - overall_risk) - overall_risk
}

# Productivity gap and its value (Rasmussen et al. 2022, eqs. 17-22, Tables
# 9-10), used only by reproduce_rasmussen_2022(). With aggregate L (the sum
# of contributions, divided by 100 for percent impacts) and observed mean x:
# for proportional effects the disease-free value is x / (1 - L) (decrease)
# or x / (1 + L) (increase); for absolute effects x + L or x - L. Each
# disease's gap is its contribution times the same factor.
legacy_gap <- function(result, observed, direction = c("decrease", "increase"),
                       effect = c("proportion", "percent", "absolute")) {
  direction <- match.arg(direction)
  effect <- match.arg(effect)
  ct <- result$contributions
  scale <- if (effect == "percent") 100 else 1
  L <- sum(ct$total) / scale
  x <- observed
  if (effect == "absolute") {
    factor <- 1
    xh <- if (direction == "decrease") x + L else x - L
  } else if (direction == "decrease") {
    if (!(L < 1)) cm_abort(sprintf("The aggregate proportional loss is %g; it must be below 1.", L))
    factor <- x / (1 - L)
    xh <- x / (1 - L)
  } else {
    if (!(L > -1)) cm_abort(sprintf("The aggregate proportional increase is %g; it must be above -1.", L))
    factor <- x / (1 + L)
    xh <- x / (1 + L)
  }
  gap <- if (direction == "decrease") xh - x else x - xh
  list(
    summary = data.frame(observed = x, disease_free = xh, gap = gap, aggregate = L,
                         direction = direction, effect = effect, stringsAsFactors = FALSE),
    attribution = data.frame(disease = ct$disease, gap = factor * ct$total / scale,
                             gap_main = factor * ct$main / scale,
                             gap_interaction = factor * ct$interaction / scale,
                             stringsAsFactors = FALSE)
  )
}

legacy_value <- function(gap, unit_value) {
  att <- gap$attribution
  list(by_disease = data.frame(disease = att$disease, gap = att$gap, value = att$gap * unit_value,
                               stringsAsFactors = FALSE),
       value = gap$summary$gap * unit_value)
}
