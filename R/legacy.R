# Legacy conversions used only to reproduce published results ----------------
#
# They are not supported for new analyses (see deconflate_hr() and
# attributable_risk()) and are not exported. reproduce_rasmussen_2022() uses
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
