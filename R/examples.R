#' Worked example from the Supplementary File of Rasmussen et al. (2022)
#'
#' Three hypothetical diseases with prevalences 0.10, 0.15 and 0.20, yield
#' impacts of 2.5%, 5% and 7.5%, and odds ratios of 2 (d1:d2), 1 (d1:d3) and
#' 3 (d2:d3). The observed mean yield in the example is 10,000 units per
#' animal per year.
#'
#' Note: the published productivity gap (10,225; 20, 61 and 143 units) rounds
#' the adjusted impacts to 0.02, 0.04 and 0.07 before computing the gap.
#' Without rounding, the published method gives 10,215.7 and 21.1, 56.7 and
#' 137.9 units.
#'
#' @return A [cm_model()].
#' @export
example_supplement <- function() {
  cm_model(
    diseases = cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20),
                           type = "prevalence"),
    associations = cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"),
                                   c(2, 1, 3), measure = "OR"),
    impacts = cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5),
                         outcome = "yield", scale = "percent",
                         units = "units/animal/year",
                         source = "Rasmussen et al. 2022, Table S1")
  )
}

#' Global dairy inputs from Rasmussen et al. (2024), at their means
#'
#' Global (herd-size weighted) lactational incidence (prevalence for PTB) from
#' Table 2, pooled inter-disease odds ratios from Table 3, and raw yield
#' (% decrease) and fertility (% increase in calving interval) impacts from
#' Table 4. Pairs absent from Table 3 are treated as independent, as in the
#' paper. Culling impacts (hazard ratios) are omitted; see [hr_to_risk()].
#'
#' The paper's Table 5 reports means of adjusted impacts over Monte Carlo
#' draws, which differ from adjusting the input means (e.g. displaced
#' abomasum yield: 1.18 vs 0.79).
#'
#' @return A [cm_model()] with outcomes `"yield"` and `"fertility"`.
#' @export
example_global_dairy <- function() {
  ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")
  inc <- c(CK = 3.06, CM = 30.49, DA = 2.16, DYS = 5.99, LAM = 25.45, MET = 9.57,
           MF = 2.41, OC = 11.46, PTB = 10.01, RP = 12.35, SCK = 47.89, SCM = 40.97)
  type <- ifelse(ids == "PTB", "prevalence", "incidence_rate")
  diseases <- cm_diseases(ids, inc[ids] / 100, type = type, time_horizon = "lactation",
                          source = "Rasmussen et al. 2024, Table 2")
  or <- c("CK:CM" = 2.13, "CK:LAM" = 1.65, "CK:MF" = 1.60, "CK:OC" = 1.97,
          "CK:RP" = 1.55, "CK:SCK" = 6.95, "CK:SCM" = 2.40, "CM:LAM" = 2.10,
          "CM:PTB" = 1.89, "CM:RP" = 2.70, "CM:SCK" = 1.64, "CM:SCM" = 3.05,
          "DA:CM" = 3.45, "DA:MF" = 2.50, "DA:RP" = 3.50, "DA:SCK" = 3.87,
          "DA:SCM" = 3.60, "DYS:LAM" = 2.09, "DYS:OC" = 0.40, "DYS:RP" = 2.74,
          "LAM:OC" = 2.63, "LAM:PTB" = 2.70, "LAM:RP" = 1.50, "LAM:SCK" = 2.01,
          "MET:CK" = 2.42, "MET:CM" = 2.30, "MET:DA" = 3.40, "MET:DYS" = 2.95,
          "MET:LAM" = 6.10, "MET:MF" = 1.50, "MET:OC" = 1.94, "MET:RP" = 3.53,
          "MET:SCK" = 1.94, "MF:DYS" = 9.70, "MF:LAM" = 3.60, "MF:RP" = 2.40,
          "RP:OC" = 2.18, "SCK:RP" = 1.52)
  pr <- do.call(rbind, strsplit(names(or), ":", fixed = TRUE))
  associations <- cm_associations(pr[, 1], pr[, 2], unname(or), measure = "OR",
                                   source = "Rasmussen et al. 2024, Table 3")
  yield <- c(CK = 0.43, CM = 3.25, DA = 2.84, DYS = 4.92, LAM = 4.81, MET = 5.61,
             MF = 0.54, OC = 3.75, PTB = 4.30, RP = 4.20, SCK = 8.40, SCM = 6.29)
  fert <- c(CK = 1.45, CM = 8.42, DA = 1.08, DYS = 2.40, LAM = 3.30, MET = 14.67,
            MF = 2.41, OC = 9.69, PTB = 5.35, RP = 6.76, SCK = 1.12, SCM = 0.26)
  impacts <- cm_impacts(
    disease = c(ids, ids), value = unname(c(yield[ids], fert[ids])),
    outcome = rep(c("yield", "fertility"), each = length(ids)),
    scale = "percent",
    direction = rep(c("decrease", "increase"), each = length(ids)),
    source = "Rasmussen et al. 2024, Table 4"
  )
  cm_model(diseases, associations, impacts, missing_associations = "independent")
}
