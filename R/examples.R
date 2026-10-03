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

#' UK dairy example from Rasmussen et al. (2022)
#'
#' Thirteen endemic diseases and conditions of UK dairy cattle, with cow-level
#' prevalence (Table 2), inter-disease odds ratios (Table 3, all other pairs
#' independent), and impacts on milk yield (Table 4, % decrease), calving
#' interval (Table 5, % increase) and culling (Table 6, hazard ratios
#' converted to excess annual culling risk). Use [uk_dairy_2022_economics()]
#' for the observed means and unit values (Table 1).
#'
#' Disease ids: CO cystic ovary, DA displaced abomasum, DYS dystocia, FAS
#' fasciolosis, GIN gastrointestinal nematodes, LAM lameness, MAS mastitis,
#' MET metritis, MF milk fever, NEO neosporosis, PTB paratuberculosis, RP
#' retained placenta, SCK subclinical ketosis.
#'
#' @section Reproducibility of the published tables:
#' With the defaults (`culling_method = "or_approx"`, `culling_scale =
#' "proportion"`) and `method = "published"`:
#' * fertility (calving interval) reproduces Tables 8-10 (disease-free value
#'   375.09 days; gap value GBP 101.79 per cow per year);
#' * culling reproduces the adjusted hazard ratios of Table 8 to within 0.01
#'   and the gap to within 0.5% (22.55% vs 22.56%; GBP 59.48 vs 59.27);
#' * yield does not reproduce Tables 8-10 exactly from the printed Table 4.
#'   The tables imply different yield inputs, e.g. subclinical ketosis
#'   340 kg / 8737 kg = 3.89% (note i of Table 4) rather than the printed
#'   3.05%, and Table 10 implies further differences. Use
#'   `yield_sck = 100 * 340 / 8737` to apply the 3.89% value.
#'
#' @param culling_method Hazard-ratio conversion: `"or_approx"` (published)
#'   or `"proportional_hazards"`.
#' @param culling_scale `"proportion"` (published) or `"absolute"`; see
#'   [as_impacts()].
#' @param yield_sck Raw yield impact of subclinical ketosis, in percent.
#' @return A [cm_model()] with outcomes `"yield"`, `"fertility"` and
#'   `"culling"`. The hazard-ratio conversion is attached as attribute
#'   `"culling_conversion"`, for [adjusted_hr()].
#' @export
#' @examples
#' m <- example_uk_dairy_2022()
#' eco <- uk_dairy_2022_economics()
#' res <- deconflate(m, method = "published")
#' gaps <- productivity_gap(res, eco$observed)
#' gaps$summary
#' value_losses(gaps, eco$unit_value, eco$additional)$total
#' adjusted_hr(res, attr(m, "culling_conversion"), method = "published")
example_uk_dairy_2022 <- function(culling_method = c("or_approx", "proportional_hazards"),
                                  culling_scale = c("proportion", "absolute"),
                                  yield_sck = 3.05) {
  culling_method <- match.arg(culling_method)
  culling_scale <- match.arg(culling_scale)
  ids <- c("CO", "DA", "DYS", "FAS", "GIN", "LAM", "MAS", "MET", "MF", "NEO", "PTB", "RP", "SCK")
  prev <- c(CO = 0.09, DA = 0.03, DYS = 0.02, FAS = 0.10, GIN = 0.21, LAM = 0.30, MAS = 0.30,
            MET = 0.10, MF = 0.08, NEO = 0.15, PTB = 0.07, RP = 0.05, SCK = 0.22)
  diseases <- cm_diseases(ids, prev[ids], type = "prevalence", time_horizon = "year",
                          reference_population = "UK dairy cows",
                          source = "Rasmussen et al. 2022, Table 2")
  or <- c("RP:MET" = 6.20, "DA:SCK" = 4.25, "RP:DYS" = 4.10, "MET:DA" = 3.40,
          "MET:DYS" = 3.20, "LAM:PTB" = 2.70, "MF:DA" = 2.50, "MAS:MET" = 2.30,
          "DA:RP" = 2.20, "MAS:DA" = 2.10, "SCK:MF" = 2.10, "LAM:SCK" = 2.01,
          "MAS:MF" = 1.90, "MAS:PTB" = 1.89, "MAS:CO" = 1.65, "MAS:SCK" = 1.64,
          "SCK:CO" = 1.60, "MET:SCK" = 1.40, "SCK:RP" = 1.20)
  pr <- do.call(rbind, strsplit(names(or), ":", fixed = TRUE))
  associations <- cm_associations(pr[, 1], pr[, 2], unname(or), measure = "OR",
                                  source = "Rasmussen et al. 2022, Table 3")
  yield <- c(CO = 0, DA = 4.04, DYS = 4.05, FAS = 7.33, GIN = 3.28, LAM = 5.54, MAS = 4.57,
             MET = 3.95, MF = 0.41, NEO = 4.20, PTB = 5.90, RP = 7.38, SCK = yield_sck)
  fert <- c(CO = 11.26, DA = 0, DYS = 6.96, FAS = 0, GIN = 1.20, LAM = 12.47, MAS = 0,
            MET = 4.74, MF = 0, NEO = 7.21, PTB = 5.79, RP = 2.74, SCK = 1.50)
  hr <- c(DA = 3.83, LAM = 3.40, MAS = 2.78, MF = 2.50, PTB = 2.40, MET = 2.20,
          SCK = 2.10, DYS = 1.90, NEO = 1.60)
  conv <- hr_conversion(diseases, hr, overall_risk = 0.27, method = culling_method)
  impacts <- combine_impacts(
    cm_impacts(ids, unname(yield[ids]), outcome = "yield", scale = "percent",
               units = "kg/cow/year", direction = "decrease",
               source = "Rasmussen et al. 2022, Table 4"),
    cm_impacts(ids, unname(fert[ids]), outcome = "fertility", scale = "percent",
               units = "days", direction = "increase",
               source = "Rasmussen et al. 2022, Table 5"),
    as_impacts(conv, outcome = "culling", scale = culling_scale,
               source = "Rasmussen et al. 2022, Table 6")
  )
  m <- cm_model(diseases, associations, impacts, missing_associations = "independent")
  attr(m, "culling_conversion") <- conv
  m
}

#' Economic inputs for the UK dairy example (Rasmussen et al. 2022, Table 1)
#'
#' Observed means and unit values for [productivity_gap()] and
#' [value_losses()]: milk yield 8737 kg/cow/year valued at GBP 0.3022/kg;
#' calving interval 401 days, each day valued at lifetime daily yield
#' (13 kg) times the milk price; culling rate 27% per year valued at the
#' replacement price (GBP 1335.36); and veterinary expenditure of GBP 71.09
#' per cow per year added as a lump sum.
#'
#' @param culling_scale Must match the `culling_scale` used in
#'   [example_uk_dairy_2022()]: with `"proportion"` the culling rate is in
#'   percent (27) and valued per percentage point; with `"absolute"` it is a
#'   proportion (0.27) valued per unit.
#' @return A list with `observed`, `unit_value` and `additional`.
#' @export
uk_dairy_2022_economics <- function(culling_scale = c("proportion", "absolute")) {
  culling_scale <- match.arg(culling_scale)
  price <- 30.22 / 100
  if (culling_scale == "proportion") {
    observed <- c(yield = 8737, fertility = 401, culling = 27)
    cull_value <- 1335.36 / 100
  } else {
    observed <- c(yield = 8737, fertility = 401, culling = 0.27)
    cull_value <- 1335.36
  }
  list(observed = observed,
       unit_value = c(yield = price, fertility = 13 * price, culling = cull_value),
       additional = c(veterinary = 71.09))
}

#' Monte Carlo sampler for the global dairy inputs (Rasmussen et al. 2024)
#'
#' Input distributions from Tables 2-4 of Rasmussen et al. (2024) for
#' [example_global_dairy()]: global beta distributions for incidence, PERT
#' distributions (the reported central value is used as the mode), normal
#' distributions for odds ratios (truncated at zero) and impacts.
#'
#' Running [cm_monte_carlo()] with `method = "published"` and averaging the
#' adjusted yield impacts over draws reproduces Table 5 closely. Fertility
#' impacts whose raw distributions extend below zero (e.g. metritis,
#' subclinical ketosis, subclinical mastitis) are unstable under the
#' published approximation, and their means are sensitive to how such draws
#' are handled.
#'
#' @return A [cm_sampler()].
#' @export
#' @examples
#' \donttest{
#' mc <- cm_monte_carlo(sampler_global_dairy(), 200, method = "published", seed = 1)
#' s <- summary(mc)
#' s[s$outcome == "yield", c("disease", "mean")]
#' }
sampler_global_dairy <- function() {
  b <- function(a, bb) dist_beta(a, bb)
  diseases <- list(
    CK = b(13.95, 441.59), CM = b(8.55, 18.20), DA = b(27.50, 1245.11),
    DYS = dist_pert(0.0190, 0.0599, 0.1080), LAM = b(78.29, 227.42),
    MET = b(17.73, 167.35), MF = b(7.76, 313.41), OC = dist_pert(0.0270, 0.1146, 0.1907),
    PTB = dist_pert(0.0119, 0.1001, 0.2108), RP = b(33.75, 239.46),
    SCK = b(178.14, 193.73), SCM = b(116.23, 167.02)
  )
  n0 <- function(m, s) dist_normal(m, s, lower = 0)
  pe <- function(mode, mn, mx) dist_pert(mn, mode, mx)
  associations <- list(
    "CK:CM" = pe(2.13, 1.20, 3.40), "CK:LAM" = pe(1.65, 1.20, 2.40), "CK:MF" = n0(1.60, 0.13),
    "CK:OC" = pe(1.97, 1.30, 4.10), "CK:RP" = pe(1.55, 1.00, 1.90), "CK:SCK" = n0(6.95, 1.28),
    "CK:SCM" = n0(2.40, 0.41), "CM:PTB" = n0(1.89, 0.20), "CM:RP" = n0(2.70, 0.33),
    "CM:SCK" = n0(1.64, 0.20), "CM:SCM" = pe(3.05, 1.30, 6.50), "DA:CM" = pe(3.45, 1.40, 4.80),
    "DA:MF" = n0(2.50, 0.48), "DA:RP" = pe(3.50, 1.60, 4.60), "DA:SCK" = n0(3.87, 0.34),
    "DA:SCM" = n0(3.60, 1.35), "DYS:LAM" = n0(2.09, 0.26), "DYS:RP" = pe(2.74, 1.25, 5.96),
    "LAM:OC" = n0(2.63, 1.44), "LAM:PTB" = n0(2.70, 1.22), "LAM:RP" = n0(1.50, 0.31),
    "LAM:SCK" = n0(2.01, 0.20), "MET:CK" = pe(2.42, 1.20, 10.40), "MET:CM" = pe(2.30, 1.20, 3.80),
    "MET:DA" = pe(3.40, 1.60, 7.60), "MET:DYS" = pe(2.95, 0.98, 9.72), "MET:LAM" = n0(6.10, 1.45),
    "MET:MF" = n0(1.50, 0.15), "MET:OC" = pe(1.94, 1.20, 3.00), "MET:RP" = pe(3.53, 1.80, 6.52),
    "MET:SCK" = n0(1.94, 0.09), "MF:DYS" = n0(9.70, 1.30), "MF:RP" = n0(2.40, 0.20),
    "RP:OC" = pe(2.18, 1.78, 2.57), "SCK:RP" = n0(1.52, 0.19)
  )
  nn <- function(m, s) dist_normal(m, s)
  impacts <- list(
    "yield:CK" = pe(0.43, 0.24, 1.04), "yield:CM" = nn(3.25, 0.76), "yield:DA" = pe(2.84, -1.45, 9.19),
    "yield:DYS" = nn(4.92, 0.97), "yield:LAM" = nn(4.81, 0.87), "yield:MET" = nn(5.61, 1.35),
    "yield:OC" = pe(3.75, 1.71, 4.33), "yield:PTB" = nn(4.30, 0.67), "yield:RP" = nn(4.20, 1.15),
    "yield:SCK" = nn(8.40, 1.19), "yield:SCM" = nn(6.29, 1.20),
    "fertility:CK" = nn(1.45, 0.36), "fertility:CM" = nn(8.42, 2.42), "fertility:DA" = nn(1.08, 2.04),
    "fertility:DYS" = nn(2.40, 0.93), "fertility:LAM" = pe(3.30, 1.19, 10.71),
    "fertility:MET" = nn(14.67, 8.54), "fertility:MF" = pe(2.41, 2.03, 3.10),
    "fertility:OC" = pe(9.69, 5.04, 21.43), "fertility:PTB" = nn(5.35, 2.53),
    "fertility:RP" = nn(6.76, 1.56), "fertility:SCK" = nn(1.12, 1.82),
    "fertility:SCM" = pe(0.26, -0.12, 5.68)
  )
  cm_sampler(example_global_dairy(), diseases = diseases, associations = associations,
             impacts = impacts)
}
