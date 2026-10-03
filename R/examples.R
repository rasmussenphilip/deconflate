#' Worked example from the Supplementary File of Rasmussen et al. (2022)
#'
#' Three hypothetical diseases with prevalences 0.10, 0.15 and 0.20, yield
#' impacts of 2.5%, 5% and 7.5% (entered in percent, so results are in
#' percent), and odds ratios of 2 (d1:d2), 1 (d1:d3) and 3 (d2:d3). The
#' observed mean yield in the example is 10,000 units per animal per year.
#'
#' Note: the published productivity gap (10,225; 20, 61 and 143 units) rounds
#' the adjusted impacts to 2%, 4% and 7% before computing the gap. Without
#' rounding, the published method gives 10,215.7 and 21.1, 56.7 and 137.9
#' units.
#'
#' @return A [cm_model()].
#' @export
#' @examples
#' deconflate(example_supplement(), method = "published")
example_supplement <- function() {
  pop <- cm_population(
    cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20), type = "prevalence"),
    cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3), measure = "OR")
  )
  cm_model(pop, cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%",
                           source = "Rasmussen et al. 2022, Table S1"))
}

global_dairy_inputs <- function(inputs) {
  ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")
  if (inputs == "analysis") {
    inc <- c(CK = 3.06287320599544, CM = 30.5081708702075, DA = 2.16266798689032,
             DYS = 6.11107185570655, LAM = 25.4288506511192, MET = 9.58013664473386,
             MF = 2.41829408262443, OC = 11.2683551256031, PTB = 10.3964073833149,
             RP = 12.346845541711, SCK = 47.9052187708235, SCM = 40.9414831921559)
    type <- ifelse(ids == "PTB", "prevalence",
                   ifelse(ids == "SCM", "probability", "incidence_rate"))
    src <- "Rasmussen et al. 2024, analysis code (1st revision)"
    yield <- c(CK = 0.4321944, CM = 3.2499, DA = 2.83693, DYS = 4.919088, LAM = 4.8061,
               MET = 5.613085, MF = 0.5365130859025376, OC = 3.747839, PTB = 4.3,
               RP = 4.198664, SCK = 8.396472, SCM = 6.293184)
    fert <- c(CK = 1.445122, CM = 8.42, DA = 1.082641, DYS = 2.399177, LAM = 3.304898,
              MET = 14.67308, MF = 2.414949, OC = 9.685465, PTB = 5.349124,
              RP = 6.760971, SCK = 1.122037, SCM = 0.2636645)
    hr <- c(CK = 1.5001, CM = 2.3, DA = 2.851179, DYS = 1.258143, LAM = 1.744976,
            MET = 1.116444, MF = 2.999886, OC = 1.62, PTB = 2.310508, RP = 1.599928,
            SCK = 1.92, SCM = 1.449996)
  } else {
    inc <- c(CK = 3.06, CM = 30.49, DA = 2.16, DYS = 5.99, LAM = 25.45, MET = 9.57,
             MF = 2.41, OC = 11.46, PTB = 10.01, RP = 12.35, SCK = 47.89, SCM = 40.97)
    type <- ifelse(ids == "PTB", "prevalence", "incidence_rate")
    src <- "Rasmussen et al. 2024, Tables 2-4"
    yield <- c(CK = 0.43, CM = 3.25, DA = 2.84, DYS = 4.92, LAM = 4.81, MET = 5.61,
               MF = 0.54, OC = 3.75, PTB = 4.30, RP = 4.20, SCK = 8.40, SCM = 6.29)
    fert <- c(CK = 1.45, CM = 8.42, DA = 1.08, DYS = 2.40, LAM = 3.30, MET = 14.67,
              MF = 2.41, OC = 9.69, PTB = 5.35, RP = 6.76, SCK = 1.12, SCM = 0.26)
    hr <- c(CK = 1.50, CM = 2.30, DA = 2.85, DYS = 1.26, LAM = 1.74, MET = 1.05,
            MF = 3.00, OC = 1.62, PTB = 2.31, RP = 1.60, SCK = 1.92, SCM = 1.45)
  }
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
  pop <- cm_population(
    cm_diseases(ids, unname(inc[ids]) / 100, type = type, time_horizon = "lactation", source = src),
    cm_associations(pr[, 1], pr[, 2], unname(or), measure = "OR",
                    source = "Rasmussen et al. 2024, Table 3")
  )
  list(ids = ids, pop = pop, src = src, yield = yield, fert = fert, hr = hr)
}

#' Global dairy inputs from Rasmussen et al. (2024)
#'
#' The 12 diseases of Rasmussen et al. (2024): global (herd-size weighted)
#' lactational incidence (prevalence for PTB), pooled inter-disease odds
#' ratios (Table 3; pairs not listed are independent), and raw impacts on
#' yield (% decrease) and fertility (% increase in calving interval), at
#' their central values (means of normal distributions, modes of PERT
#' distributions). Each impact type is one analysis.
#'
#' @section Two versions of the inputs:
#' * `inputs = "analysis"` (default) uses the inputs of the published
#'   analysis code (1st revision), which reproduce Table 5: disease
#'   probabilities fixed at `1 - exp(-incidence)` at the unrounded global
#'   means, except subclinical mastitis (entered without conversion, 0.4094
#'   rather than 0.3360, as in the analysis), and unrounded impacts.
#' * `inputs = "tables"` uses Tables 2-4 as printed, with incidence converted
#'   for every disease except PTB.
#'
#' @section Culling:
#' The culling hazard ratios are available as a hazard-ratio model from
#' [example_global_dairy_hr()] (see [deconflate_hr()]). With `culling = TRUE`
#' the analyses also include `culling_hr_minus_1`: hazard ratios minus 1
#' treated as additive impacts, which is how the 2024 analysis adjusted them
#' (with the published method). It is included only to reproduce Table 5.
#'
#' @param inputs `"analysis"` or `"tables"`.
#' @param culling Include the legacy `culling_hr_minus_1` analysis?
#' @return A [cm_analyses()] object.
#' @export
#' @examples
#' gd <- example_global_dairy()
#' deconflate(gd, method = "published")$yield$adjusted
example_global_dairy <- function(inputs = c("analysis", "tables"), culling = FALSE) {
  inputs <- match.arg(inputs)
  g <- global_dairy_inputs(inputs)
  yield <- cm_impacts(g$ids, unname(g$yield[g$ids]), label = "milk yield loss",
                      units = "% decrease", source = g$src)
  fert <- cm_impacts(g$ids, unname(g$fert[g$ids]), label = "calving interval increase",
                     units = "% increase", source = g$src)
  if (culling) {
    cull <- cm_impacts(g$ids, unname(g$hr[g$ids]) - 1,
                       label = "culling hazard ratio minus 1 (2024 analysis; legacy)",
                       units = "hazard ratio - 1", source = g$src)
    cm_analyses(g$pop, yield = yield, fertility = fert, culling_hr_minus_1 = cull)
  } else {
    cm_analyses(g$pop, yield = yield, fertility = fert)
  }
}

#' Global dairy culling hazard ratios (Rasmussen et al. 2024)
#'
#' The culling hazard ratios of Rasmussen et al. (2024), Table 4 (or the
#' analysis inputs; see [example_global_dairy()]) as a hazard-ratio model for
#' [deconflate_hr()] and [attributable_risk()]. With `inputs = "analysis"`,
#' metritis uses the analysis value (PERT mode 1.116) rather than Table 4
#' (normal, mean 1.05).
#'
#' @param inputs `"analysis"` or `"tables"`.
#' @return A [cm_hr_model()].
#' @export
#' @examples
#' hr <- example_global_dairy_hr()
#' deconflate_hr(hr, method = "published")$adjusted
example_global_dairy_hr <- function(inputs = c("analysis", "tables")) {
  inputs <- match.arg(inputs)
  g <- global_dairy_inputs(inputs)
  cm_hr_model(g$pop, cm_hazard_ratios(g$ids, unname(g$hr[g$ids]), source = g$src))
}

#' UK dairy example from Rasmussen et al. (2022)
#'
#' Thirteen endemic diseases and conditions of UK dairy cattle, with cow-level
#' prevalence (Table 2), inter-disease odds ratios (Table 3, all other pairs
#' independent), and three analyses: milk yield (Table 4, % decrease),
#' calving interval (Table 5, % increase) and culling (Table 6). Use
#' [uk_dairy_2022_economics()] for the observed means and unit values
#' (Table 1), and [reproduce_rasmussen_2022()] for the published tables.
#'
#' The culling analysis reproduces the paper's approach, which converted
#' hazard ratios to excess annual culling risks by treating them as odds
#' ratios. That conversion is not supported for new analyses (use
#' [deconflate_hr()]); the hazard ratios are attached as attribute
#' `"hazard_ratios"` for that purpose.
#'
#' Disease ids: CO cystic ovary, DA displaced abomasum, DYS dystocia, FAS
#' fasciolosis, GIN gastrointestinal nematodes, LAM lameness, MAS mastitis,
#' MET metritis, MF milk fever, NEO neosporosis, PTB paratuberculosis, RP
#' retained placenta, SCK subclinical ketosis.
#'
#' @param yield_sck Raw yield impact of subclinical ketosis, in percent (the
#'   printed 3.05, or 100 * 340 / 8737 = 3.89 from note i of Table 4).
#' @return A [cm_analyses()] object with analyses `yield`, `fertility` and
#'   `culling`, and attribute `"hazard_ratios"` (a [cm_hazard_ratios()]).
#' @export
#' @examples
#' uk <- example_uk_dairy_2022()
#' deconflate(uk, method = "published")
example_uk_dairy_2022 <- function(yield_sck = 3.05) {
  ids <- c("CO", "DA", "DYS", "FAS", "GIN", "LAM", "MAS", "MET", "MF", "NEO", "PTB", "RP", "SCK")
  prev <- c(CO = 0.09, DA = 0.03, DYS = 0.02, FAS = 0.10, GIN = 0.21, LAM = 0.30, MAS = 0.30,
            MET = 0.10, MF = 0.08, NEO = 0.15, PTB = 0.07, RP = 0.05, SCK = 0.22)
  or <- c("RP:MET" = 6.20, "DA:SCK" = 4.25, "RP:DYS" = 4.10, "MET:DA" = 3.40,
          "MET:DYS" = 3.20, "LAM:PTB" = 2.70, "MF:DA" = 2.50, "MAS:MET" = 2.30,
          "DA:RP" = 2.20, "MAS:DA" = 2.10, "SCK:MF" = 2.10, "LAM:SCK" = 2.01,
          "MAS:MF" = 1.90, "MAS:PTB" = 1.89, "MAS:CO" = 1.65, "MAS:SCK" = 1.64,
          "SCK:CO" = 1.60, "MET:SCK" = 1.40, "SCK:RP" = 1.20)
  pr <- do.call(rbind, strsplit(names(or), ":", fixed = TRUE))
  pop <- cm_population(
    cm_diseases(ids, unname(prev[ids]), type = "prevalence", time_horizon = "year",
                reference_population = "UK dairy cows", source = "Rasmussen et al. 2022, Table 2"),
    cm_associations(pr[, 1], pr[, 2], unname(or), measure = "OR",
                    source = "Rasmussen et al. 2022, Table 3")
  )
  yield <- c(CO = 0, DA = 4.04, DYS = 4.05, FAS = 7.33, GIN = 3.28, LAM = 5.54, MAS = 4.57,
             MET = 3.95, MF = 0.41, NEO = 4.20, PTB = 5.90, RP = 7.38, SCK = yield_sck)
  fert <- c(CO = 11.26, DA = 0, DYS = 6.96, FAS = 0, GIN = 1.20, LAM = 12.47, MAS = 0,
            MET = 4.74, MF = 0, NEO = 7.21, PTB = 5.79, RP = 2.74, SCK = 1.50)
  hr <- c(DA = 3.83, LAM = 3.40, MAS = 2.78, MF = 2.50, PTB = 2.40, MET = 2.20,
          SCK = 2.10, DYS = 1.90, NEO = 1.60)
  hr_all <- stats::setNames(rep(1, length(ids)), ids)
  hr_all[names(hr)] <- hr
  excess <- legacy_hr_as_or_excess(unname(hr_all), unname(prev[ids]), 0.27)
  out <- cm_analyses(
    pop,
    yield = cm_impacts(ids, unname(yield[ids]), label = "milk yield loss", units = "% decrease",
                       source = "Rasmussen et al. 2022, Table 4"),
    fertility = cm_impacts(ids, unname(fert[ids]), label = "calving interval increase",
                           units = "% increase", source = "Rasmussen et al. 2022, Table 5"),
    culling = cm_impacts(ids, excess,
                         label = "excess annual culling risk (hazard ratio treated as an odds ratio; legacy)",
                         units = "proportional increase in the culling rate",
                         source = "Rasmussen et al. 2022, Table 6 (legacy conversion)")
  )
  attr(out, "hazard_ratios") <- cm_hazard_ratios(ids, unname(hr_all[ids]),
                                                 source = "Rasmussen et al. 2022, Table 6")
  out
}

#' Economic inputs for the UK dairy example (Rasmussen et al. 2022, Table 1)
#'
#' Valuation inputs for each analysis of [example_uk_dairy_2022()], for
#' [productivity_gap()] and [value_losses()] (or the `valuation` argument of
#' [summary.cm_result()] and [compare_methods()]): milk yield 8737 kg/cow/year
#' valued at GBP 0.3022/kg; calving interval 401 days, each day valued at
#' lifetime daily yield (13 kg) times the milk price; culling rate 27% per
#' year (impacts as proportional increases, as in the paper) valued at the
#' replacement price per percentage point (GBP 13.3536). Veterinary
#' expenditure of GBP 71.09 per cow per year is a separate lump sum.
#'
#' @return A list with `valuation` (one valuation list per analysis) and
#'   `additional` (lump-sum costs).
#' @export
uk_dairy_2022_economics <- function() {
  price <- 30.22 / 100
  list(valuation = list(
    yield = list(observed = 8737, direction = "decrease", effect = "percent", unit_value = price),
    fertility = list(observed = 401, direction = "increase", effect = "percent",
                     unit_value = 13 * price),
    culling = list(observed = 27, direction = "increase", effect = "proportion",
                   unit_value = 1335.36 / 100)),
    additional = c(veterinary = 71.09))
}

#' Monte Carlo sampler for the global dairy inputs (Rasmussen et al. 2024)
#'
#' A [cm_batch_sampler()] for [example_global_dairy()]: shared draws of the
#' odds ratios (and, with `inputs = "tables"`, the incidences), and each
#' analysis's impact distributions. PERT distributions use the reported
#' central value as the mode (shape 4), and normal distributions of odds
#' ratios are truncated at zero.
#'
#' With `inputs = "analysis"` the disease probabilities are fixed and the
#' impact distributions use the unrounded parameters of the analysis code.
#' For `culling_hr_minus_1`, the analysis entered HR - 1 and scaled the
#' standard deviations of normal distributions by (HR - 1) / HR; this is kept
#' for reproduction. With `inputs = "tables"`, incidences are drawn from the
#' global distributions of Table 2 and Tables 3-4 are used as printed.
#'
#' With `method = "published"`, Monte Carlo means reproduce Table 5 (see
#' [reproduce_rasmussen_2024()] and `vignette("reproducing-published")`).
#'
#' @param inputs `"analysis"` or `"tables"`.
#' @param culling Include the legacy `culling_hr_minus_1` analysis?
#' @return A [cm_batch_sampler()].
#' @export
#' @examples
#' \donttest{
#' mc <- cm_monte_carlo(sampler_global_dairy(), 200, method = "published", seed = 1)
#' s <- summary(mc, diagnose = FALSE)
#' s[s$analysis == "yield", c("disease", "mean")]
#' }
sampler_global_dairy <- function(inputs = c("analysis", "tables"), culling = FALSE) {
  inputs <- match.arg(inputs)
  n0 <- function(m, s) dist_normal(m, s, lower = 0)
  pe <- function(mode, mn, mx) dist_pert(mn, mode, mx)
  nn <- function(m, s) dist_normal(m, s)
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
  if (inputs == "analysis") {
    diseases <- list()
    yield <- list(
      CK = pe(0.4321944, 0.2351676, 1.043482), CM = nn(3.2499, 0.7584468),
      DA = pe(2.83693, -1.451884, 9.190728), DYS = nn(4.919088, 0.9688162),
      LAM = nn(4.8061, 0.8651519), MET = nn(5.613085, 1.350757),
      OC = pe(3.747839, 1.711971, 4.326951), PTB = nn(4.3, 0.6683673),
      RP = nn(4.198664, 1.154555), SCK = nn(8.396472, 1.185384), SCM = nn(6.293184, 1.200231))
    fert <- list(
      CK = nn(1.445122, 0.3624052), CM = nn(8.42, 2.424912), DA = nn(1.082641, 2.036204),
      DYS = nn(2.399177, 0.9317402), LAM = pe(3.304898, 1.190476, 10.71429),
      MET = nn(14.67308, 8.535001), MF = pe(2.414949, 2.032968, 3.095238),
      OC = pe(9.685465, 5.043478, 21.42857), PTB = nn(5.349124, 2.527459),
      RP = nn(6.760971, 1.558483), SCK = nn(1.122037, 1.822338),
      SCM = pe(0.2636645, -0.1242236, 5.681529))
    cull <- list(
      CK = nn(0.5001, 0.10000976530231317), CM = nn(1.3, 0.17342740434782608),
      DA = pe(1.851179, 0, 6.9), DYS = pe(0.258143, -0.4, 1.1),
      LAM = nn(0.744976, 0.07449644729921788), MET = pe(0.116444, -0.4, 0.5),
      MF = nn(1.999886, 0.6012380424926813), OC = nn(0.62, 0.1591395716049383),
      PTB = nn(1.310508, 0.2011405862857865), RP = nn(0.599928, 0.11954628295273287),
      SCK = nn(0.92, 0.0855654625), SCM = nn(0.449996, 0.07758647609400302))
  } else {
    b <- function(a, bb) dist_beta(a, bb)
    diseases <- list(
      CK = b(13.95, 441.59), CM = b(8.55, 18.20), DA = b(27.50, 1245.11),
      DYS = dist_pert(0.0190, 0.0599, 0.1080), LAM = b(78.29, 227.42),
      MET = b(17.73, 167.35), MF = b(7.76, 313.41), OC = dist_pert(0.0270, 0.1146, 0.1907),
      PTB = dist_pert(0.0119, 0.1001, 0.2108), RP = b(33.75, 239.46),
      SCK = b(178.14, 193.73), SCM = b(116.23, 167.02))
    yield <- list(
      CK = pe(0.43, 0.24, 1.04), CM = nn(3.25, 0.76), DA = pe(2.84, -1.45, 9.19),
      DYS = nn(4.92, 0.97), LAM = nn(4.81, 0.87), MET = nn(5.61, 1.35),
      OC = pe(3.75, 1.71, 4.33), PTB = nn(4.30, 0.67), RP = nn(4.20, 1.15),
      SCK = nn(8.40, 1.19), SCM = nn(6.29, 1.20))
    fert <- list(
      CK = nn(1.45, 0.36), CM = nn(8.42, 2.42), DA = nn(1.08, 2.04), DYS = nn(2.40, 0.93),
      LAM = pe(3.30, 1.19, 10.71), MET = nn(14.67, 8.54), MF = pe(2.41, 2.03, 3.10),
      OC = pe(9.69, 5.04, 21.43), PTB = nn(5.35, 2.53), RP = nn(6.76, 1.56),
      SCK = nn(1.12, 1.82), SCM = pe(0.26, -0.12, 5.68))
    cull <- list(
      CK = nn(0.50, 0.30), CM = nn(1.30, 0.31), DA = pe(1.85, 0, 6.90), DYS = pe(0.26, -0.40, 1.10),
      LAM = nn(0.74, 0.17), MET = nn(0.05, 0.15), MF = nn(2.00, 0.90), OC = nn(0.62, 0.42),
      PTB = nn(1.31, 0.35), RP = nn(0.60, 0.32), SCK = nn(0.92, 0.18), SCM = nn(0.45, 0.25))
  }
  imp <- list(yield = yield, fertility = fert)
  if (culling) imp$culling_hr_minus_1 <- cull
  cm_batch_sampler(example_global_dairy(inputs = inputs, culling = culling),
                   diseases = diseases, associations = associations, impacts = imp)
}
