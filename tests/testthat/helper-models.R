# Shared test inputs.

ids3 <- c("d1", "d2", "d3")

# Supplementary File population (Rasmussen et al. 2022): P = 0.10, 0.15,
# 0.20; OR d1:d2 = 2, d1:d3 = 1, d2:d3 = 3.
supp_population <- function(three_way = NULL) {
  cm_population(
    cm_diseases(ids3, c(0.10, 0.15, 0.20)),
    cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)),
    three_way = three_way
  )
}

# The supplement population without the d1:d3 row: that pair is unknown, so
# deconflate() uses the global method (inst/validation/reference_v040.py, 6b).
supp_unknown_population <- function(three_way = NULL) {
  cm_population(
    cm_diseases(ids3, c(0.10, 0.15, 0.20)),
    cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)),
    three_way = three_way
  )
}

# A model on the supplement population with the given additive impacts.
supp_model <- function(values, estimand = "crude", adjusted_for = NA_character_,
                       interactions = NULL, units = NULL) {
  cm_model(supp_population(),
           cm_impacts(ids3, values, estimand = estimand, adjusted_for = adjusted_for,
                      units = units),
           interactions)
}

# The supplement population with event impacts: by default culling hazard
# ratios 1.5, 2.0 and 1.3 (snapshot_crude). Adjust with
# deconflate(..., event_model = TRUE, overall_risk = ...).
supp_hr_model <- function(values = c(1.5, 2.0, 1.3), measure = "HR",
                          estimand = "snapshot_crude", adjusted_for = NA_character_) {
  cm_model(supp_population(),
           cm_impacts(ids3, values, measure = measure, estimand = estimand,
                      adjusted_for = adjusted_for))
}

# Bob's infeasible triple: p = 0.5, ORs 20 (a:b), 20 (a:c), 0.05 (b:c).
bad_population <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

# A model where the published approximation (internal adjust_impacts(method =
# "published"), used by compare_methods() and the reproduce_*() functions)
# divides by zero for d1: b1 = m1^2 / (m1 + A12 m2 + A13 m3) with
# m = (-A12, 1, 0).
pole_model <- function() {
  e <- deconflate(supp_model(c(1, 1, 1)), n_draws = 0)$conflation$A["d1", "d2"]
  supp_model(c(-e, 1, 0))
}

# Path of a shipped example file (inst/extdata/<folder>/<file>).
extdata <- function(folder, file = NULL) {
  dir <- system.file("extdata", folder, package = "deconflate")
  if (is.null(file)) dir else file.path(dir, file)
}

# The five-disease example (inst/extdata/five_diseases) with one impact
# table, read with cm_read_inputs(); `three_way = TRUE` also reads
# three_way.csv, and `interactions` names an interaction file.
read_five <- function(impacts = "impacts_yield.csv", three_way = FALSE, interactions = NULL) {
  f <- function(x) extdata("five_diseases", x)
  cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                 impacts = f(impacts),
                 interactions = if (!is.null(interactions)) f(interactions),
                 three_way = if (three_way) f("three_way.csv"))
}
