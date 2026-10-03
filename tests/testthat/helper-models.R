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

# A model on the supplement population with the given impacts.
supp_model <- function(values, estimand = "crude", adjusted_for = NA_character_,
                       interactions = NULL, units = NULL) {
  cm_model(supp_population(),
           cm_impacts(ids3, values, estimand = estimand, adjusted_for = adjusted_for,
                      units = units),
           interactions)
}

# Supplement population with culling hazard ratios 1.5, 2.0 and 1.3.
supp_hr_model <- function() {
  cm_hr_model(supp_population(), cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = "snapshot_crude"))
}

# Bob's infeasible triple: p = 0.5, ORs 20 (a:b), 20 (a:c), 0.05 (b:c).
bad_population <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

# A model where the published approximation divides by zero for d1:
# b1 = m1^2 / (m1 + A12 m2 + A13 m3) with m = (-A12, 1, 0).
pole_model <- function() {
  e <- deconflate(supp_model(c(1, 1, 1)))$conflation$A["d1", "d2"]
  supp_model(c(-e, 1, 0))
}
