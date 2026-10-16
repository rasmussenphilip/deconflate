# Shared test inputs.

ids3 <- c("d1", "d2", "d3")

# Supplementary File population (Rasmussen et al. 2022): P = 0.10, 0.15,
# 0.20; OR d1:d2 = 2, d1:d3 = 1, d2:d3 = 3. `extra`: more association rows
# (a cm_associations() table), appended.
supp_population <- function(extra = NULL) {
  a <- cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3))
  if (!is.null(extra)) {
    extra$row <- nrow(a) + seq_len(nrow(extra))
    a <- rbind(a, extra)
  }
  cm_population(cm_diseases(ids3, c(0.10, 0.15, 0.20)), a)
}

# The supplement population without the d1:d3 row: that pair is unknown, and
# the model fills it in (inst/validation/reference_v040.py, 6b).
supp_unknown_population <- function(extra = NULL) {
  a <- cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
  if (!is.null(extra)) {
    extra$row <- nrow(a) + seq_len(nrow(extra))
    a <- rbind(a, extra)
  }
  cm_population(cm_diseases(ids3, c(0.10, 0.15, 0.20)), a)
}

# A model on the supplement population with the given additive impacts
# (one row per disease; `adjusted_for` per row), plus optional
# interactions: a named vector of values by combination ("d1;d2"), entered
# as rows adjusted for all (the interaction itself).
supp_model <- function(values, adjusted_for = NA_character_, interactions = NULL, units = NULL,
                       population = supp_population()) {
  d <- c(ids3, names(interactions))
  v <- c(values, unname(interactions))
  adj <- c(rep_len(adjusted_for, 3), rep("all", length(interactions)))
  cm_model(population, cm_impacts(d, v, adjusted_for = adj, units = units))
}

# The supplement population with event impacts: by default culling hazard
# ratios 1.5, 2.0 and 1.3 (crude). Adjust with
# deconflate(..., event_model = TRUE, overall_risk = ...).
supp_hr_model <- function(values = c(1.5, 2.0, 1.3), measure = "HR", adjusted_for = NA_character_,
                          population = supp_population()) {
  cm_model(population, cm_impacts(ids3, values, measure = measure, adjusted_for = adjusted_for))
}

# Bob's infeasible triple: p = 0.5, ORs 20 (a:b), 20 (a:c), 0.05 (b:c).
bad_population <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

# Path of a shipped example file (inst/extdata/<folder>/<file>).
extdata <- function(folder, file = NULL) {
  dir <- system.file("extdata", folder, package = "deconflate")
  if (is.null(file)) dir else file.path(dir, file)
}

# The five-disease example (inst/extdata/five_diseases) with one impact
# table, read with cm_read_inputs().
read_five <- function(impacts = "impacts_yield.csv") {
  f <- function(x) extdata("five_diseases", x)
  suppressMessages(cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                                  impacts = f(impacts)))
}

# The five-disease example with only its ten crude pair rows (row 1 of the
# two LAM-MAS rows, no OR_ratio row): the layout of earlier versions.
read_five_pairs <- function(impacts = "impacts_yield.csv") {
  f <- function(x) extdata("five_diseases", x)
  a <- utils::read.csv(f("associations.csv"), stringsAsFactors = FALSE)
  a <- a[-c(2, 12), ]
  suppressMessages(cm_read_inputs(diseases = f("diseases.csv"), associations = a, impacts = f(impacts)))
}

# Give every pair of a population without an association an odds ratio of 1
# (as the published analyses did with pairs they had no estimate for).
with_independent_pairs <- function(pop) {
  pt <- pair_tables(pop)
  u <- pt[pt$status == "unknown", , drop = FALSE]
  if (!nrow(u)) return(pop)
  add <- cm_associations(u$disease1, u$disease2, 1, measure = "OR")
  a <- pop$associations
  if (!is.null(a)) add$row <- max(a$row) + seq_len(nrow(add))
  pop$associations <- if (is.null(a)) add else rbind(a, add)
  pop
}

# A value of the long totals table of summary(): total_value(s, "reduction").
total_value <- function(s, quantity) s$totals$value[s$totals$quantity == quantity]

# Printed output as one line with single spaces (printing wraps long lines).
printed <- function(x, ...) gsub("\\s+", " ", paste(utils::capture.output(print(x, ...)), collapse = " "))

# The five-disease culling model, adjusted without draws.
cull_five <- function(...) deconflate(read_five("culling.csv"), event_model = TRUE, n_draws = 0, ...)
