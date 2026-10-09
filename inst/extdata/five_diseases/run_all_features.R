# deconflate: a tour of every feature with the five-disease example
#
# The CSV files in this folder describe an illustrative dairy population
# with five diseases (LAM lameness, MAS mastitis, MET metritis, SCK
# subclinical ketosis, RP retained placenta) and four impact tables: yield
# (% of yield), calving interval (days), welfare (score points, with
# interactions) and culling (event impacts: hazard ratios, a risk ratio, an
# odds ratio and a risk difference). The values are made up to exercise the
# package; they are not estimates. See README.md for what each file
# contains.
#
# Step through the script section by section, or run it all and see every
# result with
#   source(system.file("extdata", "five_diseases", "run_all_features.R",
#                      package = "deconflate"), echo = TRUE, max.deparse.length = Inf)
# (plain source() runs it silently). Sections with draws use small numbers of
# draws so that the script runs in a few minutes; the default of
# deconflate() is 1000.
#
# To keep the results, set a folder before running the script:
#   out_dir <- "C:/Users/me/Desktop/five_diseases_output"
#   source(..., echo = TRUE, max.deparse.length = Inf)
# The printed output is then written to run_all_features_output.txt (and
# still shown in the console) and all plots to run_all_features_plots.pdf in
# that folder. Messages and warnings appear in the console only. If the
# script stops with an error, run sink() to close the output file.

library(deconflate)

five <- system.file("extdata", "five_diseases", package = "deconflate")
if (!nzchar(five)) five <- "."  # running from a copy of this folder
f <- function(file) file.path(five, file)

# Output folder: set `out_dir` before running the script to save the results
# (see above). Without it, plots go to the screen, or to a PDF in tempdir()
# when the script runs non-interactively.
if (!exists("out_dir")) out_dir <- NULL
if (!is.null(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  text_file <- file.path(out_dir, "run_all_features_output.txt")
  sink(text_file, split = TRUE)
}
save_plots <- !is.null(out_dir) || !interactive()
if (save_plots) {
  plot_file <- file.path(if (is.null(out_dir)) tempdir() else out_dir, "run_all_features_plots.pdf")
  grDevices::pdf(plot_file, width = 9, height = 6)
}

section <- function(title) cat("\n\n==== ", title, " ====\n\n", sep = "")

# Run an expression and report a deconflate error instead of stopping.
# Several calls below fail on purpose, to show which combinations are not
# supported.
try_show <- function(expr) {
  tryCatch(expr, deconflate_error = function(e) {
    cat("Error (", paste(setdiff(class(e), c("error", "condition")), collapse = ", "), "): ",
        conditionMessage(e), "\n", sep = "")
    invisible(NULL)
  })
}

# Draw a plot, or say why it could not be drawn (e.g. a plot pane too small
# for several panels).
try_plot <- function(expr) {
  tryCatch(expr, error = function(e) message("Plot skipped: ", conditionMessage(e)))
}

# ---- 1. Read and check the inputs -------------------------------------------
section("1. Reading the inputs")

# Every table is its own argument, with any file name. One impact table per
# model: the disease, association and three-way tables are shared.
read_five <- function(impacts, interactions = NULL, three_way = TRUE,
                      associations = f("associations.csv")) {
  cm_read_inputs(diseases = f("diseases.csv"), associations = associations,
                 impacts = f(impacts),
                 interactions = if (!is.null(interactions)) f(interactions),
                 three_way = if (three_way) f("three_way.csv"))
}

# No problems are expected (a note is printed if any).
print(cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                      impacts = f("impacts_yield.csv"), three_way = f("three_way.csv")))
m_yield <- read_five("impacts_yield.csv", three_way = FALSE)
m_yield
m_yield$diseases                # MAS is an incidence rate, converted to a probability
pair_tables(m_yield)            # the 2x2 table of every pair, from each measure
names(m_yield$distributions)    # the distributions given in the tables

m_yield3 <- read_five("impacts_yield.csv")                            # with the three-way term
m_ci <- read_five("impacts_calving_interval.csv", three_way = FALSE)
m_welfare <- read_five("impacts_welfare.csv", interactions = "interactions_welfare.csv")
m_cull <- read_five("culling.csv")
m_yield3$three_way              # a three-way scenario
m_cull$impacts                  # event impacts: a measure column

# ---- 2. Feasibility and the joint distribution ------------------------------
section("2. Feasibility and the joint distribution")

feas <- check_feasibility(m_yield3)  # triple screen, plus an LP check if lpSolve is installed
feas$feasible
joint <- fit_joint(m_yield3)          # exact maximum-entropy joint (32 combinations), with the three-way term
joint
head(combination_probs(joint), 10)
sum(combination_probs(joint)$prob[combination_probs(joint)$n_diseases >= 3])

# The sampled backend (used automatically above 20 diseases), compared with
# the exact one
joint_s <- fit_joint(m_yield3, backend = "sampled", n_samples = 20000, n_chains = 500, seed = 1)
joint_s
joint_s$diagnostics$summary
n_dis <- function(j) tapply(j$prob, rowSums(j$cells), sum)
rbind(exact = n_dis(joint), sampled = n_dis(joint_s))

# ---- 3. Adjustment -----------------------------------------------------------
section("3. Adjustment")

# yield: crude impacts, every pair has an association, no interactions or
# three-way terms: method = "auto" uses the simultaneous method. n_draws = 0
# gives point estimates only (uncertainty is in section 5).
res_yield <- deconflate(m_yield, n_draws = 0)
res_yield
# With the three-way term, "auto" switches to the global method (a note says
# why). With every pair given and no interactions, the result is the same.
res_yield3 <- deconflate(m_yield3, n_draws = 0, joint = joint)
res_yield3$notes
all.equal(res_yield3$adjusted$adjusted, res_yield$adjusted$adjusted)
# The published approximation (Rasmussen et al. 2022, eq. 16) is not a method
# of deconflate(); compare_methods() shows it beside the others, for
# comparison only
compare_methods(m_yield)$impacts

# calving interval: MAS and SCK are regression coefficients (adjusted_linear)
res_ci <- deconflate(m_ci, n_draws = 0)
res_ci$adjusted
compare_methods(m_ci)$failed   # the published approximation: crude estimates only

# welfare: interactions switch "auto" to the global method
res_welfare <- deconflate(m_welfare, n_draws = 0, joint = joint)
res_welfare
res_welfare$contributions                            # main and interaction parts

# Exact feasibility check instead of the triple screen (needs lpSolve)
if (requireNamespace("lpSolve", quietly = TRUE)) {
  deconflate(m_yield, n_draws = 0, feasibility = "lp")$diagnostics$feasibility
}

# ---- 4. Event impacts (culling) ---------------------------------------------
section("4. Event impacts")

# culling.csv mixes hazard ratios (LAM, MAS stratified by LAM and SCK), a risk
# ratio (MET), an odds ratio (SCK) and a risk difference (RP stratified by
# all others). All are mapped onto the snapshot hazard model; the overall
# risk (here 25% of animals culled in the lactation) is required, and the
# risk-based estimates must refer to the same period.
res_cull <- deconflate(m_cull, event_model = TRUE, overall_risk = 0.25, n_draws = 0, joint = joint)
res_cull
res_cull$adjusted                  # adjusted hazard ratios: each disease's own effect
res_cull$attributable$summary      # overall, disease-free and attributable risk
res_cull$attributable$by_disease   # the attributable risk split by disease (Shapley)
# event_model must match the table
try_show(deconflate(m_cull, n_draws = 0))
try_show(deconflate(m_yield, event_model = TRUE, overall_risk = 0.25))
# compare_methods(): the first-order and published (2024) approaches need
# hazard (or rate) ratios, so with this table only the snapshot model runs
compare_methods(m_cull, event_model = TRUE, overall_risk = 0.25)$failed
hr_only <- cm_model(m_cull, cm_impacts(m_cull$diseases$id, c(1.74, 1.6, 1.3, 1.4, 1.5),
                                       measure = "HR", estimand = "snapshot_crude",
                                       label = "culling"))
compare_methods(hr_only, event_model = TRUE, overall_risk = 0.25)

# ---- 5. Uncertainty -----------------------------------------------------------
section("5. Uncertainty")

# Every distribution in the tables is used: beta, pert, uniform and pert_mean
# probabilities; lognormal_ci, normal, lognormal, pert, beta and uniform
# association measures; a lognormal_ci three-way ratio; normal, pert,
# pert_mean, fixed, lognormal and uniform impacts; normal and uniform
# interactions; lognormal_ci, pert and normal event impacts. Values without a
# dist are fixed. deconflate() draws them (n_draws, 1000 by default) and
# reports 95% intervals next to the central estimates.
unc_yield <- deconflate(m_yield, n_draws = 200, seed = 1)
unc_yield
unc_yield$adjusted                 # central estimate, lower and upper
unc_yield$totals
head(unc_yield$draws$summary)      # mean, sd, Monte Carlo error and stability per quantity
unc_yield$draws$rejections         # draws that could not be used, and why
# Latin hypercube sampling: Monte Carlo errors from replicate blocks
lhs <- deconflate(m_yield, n_draws = 200, seed = 2, sampling = "lhs", lhs_replicates = 10)
lhs$draws$summary[lhs$draws$summary$quantity == "total", c("quantity", "mean", "mcse")]
# Interactions and the three-way ratio are drawn too (global method)
unc_welfare <- deconflate(m_welfare, n_draws = 100, seed = 3)
unc_welfare$totals
# Event impacts, with an uncertain overall risk
unc_cull <- deconflate(m_cull, event_model = TRUE, overall_risk = dist_beta(250, 750),
                       n_draws = 100, seed = 4)
unc_cull$attributable$summary
# Methods compared on the same draws
compare_methods(m_yield, methods = c("published", "simultaneous"), n_draws = 200, seed = 5)
# Without distributions, no draws are run
deconflate(example_supplement())$notes

# ---- 6. Reporting, gaps and values ------------------------------------------
section("6. Reporting, gaps and values")

contribution_table(res_yield)
contribution_table(unc_yield)      # with intervals
summary(res_yield)
attribute_burden(res_welfare)
contribution_table(res_cull)

# Gaps and values are computed in base R. yield: percent losses of an
# observed 10,000 kg, valued at 0.35 per kg
ry <- res_yield
disease_free <- 10000 / (1 - ry$totals$adjusted_total / 100)
gap <- disease_free - 10000
data.frame(disease = ry$contributions$disease, gap = gap * ry$contributions$share,
           value = gap * ry$contributions$share * 0.35)
# The same from every draw gives an interval for the gap
tot <- unc_yield$draws$values[, "total"]
quantile(10000 / (1 - tot / 100) - 10000, c(0.025, 0.5, 0.975))
# calving interval: absolute effects (days), so each contribution is its part
# of the gap; valued at 3 per day
data.frame(disease = res_ci$contributions$disease, days = res_ci$contributions$total,
           value = res_ci$contributions$total * 3)

# ---- 7. Sensitivity and screens -----------------------------------------------
section("7. Sensitivity and screens")

oat <- sensitivity_oat(m_yield, variation = 0.2)
head(oat, 10)
scr_a <- screen_associations(m_yield, or_values = c(0.5, 3))
scr_a
scr_i <- screen_interactions(m_yield, values = c(-1, 1))
scr_i
scr_t <- screen_three_way(m_welfare, ratios = c(0.5, 2))
scr_t
# Event impacts
head(sensitivity_oat(m_cull, event_model = TRUE, overall_risk = 0.25), 10)
# Without any association estimates: which associations would matter?
m_none <- cm_read_inputs(diseases = f("diseases.csv"), impacts = f("impacts_yield.csv"))
try_show(deconflate(m_none))
screen_associations(m_none, or_values = c(0.5, 3))

# Scenarios built with the set_* functions, adjusted like any model
strong <- set_association(m_yield, "LAM", "RP", 3)                 # an unrelated pair made associated
no_three <- set_three_way(m_welfare, "LAM", "MAS", "SCK", 1)       # no three-way association
syn <- set_interaction(m_welfare, "MAS", "SCK", 1.5)               # an extra interaction
sapply(list(base = m_yield, strong_LAM_RP = strong),
       function(m) deconflate(m, n_draws = 0)$totals$adjusted_total)
sapply(list(base = m_welfare, no_three_way = no_three, extra_interaction = syn),
       function(m) deconflate(m, n_draws = 0)$totals$adjusted_total)

# ---- 8. Thresholds: when does a conclusion change? -------------------------
section("8. Thresholds")

# rank: at what LAM-SCK odds ratio do two diseases swap places?
th_rank <- cm_threshold(m_yield, "assoc:LAM:SCK", c(0.5, 10), conclusion = "rank")
th_rank
# sign: how small a raw SCK yield impact makes its adjusted impact negative?
th_sign <- cm_threshold(m_yield, "impact:SCK", c(0, 2.5), conclusion = "sign")
th_sign$thresholds
# total: what raw MET impact gives an adjusted aggregate of 3% of yield?
th_total <- cm_threshold(m_yield, "impact:MET", c(0, 15), conclusion = "total", target = 3)
th_total$thresholds
# change: what LAM-MAS interaction lowers the welfare aggregate by 10%? (The
# raw impacts stay fixed, so a larger interaction leaves less to the main
# effects.)
th_change <- cm_threshold(m_welfare, "inter:LAM:MAS", c(-5, 20), conclusion = "change",
                          target = -0.1)
th_change$thresholds
# a three-way ratio, searched on a log scale (global method)
th_three <- cm_threshold(m_welfare, "three:LAM:MAS:SCK", c(0.2, 5), conclusion = "rank")
th_three$summary
# a disease probability: where would SCK's yield contribution overtake MET's?
th_prob <- cm_threshold(m_yield, "prob:SCK", c(0.1, 0.6), conclusion = "rank",
                        diseases = c("MET", "SCK"))
th_prob$thresholds
# Ranges that cannot be evaluated are listed in `regions`; a jump through a
# pole or a nearly singular system is reported as a discontinuity, never as a
# threshold
th_prob$regions
# event impacts: at what overall risk does the attributable risk reach 0.08?
th_risk <- cm_threshold(m_cull, "risk", c(0.1, 0.5), conclusion = "total", target = 0.08,
                        event_model = TRUE, overall_risk = 0.25)
th_risk$thresholds

# ---- 9. Simulation and non-additive losses ---------------------------------
section("9. Simulation and Shapley values")

# Simulate what studies would report from known true impacts and
# interactions, and check that the global method recovers them.
truth <- c(LAM = 8, MAS = 4, MET = 3, SCK = 1, RP = 2)
raw <- simulate_raw_impacts(m_yield3, truth,
                            interactions = cm_interactions("LAM", "MAS", 2), joint = joint,
                            units = "score points")
raw
rec <- deconflate(cm_model(m_yield3, raw, cm_interactions("LAM", "MAS", 2)), n_draws = 0,
                  joint = joint)
cbind(true = truth, recovered = rec$adjusted$adjusted)

# Coefficients of an additive regression adjusted for all other diseases
raw_adj <- simulate_raw_impacts(m_yield3, truth, estimand = "adjusted_linear",
                                adjusted_for = "all", joint = joint)
raw_adj$value  # without interactions these equal the true impacts
# A finite study of 5000 animals adds sampling error
simulate_raw_impacts(m_yield3, truth, n = 5000, joint = joint)$value

# Shapley values cell by cell: additive with a pairwise term, and
# multiplicative (proportional losses that cannot exceed 100%)
shapley_by_cell(joint, loss_additive(truth, terms = list(list(diseases = c("LAM", "MAS"),
                                                             value = 2))))
shapley_by_cell(joint, loss_multiplicative(c(LAM = 0.05, MAS = 0.03, MET = 0.06,
                                             SCK = 0.02, RP = 0.04)))

# ---- 10. Variants: unknown pairs and adjusted associations ------------------
section("10. Variants")

assoc <- utils::read.csv(f("associations.csv"), stringsAsFactors = FALSE)
# An unknown association: leave the MAS-RP row out. "auto" uses the global
# method, which fills the pair in; the result lists the odds ratio it got.
r <- which(assoc$disease1 == "MAS" & assoc$disease2 == "RP")
m_u <- read_five("impacts_yield.csv", three_way = FALSE, associations = assoc[-r, ])
res_u <- deconflate(m_u, n_draws = 0)
res_u$notes
res_u$unknown_pairs
# compare_methods() does not switch: the pairwise methods cannot run
compare_methods(m_u)$failed
# With an unknown pair, a three-way term can change additive results
screen_three_way(m_u, ratios = c(0.5, 2))
# An odds ratio of 1 states that two diseases are unrelated (LAM-RP in
# associations.csv); the measures "independent", "unknown" and "table" are
# no longer accepted, and the check says what to write instead
old <- rbind(assoc[, c("disease1", "disease2", "value", "measure")],
             data.frame(disease1 = "MAS", disease2 = "RP", value = NA, measure = "unknown"))
old <- old[-r, ]
old$measure[old$disease1 == "LAM" & old$disease2 == "RP"] <- "independent"
cm_check_inputs(diseases = f("diseases.csv"), associations = old, impacts = f("impacts_yield.csv"))

# A covariate-adjusted risk ratio (MET-SCK): rejected unless used as a
# marginal approximation
a_adj <- assoc
r2 <- which(a_adj$disease1 == "MET" & a_adj$disease2 == "SCK")
a_adj$adjusted[r2] <- TRUE
a_adj$adjusted_for[r2] <- "LAM"
cm_check_inputs(diseases = f("diseases.csv"), associations = a_adj, impacts = f("impacts_yield.csv"))
m_adj <- cm_read_inputs(diseases = f("diseases.csv"), associations = a_adj,
                        impacts = f("impacts_yield.csv"), adjusted_associations = "use_as_marginal")
deconflate(m_adj, n_draws = 0)$adjusted

# ---- 11. Plots ---------------------------------------------------------------
section("11. Plots")

try_plot(plot(res_yield))
try_plot(plot(unc_yield))
try_plot(plot_burden(res_yield))
try_plot(plot(res_cull))
try_plot(plot_burden(unc_cull))
try_plot(plot(oat))
try_plot(plot(scr_a))
try_plot(plot(th_rank))
try_plot(plot(th_change))

if (save_plots) {
  grDevices::dev.off()
  cat("\nPlots written to", plot_file, "\n")
}
if (!is.null(out_dir)) {
  cat("Printed output written to", text_file, "\n")
  sink()
}
