# deconflate: a tour of every feature with the five-disease example
#
# The CSV files in this folder describe an illustrative dairy population
# with five diseases (LAM lameness, MAS mastitis, MET metritis, SCK
# subclinical ketosis, RP retained placenta) and four impact tables: yield
# (% of yield), calving interval (days, with several rows for one disease and
# a subgroup row), welfare (score points, with interactions) and culling
# (event impacts: hazard ratios, a risk ratio, an odds ratio, a risk
# difference and an interaction). The values are made up to exercise the
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

# Three tables, each its own argument, with any file name. One impact table
# per model: the disease and association tables are shared.
read_five <- function(impacts, associations = f("associations.csv")) {
  cm_read_inputs(diseases = f("diseases.csv"), associations = associations, impacts = f(impacts))
}

# No problems are expected (a note is printed if any).
print(cm_check_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                      impacts = f("impacts_yield.csv")))
m_yield <- read_five("impacts_yield.csv")
m_yield
m_yield$diseases                # MAS is an incidence rate, converted to a probability
m_yield$associations            # every kind of association row (see the note column)
names(m_yield$distributions)    # the distributions given in the tables

m_ci <- read_five("impacts_calving_interval.csv")
m_welfare <- read_five("impacts_welfare.csv")
m_cull <- read_five("culling.csv")
m_ci$impacts                    # two MAS rows, an adjusted row, a subgroup row
m_welfare$impacts               # combination rows: interactions
m_cull$impacts                  # event impacts: a measure column

# The yield model with only the ten crude pair rows (the layout of earlier
# versions): the pairwise tables then give the answer directly
assoc <- utils::read.csv(f("associations.csv"), stringsAsFactors = FALSE)
m_pairs <- read_five("impacts_yield.csv", associations = assoc[-c(2, 12), ])

# ---- 2. Feasibility and the joint distribution ------------------------------
section("2. Feasibility and the joint distribution")

feas <- check_feasibility(m_yield)  # triple screen, plus an LP check if lpSolve is installed
feas$feasible
joint <- fit_joint(m_yield)          # exact maximum-entropy joint (32 combinations), matching every row
joint
cmb <- data.frame(joint$cells, prob = joint$prob)
head(cmb[order(-cmb$prob), ], 10)    # the most common combinations
sum(joint$prob[rowSums(joint$cells) >= 3])  # share of cows with three or more diseases
# The OR_ratio row: the LAM-MAS odds ratio among cows with and without SCK
or_in <- function(j, among) {
  w <- j$prob * (j$cells[, "SCK"] == among)
  c11 <- sum(w[j$cells[, "LAM"] == 1 & j$cells[, "MAS"] == 1])
  c10 <- sum(w[j$cells[, "LAM"] == 1 & j$cells[, "MAS"] == 0])
  c01 <- sum(w[j$cells[, "LAM"] == 0 & j$cells[, "MAS"] == 1])
  c00 <- sum(w[j$cells[, "LAM"] == 0 & j$cells[, "MAS"] == 0])
  c11 * c00 / (c10 * c01)
}
or_in(joint, 1) / or_in(joint, 0)    # 1.5, as the row says

# The sampled backend (used automatically above 20 diseases) handles crude
# pair rows only; compared with the exact one on the pair-only model
joint_p <- fit_joint(m_pairs)
joint_s <- fit_joint(m_pairs, backend = "sampled", n_samples = 20000, n_chains = 500, seed = 1)
joint_s
joint_s$diagnostics$summary
n_dis <- function(j) tapply(j$prob, rowSums(j$cells), sum)
rbind(exact = n_dis(joint_p), sampled = n_dis(joint_s))
try_show(fit_joint(m_yield, backend = "sampled", n_samples = 2000, n_chains = 100, seed = 1))

# ---- 3. Adjustment -----------------------------------------------------------
section("3. Adjustment")

# yield: crude impacts, one row per disease. The association rows include a
# pooled pair and an OR_ratio row, so the full model is fitted.
# n_draws = 0 gives point estimates only (uncertainty is in section 5).
res_yield <- deconflate(m_yield, n_draws = 0)
res_yield
res_yield$diagnostics$method
res_yield$heterogeneity          # the two LAM-MAS studies, pooled (random effects)
res_yield$estimates              # every row with its fitted value
deconflate(m_yield, n_draws = 0, pool = "fixed")$heterogeneity
# Only crude pair rows: the pairwise tables give the answer directly
res_pairs <- deconflate(m_pairs, n_draws = 0)
res_pairs$diagnostics$method
# The OR_ratio row keeps the pairwise tables fixed: with the same pair rows,
# it changes additive results without interactions not at all
m_pairs_or <- read_five("impacts_yield.csv", associations = assoc[-2, ])
all.equal(deconflate(m_pairs_or, n_draws = 0)$adjusted$adjusted, res_pairs$adjusted$adjusted)

# calving interval: two MAS rows (one adjusted for SCK, one crude) fitted
# together by weighted least squares, SCK adjusted for all (exact), and RP
# among cows with no other disease
res_ci <- deconflate(m_ci, n_draws = 0)
res_ci$adjusted
res_ci$estimates[res_ci$estimates$table == "impacts", ]
res_ci$heterogeneity

# welfare: interactions (LAM;MAS and MET;RP adjusted for all, LAM;SCK crude)
res_welfare <- deconflate(m_welfare, n_draws = 0, joint = joint)
res_welfare
res_welfare$contributions                            # main and interaction parts
res_welfare$conflation$A                             # what each row contains of each term

# Diagnosing the inputs: with a raw SCK yield impact of 1.0 instead of 2.5,
# SCK's adjusted impact changes sign (a warning). diagnose = TRUE shows
# which inputs drive it: the change per SD of each uncertain input, the shift
# that alone would remove the flag, and the value of each exact input at
# which it disappears. It reruns the model about twice per uncertain input,
# so it is off by default.
m_low <- m_yield
m_low$impacts$value[m_low$impacts$disease == "SCK"] <- 1.0
res_low <- deconflate(m_low, n_draws = 0, diagnose = TRUE)
res_low$problems
res_low$problems$flagged                             # every disease: raw, adjusted, flagged?
head(res_low$problems$inputs)                        # every uncertain input, per flagged disease
# Rows that disagree with other rows: the crude MAS calving-interval row at
# 20 days instead of 9
m_ci20 <- m_ci
m_ci20$impacts$value[3] <- 20
deconflate(m_ci20, n_draws = 0, diagnose = TRUE)$problems

# ---- 4. Event impacts (culling) ---------------------------------------------
section("4. Event impacts")

# culling.csv mixes hazard ratios (LAM, MAS stratified by LAM and SCK), a risk
# ratio (MET), an odds ratio (SCK), a risk difference (RP stratified by all
# others) and a crude hazard ratio for cows with both LAM and MAS (an
# interaction). All are mapped onto the snapshot hazard model; the overall
# risk (here 25% of animals culled in the lactation) is required, and the
# risk-based estimates must refer to the same period.
res_cull <- deconflate(m_cull, event_model = TRUE, overall_risk = 0.25, n_draws = 0, joint = joint)
res_cull
res_cull$adjusted                  # adjusted hazard ratios, and each term alone vs
                                   # disease-free on the risk scale (risk, RR, RD, OR)
res_cull$attributable$summary      # overall, disease-free, attributable and raw
                                   # attributable risk (the diseases unrelated)
res_cull$attributable$by_disease   # split by disease (Shapley): own part and share
                                   # of the LAM;MAS interaction
round(res_cull$conflation$A, 3)    # how each row's value moves with each term
# The risk of disease-free animals instead of the overall risk: the model
# then gives the overall risk
res_base <- deconflate(m_cull, event_model = TRUE, baseline_risk = 0.18, n_draws = 0, joint = joint)
res_base$attributable$summary
# Diagnosing the inputs for event impacts: with an SCK odds ratio of 1.2, the
# adjusted odds ratio of a cow with only SCK against a disease-free cow falls
# below 1, and LAM's hazard ratio shrinks by more than 90% (on the log scale);
# results are compared on each raw estimate's own measure
m_cull_low <- m_cull
m_cull_low$impacts$value[m_cull_low$impacts$disease == "SCK"] <- 1.2
deconflate(m_cull_low, event_model = TRUE, overall_risk = 0.25, n_draws = 0, diagnose = TRUE)$problems
# event_model must match the table
try_show(deconflate(m_cull, n_draws = 0))
try_show(deconflate(m_yield, event_model = TRUE, overall_risk = 0.25))

# ---- 5. Uncertainty -----------------------------------------------------------
section("5. Uncertainty")

# Every distribution in the tables is used: beta, pert, uniform and pert_mean
# probabilities; lognormal_ci, normal, lognormal, pert, beta and uniform
# association rows; normal, pert, pert_mean, fixed, lognormal and uniform
# impact rows; lognormal_ci, pert and normal event impacts. Values without a
# dist are fixed. deconflate() draws them (n_draws, 1000 by default) and
# reports 95% intervals next to the central estimates.
unc_yield <- deconflate(m_yield, n_draws = 200, seed = 1)
unc_yield
unc_yield$adjusted                 # central estimate, lower and upper
unc_yield$totals
head(unc_yield$draws$summary)      # mean, sd, Monte Carlo error and stability per quantity
unc_yield$draws$rejections         # draws that could not be used, and why
# Latin hypercube sampling: Monte Carlo errors from replicate blocks
lhs <- deconflate(m_pairs, n_draws = 200, seed = 2, sampling = "lhs", lhs_replicates = 10)
lhs$draws$summary[lhs$draws$summary$quantity == "total", c("quantity", "mean", "mcse")]
# Interaction rows are drawn too
unc_welfare <- deconflate(m_welfare, n_draws = 100, seed = 3)
unc_welfare$totals
# Several impact rows: each draw is fitted by weighted least squares
unc_ci <- deconflate(m_ci, n_draws = 100, seed = 6)
unc_ci$adjusted
# Event impacts, with an uncertain overall risk
unc_cull <- deconflate(m_cull, event_model = TRUE, overall_risk = dist_beta(250, 750),
                       n_draws = 100, seed = 4)
unc_cull$attributable$summary
# Without distributions, no draws are run
deconflate(example_supplement())$notes

# ---- 6. Reporting, gaps and values ------------------------------------------
section("6. Reporting, gaps and values")

contribution_table(res_yield)
contribution_table(unc_yield)      # with intervals
summary(res_yield)
res_welfare$contributions          # own part and share of the interactions
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
# Event impacts (the aggregate is the attributable risk), including an
# interaction between two diseases as a hazard ratio
head(sensitivity_oat(m_cull, event_model = TRUE, overall_risk = 0.25), 10)
screen_interactions(m_cull, values = c(0.5, 2), pairs = c("MET:SCK", "SCK:RP"), event_model = TRUE,
                    overall_risk = 0.25)
# Without any association estimates: which associations would matter?
m_none <- cm_read_inputs(diseases = f("diseases.csv"), impacts = f("impacts_yield.csv"))
try_show(deconflate(m_none))
screen_associations(m_none, or_values = c(0.5, 3))

# Scenarios built with the set_* functions, adjusted like any model
strong <- set_association(m_yield, "LAM", "RP", 3)                 # an unrelated pair made associated
no_three <- set_three_way(m_welfare, "LAM", "MAS", "SCK", 1)       # replaces the OR_ratio row
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
# the OR_ratio row, searched on a log scale
th_three <- cm_threshold(m_welfare, "three:LAM:MAS:SCK", c(0.2, 5), conclusion = "rank")
th_three$summary
# any row by its number: the second LAM-MAS study
th_row <- cm_threshold(m_yield, "associations:2", c(0.5, 10), conclusion = "change", target = -0.05)
th_row$thresholds
# a disease probability: where would SCK's yield contribution overtake MET's?
th_prob <- cm_threshold(m_yield, "prob:SCK", c(0.1, 0.6), conclusion = "rank",
                        diseases = c("MET", "SCK"))
th_prob$thresholds
# Ranges that cannot be evaluated are listed in `regions`; a jump through a
# pole or a singular system is reported as a discontinuity, never as a
# threshold
th_prob$regions
# event impacts: at what overall risk does the attributable risk reach 0.08?
th_risk <- cm_threshold(m_cull, "risk", c(0.1, 0.5), conclusion = "total", target = 0.08,
                        event_model = TRUE, overall_risk = 0.25)
th_risk$thresholds

# ---- 9. Simulation ------------------------------------------------------------
section("9. Simulation")

# Simulate what studies would report from known true impacts and
# interactions, and check that the model recovers them. By default: a crude
# row for each disease and each combination.
truth <- c(LAM = 8, MAS = 4, MET = 3, SCK = 1, RP = 2)
raw <- simulate_raw_impacts(m_yield, truth, interactions = c("LAM;MAS" = 2), joint = joint,
                            units = "score points")
raw
rec <- deconflate(cm_model(m_yield, raw), n_draws = 0, joint = joint)
cbind(true = c(truth, "LAM;MAS" = 2), recovered = rec$adjusted$adjusted)
# Rows of any kind: MAS adjusted for SCK, RP among no others, the LAM;MAS
# combination (crude: both against neither), and LAM within MAS. The
# interaction is a term of the model because a row names the combination;
# the seven rows over-determine the six terms, and agree exactly.
rows <- data.frame(disease = c("LAM", "MAS", "MET", "SCK", "RP", "LAM;MAS", "LAM"),
                   among = c(NA, NA, NA, NA, "no others", NA, "MAS"),
                   adjusted_for = c(NA, "SCK", NA, NA, NA, NA, NA))
raw_rows <- simulate_raw_impacts(m_yield, truth, interactions = c("LAM;MAS" = 2), rows = rows,
                                 joint = joint)
raw_rows
deconflate(cm_model(m_yield, raw_rows), n_draws = 0, joint = joint)$adjusted

# Coefficients of an additive regression adjusted for all other diseases
raw_adj <- simulate_raw_impacts(m_yield, truth, adjusted_for = "all", joint = joint)
raw_adj$value  # without interactions these equal the true impacts
# A finite study of 5000 animals adds sampling error
simulate_raw_impacts(m_yield, truth, n = 5000, joint = joint)$value
# Event impacts: what an HR, RR, OR and RD row would report given true hazard
# ratios, and the model recovering them
truth_hr <- c(LAM = 1.3, MAS = 1.5, MET = 1.2, SCK = 1.1, RP = 1.6)
raw_cull <- simulate_raw_impacts(m_cull, truth_hr, interactions = c("LAM;MAS" = 1.4), event_model = TRUE,
                                 measure = c("HR", "RR", "OR", "RD", "HR", "HR"), overall_risk = 0.25,
                                 joint = joint)
raw_cull
deconflate(cm_model(m_cull, raw_cull), event_model = TRUE, overall_risk = 0.25, n_draws = 0,
           joint = joint)$adjusted[, c("disease", "measure", "raw", "adjusted")]

# ---- 10. Variants: unknown pairs, adjusted associations, old layouts --------
section("10. Variants")

# An unknown association: leave the MAS-RP row out. The fit fills the pair
# in; the result lists the odds ratio it got.
r <- which(assoc$disease == "MAS" & assoc$associated_with == "RP")
m_u <- read_five("impacts_yield.csv", associations = assoc[-r, ])
res_u <- deconflate(m_u, n_draws = 0)
res_u$unknown_pairs
# With an unknown pair, an OR_ratio row can change additive results
screen_three_way(m_u, ratios = c(0.5, 2), triples = list(c("MAS", "RP", "MET")))

# An odds ratio adjusted for other diseases of the model (MET-SCK given LAM,
# from a logistic regression): matched as the regression coefficient. Ids
# that are not diseases (parity) are noted and ignored.
a_adj <- assoc
r2 <- which(a_adj$disease == "MET" & a_adj$associated_with == "SCK")
a_adj$measure[r2] <- "OR"
a_adj$value[r2] <- 1.6
a_adj[r2, c("dist", "p1", "p2", "p3")] <- NA
a_adj$adjusted_for[r2] <- "LAM; parity"
m_adj <- read_five("impacts_yield.csv", associations = a_adj)
res_adj <- deconflate(m_adj, n_draws = 0)
res_adj$notes
res_adj$adjusted

# Files of earlier versions: the check says what to write instead
old <- data.frame(disease1 = c("LAM", "MAS"), disease2 = c("MAS", "RP"), value = c(1.8, NA),
                  measure = c("OR", "unknown"))
cm_check_inputs(diseases = f("diseases.csv"), associations = old, impacts = f("impacts_yield.csv"))
old2 <- assoc[1:3, c("disease", "associated_with", "measure", "value")]
old2$measure[2] <- "cond_prob"
cm_check_inputs(diseases = f("diseases.csv"), associations = old2, impacts = f("impacts_yield.csv"))

# ---- 11. Plots ---------------------------------------------------------------
section("11. Plots")

try_plot(plot(res_yield))
try_plot(plot(unc_yield))
try_plot(plot(res_cull))
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
