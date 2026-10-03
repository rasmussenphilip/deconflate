# deconflate: a tour of every feature with the five-disease example
#
# The CSV files in this folder describe an illustrative dairy population
# with five diseases (LAM lameness, MAS mastitis, MET metritis, SCK
# subclinical ketosis, RP retained placenta) and three impact analyses:
# yield (% of yield), calving_interval (days) and welfare (score points).
# The values are made up to exercise the package; they are not estimates.
# See README.md for what each file contains.
#
# Step through the script section by section, or run it all and see every
# result with
#   source(system.file("extdata", "five_diseases", "run_all_features.R",
#                      package = "deconflate"), echo = TRUE, max.deparse.length = Inf)
# (plain source() runs it silently). Monte Carlo sections use small numbers
# of draws so that the script runs in a few minutes; use more draws for real
# analyses.

library(deconflate)

five <- system.file("extdata", "five_diseases", package = "deconflate")
if (!nzchar(five)) five <- "."  # running from a copy of this folder

# Draw plots to a PDF file instead of the screen? (Set TRUE when sourcing
# the script non-interactively.)
save_plots <- !interactive()
if (save_plots) {
  plot_file <- file.path(tempdir(), "deconflate_five_diseases.pdf")
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

# No problems are expected (a note is printed if any).
print(cm_check_inputs(dir = five))
inp <- cm_read_inputs(dir = five)
inp
names(inp$analyses$models)      # calving_interval, welfare, yield
inp$population$diseases          # MAS is an incidence rate, converted to a probability
pair_tables(inp$population)      # the 2x2 table of every pair, from each measure
inp$population$three_way         # a three-way scenario (global method only)
inp$sampler                      # a batch sampler: one per analysis, shared population draws
attr(inp$sampler$samplers$yield, "specs")  # the distributions given in the tables

m_yield <- inp$analyses$models$yield
m_ci <- inp$analyses$models$calving_interval
m_welfare <- inp$analyses$models$welfare

# ---- 2. Feasibility and the joint distribution ------------------------------
section("2. Feasibility and the joint distribution")

feas <- check_feasibility(inp$population)  # triple screen, plus an LP check if lpSolve is installed
feas$feasible
joint <- fit_joint(inp$population)          # exact maximum-entropy joint (32 combinations)
joint
head(combination_probs(joint), 10)
sum(combination_probs(joint)$prob[combination_probs(joint)$n_diseases >= 3])

# The sampled backend (needed above about 20 diseases), compared with the exact one
joint_s <- fit_joint(inp$population, backend = "sampled", n_samples = 20000,
                     n_chains = 500, seed = 1)
joint_s
joint_s$diagnostics$summary
n_dis <- function(j) tapply(j$prob, rowSums(j$cells), sum)
rbind(exact = n_dis(joint), sampled = n_dis(joint_s))

# ---- 3. Adjustment: every method, every analysis ----------------------------
section("3. Adjustment")

# yield: crude impacts, no interactions: all three methods work
res_yield <- lapply(c(simultaneous = "simultaneous", published = "published", global = "global"),
                    function(mt) deconflate(m_yield, method = mt, joint = joint))
res_yield$simultaneous
adj_yield <- sapply(res_yield, function(r) r$adjusted$adjusted)
rownames(adj_yield) <- m_yield$diseases$id
adj_yield
# With all pairs constrained, the three-way term does not change additive
# results without interactions: global equals simultaneous here.
all.equal(res_yield$global$adjusted$adjusted, res_yield$simultaneous$adjusted$adjusted)

# calving_interval: MAS and SCK are regression coefficients (adjusted_linear)
res_ci <- deconflate(m_ci)
res_ci$adjusted
try_show(deconflate(m_ci, method = "published"))   # published: crude estimates only

# welfare: interactions need the global method
res_welfare <- deconflate(m_welfare, method = "global", joint = joint)
res_welfare
res_welfare$contributions                            # main and interaction parts
try_show(deconflate(m_welfare))                      # simultaneous: unsupported with interactions

# All analyses at once (global; one joint distribution serves all)
res_all <- deconflate(inp$analyses, method = "global")
sapply(res_all, function(r) r$totals$adjusted_total)

# Exact feasibility check instead of the triple screen (needs lpSolve)
if (requireNamespace("lpSolve", quietly = TRUE)) {
  deconflate(m_yield, feasibility = "lp")$diagnostics$feasibility
}

# ---- 4. Reporting and valuation ---------------------------------------------
section("4. Reporting and valuation")

val_yield <- list(observed = 10000, direction = "decrease", effect = "percent", unit_value = 0.35)
val_ci <- list(observed = 400, direction = "increase", effect = "absolute", unit_value = 3)

contribution_table(res_yield$simultaneous, val_yield)
summary(res_yield$simultaneous, valuation = val_yield)
attribute_burden(res_welfare)

gap <- productivity_gap(res_yield$simultaneous, observed = 10000, direction = "decrease",
                        effect = "percent")
gap$summary
gap$attribution
value_losses(gap, unit_value = 0.35, additional = 2500)

gap_ci <- productivity_gap(res_ci, observed = 400, direction = "increase", effect = "absolute")
value_losses(gap_ci, unit_value = 3)

# ---- 5. Comparing methods ---------------------------------------------------
section("5. Comparing methods")

# Each analysis with every method; unsupported or undefined methods are
# reported as failed, not dropped silently.
cmp <- compare_methods(inp$analyses,
                       valuation = list(yield = val_yield, calving_interval = val_ci))
cmp$totals
cmp$failed

# ---- 6. Hazard ratios (separate snapshot model) -----------------------------
section("6. Hazard ratios")

inp$hr_model
hr_snap <- deconflate_hr(inp$hr_model, joint = joint)          # snapshot (default)
hr_snap$adjusted
deconflate_hr(inp$hr_model, method = "first_order")$adjusted
try_show(deconflate_hr(inp$hr_model, method = "published"))     # crude hazard ratios only
compare_methods(inp$hr_model, overall_risk = 0.25)
ar <- attributable_risk(hr_snap, overall_risk = 0.25, unit_value = 1300)
ar$summary
ar$by_disease

# ---- 7. Monte Carlo: all analyses on shared population draws ----------------
section("7. Monte Carlo, batch")

# Every distribution in the tables is used: beta, pert, uniform and pert_mean
# probabilities; lognormal_ci, normal, lognormal, pert, beta and uniform
# association measures; a lognormal_ci three-way ratio; normal, pert,
# pert_mean, fixed, lognormal and uniform impacts; normal and uniform
# interactions. Values without a dist are fixed.
mc_all <- cm_monte_carlo(inp$sampler, 200, method = "global", seed = 1)
mc_all
sm_all <- summary(mc_all)
sm_all[, c("analysis", "method", "disease", "mean", "q0.025", "q0.975", "stability")]

# ---- 8. Monte Carlo: one analysis, with variance reduction ------------------
section("8. Monte Carlo, single analysis")

s_yield <- inp$sampler$samplers$yield
s_yield

# Two methods on identical draws
mc <- cm_monte_carlo(s_yield, 400, method = c("published", "simultaneous"), seed = 2)
mc
summary(mc, what = "rejections")
summary(mc)
summary(mc, what = "total")
compare_methods(mc)
cm_diagnose(mc)

# Latin hypercube sampling: standard errors from replicate blocks
mc_lhs <- cm_monte_carlo(s_yield, 400, sampling = "lhs", lhs_replicates = 10, seed = 3)
summary(mc_lhs, diagnose = FALSE)[, c("disease", "mean", "mcse", "stability")]

# Importance sampling with a suggested defensive-mixture proposal (a
# demonstration: the proposal targets the input most related to SCK's most
# extreme adjusted impacts; compare the standard errors with the plain run)
prop <- cm_suggest_proposal(mc, disease = "SCK", method = "simultaneous")
mc_is <- cm_monte_carlo(s_yield, 400, proposal = prop, seed = 4)
mc_is$ess
summary(mc_is, diagnose = FALSE)[, c("disease", "mean", "mcse", "stability")]
# A point mass (SCK's yield impact is fixed(2.5)) cannot be importance-sampled
try_show(cm_monte_carlo(s_yield, 10, proposal = list("impact:SCK" = dist_normal(2.5, 1))))

# Scenarios by reweighting the existing draws (no re-run)
alt <- cm_scenario(mc, list("impact:LAM" = dist_normal(5.5, 0.6)))
alt$ess
summary(alt, diagnose = FALSE)[, c("method", "disease", "mean")]
lr <- function(p) {
  x <- p[["assoc:LAM:SCK"]]
  stats::dnorm(x, 2.2, 0.3, log = TRUE) - stats::dnorm(x, 2.0, 0.3, log = TRUE)
}
summary(cm_reweight(mc, lr), diagnose = FALSE)[, c("method", "disease", "mean")]

# Uncertain observed yield and price, carried through every draw
g <- cm_mc_gap(mc, observed = dist_normal(10000, 500), direction = "decrease",
               effect = "percent", unit_value = dist_uniform(0.30, 0.40), seed = 5)
g$summary

# A sampler built in R, with a three-way ratio among the uncertain inputs
s_three <- cm_sampler(m_welfare,
                      impacts = cm_dist_table(data.frame(key = c("LAM", "MAS"),
                                                         dist = "normal",
                                                         p1 = c(10, 6), p2 = c(2, 1))),
                      three_way = list("SCK:MAS:LAM" = dist_lognormal_ci(1.5, 0.8, 2.8)))
s_three
mc_three <- cm_monte_carlo(s_three, 100, method = "global", seed = 6)
summary(mc_three, diagnose = FALSE)[, c("disease", "mean", "q0.025", "q0.975")]

# ---- 9. Sensitivity and scenario screens ------------------------------------
section("9. Sensitivity and screens")

oat <- sensitivity_oat(m_yield, variation = 0.2, valuation = val_yield)
head(oat, 10)
scr_a <- screen_associations(m_yield, or_values = c(0.5, 3))
scr_a
scr_i <- screen_interactions(m_yield, values = c(-1, 1))
scr_i
scr_t <- screen_three_way(m_welfare, ratios = c(0.5, 2))
scr_t

# Named scenarios, built with the set_* functions (models with a three-way
# term or interactions are adjusted with the global method)
strong <- set_association(m_yield, "LAM", "RP", 3)               # an independent pair made associated
no_three <- set_three_way(m_welfare, "LAM", "MAS", "SCK", 1)     # no three-way association
syn <- set_interaction(m_welfare, "MAS", "SCK", 1.5)             # an extra interaction
compare_scenarios(base = m_yield, strong_LAM_RP = strong)$totals
compare_scenarios(base = m_welfare, no_three_way = no_three, extra_interaction = syn,
                  method = "global")$totals

# ---- 10. Thresholds: when does a conclusion change? -------------------------
section("10. Thresholds")

# rank: at what LAM-SCK odds ratio do two diseases swap places?
th_rank <- cm_threshold(m_yield, "assoc:LAM:SCK", c(0.5, 10), conclusion = "rank",
                        method = "simultaneous")
th_rank
# sign: how small a raw SCK yield impact makes its adjusted impact negative?
th_sign <- cm_threshold(m_yield, "impact:SCK", c(0, 2.5), conclusion = "sign",
                        method = "simultaneous")
th_sign$thresholds
# total: what raw MET impact gives an adjusted aggregate of 3% of yield?
th_total <- cm_threshold(m_yield, "impact:MET", c(0, 15), conclusion = "total", target = 3,
                         method = "simultaneous")
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
                        diseases = c("MET", "SCK"), method = "simultaneous")
th_prob$thresholds
# the published method has a pole where its denominator is zero (for SCK at
# a raw impact of about -1.43): it is reported as a discontinuity, never as
# a threshold
th_pub <- cm_threshold(m_yield, "impact:SCK", c(-4, 2.5), conclusion = "sign",
                       diseases = "SCK", method = "published")
th_pub$thresholds
th_pub$regions

# ---- 11. Simulation and non-additive losses ---------------------------------
section("11. Simulation and Shapley values")

# Simulate what studies would report from known true impacts and
# interactions, and check that the global method recovers them.
truth <- c(LAM = 8, MAS = 4, MET = 3, SCK = 1, RP = 2)
raw <- simulate_raw_impacts(inp$population, truth,
                            interactions = cm_interactions("LAM", "MAS", 2), joint = joint,
                            units = "score points")
raw
rec <- deconflate(cm_model(inp$population, raw, cm_interactions("LAM", "MAS", 2)),
                  method = "global", joint = joint)
cbind(true = truth, recovered = rec$adjusted$adjusted)

# Coefficients of an additive regression adjusted for all other diseases
raw_adj <- simulate_raw_impacts(inp$population, truth, estimand = "adjusted_linear",
                                adjusted_for = "all", joint = joint)
raw_adj$value  # without interactions these equal the true impacts
# A finite study of 5000 animals adds sampling error
simulate_raw_impacts(inp$population, truth, n = 5000, joint = joint)$value

# Shapley values cell by cell: additive with a pairwise term, and
# multiplicative (proportional losses that cannot exceed 100%)
shapley_by_cell(joint, loss_additive(truth, terms = list(list(diseases = c("LAM", "MAS"),
                                                             value = 2))))
shapley_by_cell(joint, loss_multiplicative(c(LAM = 0.05, MAS = 0.03, MET = 0.06,
                                             SCK = 0.02, RP = 0.04)))

# ---- 12. Variants: unknown pairs and adjusted associations ------------------
section("12. Variants")

# An unknown association (MAS-RP): only the global method can use it.
a_unknown <- inp$tables$associations
r <- which(a_unknown$disease1 == "MAS" & a_unknown$disease2 == "RP")
a_unknown$measure[r] <- "unknown"
a_unknown$value[r] <- NA
inp_u <- cm_read_inputs(dir = five, associations = a_unknown)
try_show(deconflate(inp_u$analyses$models$yield))
deconflate(inp_u$analyses$models$yield, method = "global")$adjusted
# The same with every unlisted pair unknown
a_drop <- inp$tables$associations[-r, ]
inp_u2 <- cm_read_inputs(dir = five, associations = a_drop, missing_associations = "unknown")
pair_tables(inp_u2$population)[, c("disease1", "disease2", "status")]
# With unknown pairs, a three-way term can change additive results
screen_three_way(inp_u2$analyses$models$yield, ratios = c(0.5, 2))

# A covariate-adjusted risk ratio (MET-SCK): rejected unless used as a
# marginal approximation
a_adj <- inp$tables$associations
r2 <- which(a_adj$disease1 == "MET" & a_adj$disease2 == "SCK")
a_adj$adjusted[r2] <- "TRUE"
a_adj$adjusted_for[r2] <- "LAM"
cm_check_inputs(dir = five, associations = a_adj)
inp_adj <- cm_read_inputs(dir = five, associations = a_adj,
                          adjusted_associations = "use_as_marginal")
deconflate(inp_adj$analyses$models$yield)$adjusted

# ---- 13. Plots ---------------------------------------------------------------
section("13. Plots")

try_plot(plot(res_yield$simultaneous))
try_plot(plot(res_all))
try_plot(plot_burden(res_yield$simultaneous, valuation = val_yield))
try_plot(plot(mc))
try_plot(plot(oat))
try_plot(plot(scr_a))
try_plot(plot(th_rank))
try_plot(plot(th_change))

if (save_plots) {
  grDevices::dev.off()
  cat("\nPlots written to", plot_file, "\n")
}
