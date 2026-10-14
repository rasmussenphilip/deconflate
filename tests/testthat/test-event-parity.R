# Event impacts matched to additive impacts (v0.6): raw totals with the
# diseases unrelated, the main/interaction split of contributions, the local
# row equations, results on the risk scale, a baseline risk instead of the
# overall risk, simulated raw estimates, and interactions in the sensitivity
# tools. Reference values: inst/validation/reference_v06_event_parity.py (an
# independent implementation on reference_v041.py).

cull_five <- function(...) deconflate(read_five("culling.csv"), event_model = TRUE, n_draws = 0, ...)

test_that("raw totals: the same model with the diseases unrelated", {
  expect_equal(deconflate(read_five("impacts_yield.csv"), n_draws = 0)$totals$raw_total, 3.7422998718,
               tolerance = 1e-7)
  # Several rows for MAS, a subgroup row (RP among no others), SCK adjusted
  # for all: the raw total is still defined.
  expect_equal(deconflate(read_five("impacts_calving_interval.csv"), n_draws = 0)$totals$raw_total, 8.5943353953,
               tolerance = 1e-7)
  # Combination rows: their terms are fitted with the diseases unrelated.
  expect_equal(deconflate(read_five("impacts_welfare.csv"), n_draws = 0)$totals$raw_total, 6.0346247863,
               tolerance = 1e-7)
  r <- cull_five(overall_risk = 0.25)
  s <- r$attributable$summary
  expect_equal(s$attributable, 0.0541469811, tolerance = 1e-7)
  expect_equal(s$raw_attributable, 0.0805814757, tolerance = 1e-7)
  sm <- summary(r)$totals
  expect_equal(sm$reduction, 1 - 0.0541469811 / 0.0805814757, tolerance = 1e-8)
  expect_true(any(grepl("Raw attributable risk (diseases unrelated): 0.08058", capture.output(print(r)), fixed = TRUE)))
})

test_that("contributions are split into own parts and interaction shares", {
  w <- deconflate(read_five("impacts_welfare.csv"), n_draws = 0)$contributions
  expect_equal(w$main, c(2.1860744234, 0.9931824543, 0.3197758625, 0.3799454881, 0.1290929305), tolerance = 1e-7)
  expect_equal(w$total, c(2.2400759141, 1.0837970506, 0.3357170535, 0.3433323825, 0.1450341215), tolerance = 1e-7)
  by <- cull_five(overall_risk = 0.25)$attributable$by_disease
  expect_equal(names(by), c("disease", "main", "interaction", "attributable", "share"))
  expect_equal(by$attributable, c(0.0137316223, 0.0218514738, 0.0015789582, 0.0114739881, 0.0055109386),
               tolerance = 1e-7)
  expect_equal(by$main, c(0.0028275807, 0.0109474323, 0.0015402824, 0.0111226092, 0.0054178744), tolerance = 1e-7)
  expect_equal(by$interaction, by$attributable - by$main)
  # The LAM;MAS interaction is shared equally between LAM and MAS.
  expect_equal(by$interaction[1], 0.0109040415, tolerance = 1e-7)
  expect_equal(by$interaction[2], 0.0109040415, tolerance = 1e-7)
  # Without interaction terms, the own part is the whole contribution.
  h <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.2, n_draws = 0)$attributable$by_disease
  expect_equal(h$main, h$attributable)
  expect_equal(h$interaction, c(0, 0, 0))
  ct <- contribution_table(cull_five(overall_risk = 0.25))
  expect_true(all(c("main", "interaction", "attributable", "risk_alone", "rr_alone", "rd_alone", "change") %in% names(ct)))
})

test_that("local row equations of event impacts", {
  A <- cull_five(overall_risk = 0.25)$conflation$A
  expect_equal(dimnames(A), list(paste0("impacts:", 1:6), c("LAM", "MAS", "MET", "SCK", "RP", "LAM;MAS")))
  expect_equal(unname(A[1, ]), c(1, 0.3574732053, 0.0502701553, 0.2204200934, 0.0090040591, 0.6275884151),
               tolerance = 1e-6)
  expect_equal(unname(A[3, ]), c(0.0604553528, 0.0692425509, 0.7961978146, 0.0876879499, 0.1588654324, 0.0503107931),
               tolerance = 1e-6)
  expect_equal(unname(A[5, ]), c(-0.0075432019, -0.0038734884, 0.0106076282, -0.000273039, 0.2671313274,
                                 -0.0076138512), tolerance = 1e-6)
  # Draws do not keep them (they are a detail of the central estimate).
  expect_null(adjust_event(read_five("culling.csv"), 0.25, fit_joint(read_five("culling.csv")),
                           details = FALSE)$conflation)
})

test_that("results on the risk scale: an animal with only the disease vs a disease-free animal", {
  a <- cull_five(overall_risk = 0.25)$adjusted
  expect_equal(1 - exp(-cull_five(overall_risk = 0.25)$attributable$baseline_hazard), 0.1958530189, tolerance = 1e-7)
  expect_equal(a$risk_alone[1:5], c(0.2063532168, 0.2358884055, 0.2098197199, 0.2260116771, 0.2803435096),
               tolerance = 1e-7)
  expect_equal(a$rr_alone[1:5], c(1.0536126422, 1.2044154682, 1.0713121555, 1.1539861798, 1.4313974384),
               tolerance = 1e-7)
  expect_equal(a$rd_alone[1:5], c(0.0105001978, 0.0400353866, 0.013966701, 0.0301586582, 0.0844904907),
               tolerance = 1e-7)
  expect_equal(a$or_alone[1:5], c(1.0675522706, 1.2675204376, 1.0902479565, 1.1989515543, 1.599449104),
               tolerance = 1e-7)
  # The change is on the raw estimate's measure: HR, RR (MET), OR (SCK), RD (RP).
  expect_equal(a$change[1:5], c(1.060298968 / 1.74, 1.2342866774 / 1.6, 1.0713121555 / 1.23, 1.1989515543 / 1.45,
                                0.0844904907 / 0.10) - 1, tolerance = 1e-8)
  expect_true(is.na(a$change[6]))
})

test_that("a baseline risk instead of the overall risk", {
  r <- cull_five(baseline_risk = 0.18)
  expect_equal(r$adjusted$adjusted, c(1.0609782116, 1.2342456541, 1.0694231539, 1.1773301457, 1.5418427488,
                                      2.2754858679), tolerance = 1e-8)
  s <- r$attributable$summary
  expect_equal(s$disease_free_risk, 0.18, tolerance = 1e-12)
  expect_equal(s$overall_risk, 0.2312436477, tolerance = 1e-7)
  expect_equal(s$attributable, 0.0512436477, tolerance = 1e-7)
  expect_equal(s$raw_attributable, 0.0839363859, tolerance = 1e-7)
  expect_equal(r$risk_given, "baseline")
  expect_true(any(grepl("(given: disease-free risk)", capture.output(print(r)), fixed = TRUE)))
  m <- read_five("culling.csv")
  expect_error(deconflate(m, event_model = TRUE, overall_risk = 0.25, baseline_risk = 0.18),
               class = "deconflate_invalid_input")
  expect_error(deconflate(m, event_model = TRUE), "baseline_risk", class = "deconflate_invalid_input")
  expect_error(deconflate(read_five(), baseline_risk = 0.1), class = "deconflate_unsupported")
  expect_error(deconflate(m, event_model = TRUE, baseline_risk = 1.2), "`baseline_risk`")
  # A distribution, in the draws (key "risk").
  rd <- deconflate(supp_hr_model(), event_model = TRUE, baseline_risk = dist_beta(150, 850), n_draws = 20, seed = 1)
  expect_true("risk" %in% colnames(rd$draws$params))
  expect_equal(rd$risk_given, "baseline")
  # The other tools.
  o <- cm_sensitivity(deconflate(m, event_model = TRUE, baseline_risk = 0.18, n_draws = 0), keys = "risk")
  expect_equal(o$inputs$value, 0.18)
  expect_equal(o$inputs$label, "baseline risk")
  th <- cm_threshold(m, "risk", c(0.05, 0.6), conclusion = "total", target = 0.05, event_model = TRUE,
                     baseline_risk = 0.18, n_grid = 11)
  expect_equal(th$thresholds$status, "threshold")
  h <- supp_hr_model(c(1.05, 2.0, 1.3))
  p <- suppressWarnings(deconflate(h, event_model = TRUE, baseline_risk = 0.15, n_draws = 0, diagnose = TRUE))$problems
  expect_equal(p$flagged$disease[p$flagged$flagged], "d1")
})

test_that("simulated raw estimates of event impacts", {
  h <- supp_hr_model()
  truth <- c(d1 = 1.2, d2 = 1.5, d3 = 1.1)
  sim <- simulate_raw_impacts(h, truth, interactions = c("d1;d2" = 1.4), event_model = TRUE,
                              measure = c("HR", "RR", "OR", "RD"), overall_risk = 0.2)
  expect_equal(sim$measure, c("HR", "RR", "OR", "RD"))
  expect_equal(sim$value, c(1.4266265256, 1.5543228912, 1.2233961428, 0.2203533843), tolerance = 1e-7)
  simb <- simulate_raw_impacts(h, truth, interactions = c("d1;d2" = 1.4), event_model = TRUE,
                               measure = c("HR", "RR", "OR", "RD"), baseline_risk = 0.15)
  expect_equal(simb$value, c(1.4266265256, 1.573184622, 1.2197084762, 0.1930590324), tolerance = 1e-7)
  # Adjusting them recovers the truth.
  back <- deconflate(cm_model(h, sim), event_model = TRUE, overall_risk = 0.2, n_draws = 0)
  expect_equal(back$adjusted$adjusted, c(1.2, 1.5, 1.1, 1.4), tolerance = 1e-8)
  # Sampled animals: rate ratios close to the truth (adjusted for all).
  set.seed(1)
  sn <- simulate_raw_impacts(h, truth, event_model = TRUE, measure = "HR", adjusted_for = "all",
                             overall_risk = 0.05, n = 200000)
  expect_equal(sn$value, unname(truth), tolerance = 0.1)
  expect_equal(sn$source[1], "simulated (n = 200000)")
  expect_error(simulate_raw_impacts(h, c(d1 = -1, d2 = 1, d3 = 1), event_model = TRUE, overall_risk = 0.2),
               class = "deconflate_invalid_input")
  expect_error(simulate_raw_impacts(h, truth, event_model = TRUE, measure = "HZ", overall_risk = 0.2),
               "Unknown measure")
  expect_error(simulate_raw_impacts(h, truth, event_model = TRUE, overall_risk = dist_beta(2, 8)),
               class = "deconflate_invalid_input")
})

test_that("interactions of event impacts in the sensitivity tools", {
  m <- read_five("culling.csv")
  m2 <- set_interaction(m, "MET", "SCK", 1.5)
  im <- m2$impacts
  expect_equal(im$disease[nrow(im)], "MET;SCK")
  expect_equal(im$measure[nrow(im)], "HR")
  expect_equal(im$adjusted_for[nrow(im)], "all")
  # (MET then falls below 1 on its own: the interaction takes its effect)
  r2 <- suppressWarnings(deconflate(m2, event_model = TRUE, overall_risk = 0.25, n_draws = 0))
  expect_equal(r2$adjusted$adjusted[r2$adjusted$disease == "MET;SCK"], 1.5, tolerance = 1e-10)
  expect_error(set_interaction(m, "MET", "SCK", 0), class = "deconflate_invalid_input")
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  si <- cm_sensitivity(res, keys = "interaction:MET;SCK")
  expect_equal(c(si$inputs$low, si$inputs$high), c(0.67, 1.5))
  expect_true(all(is.finite(c(si$inputs$total_low, si$inputs$total_high))))
  expect_equal(si$metric, "attributable risk")
  # A hazard ratio below zero is not a valid interaction: the run is reported as failed.
  bad <- cm_sensitivity(res, keys = "interaction:MET;SCK",
                        ranges = data.frame(input = "interaction:MET;SCK", low = -1, high = 2))
  expect_match(bad$inputs$failed, "low: invalid_input")
  expect_true(is.finite(bad$inputs$total_high))
  th <- cm_threshold(m, "interaction:MET;SCK", c(0.5, 3), conclusion = "total", target = 0.054,
                     event_model = TRUE, overall_risk = 0.25, n_grid = 11)
  expect_true(th$thresholds$status %in% c("threshold", "none"))
})
