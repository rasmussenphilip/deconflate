# deconflate(diagnose = TRUE): which inputs drive flagged adjusted
# impacts, and rows that disagree. Reference values:
# inst/validation/reference_v05_input_check.py (cases A, B, C and E; an
# independent implementation on reference_v041.py). The per-SD effects use a
# central difference (h = 0.01 here, 1e-4 in the reference) and event
# impacts are solved iteratively, hence tolerances of 1e-5 to 1e-4.

effect_of <- function(p, key, disease) {
  e <- p$inputs
  e[e$key == key & e$disease == disease, , drop = FALSE]
}
threshold_of <- function(p, key, disease) {
  th <- p$thresholds
  th$threshold[th$key == key & th$disease == disease]
}

test_that("A: a sign change in additive impacts (five diseases, yield, SCK raw 1.0)", {
  m <- read_five("impacts_yield.csv")
  m$impacts$value[4] <- 1.0
  plain <- suppressWarnings(deconflate(m, n_draws = 0))
  r <- suppressWarnings(deconflate(m, n_draws = 0, diagnose = TRUE))
  p <- r$problems
  expect_s3_class(p, "cm_diagnosis")
  # The check does not change the results.
  expect_equal(r$adjusted, plain$adjusted)
  f <- p$flagged
  expect_equal(f$disease[f$flagged], "SCK")
  expect_equal(f$adjusted[f$disease == "SCK"], -0.1869779815, tolerance = 1e-8)
  expect_equal(f$criterion[f$disease == "SCK"], -0.2869779815, tolerance = 1e-8)
  expect_equal(f$reason[f$disease == "SCK"], "other side of zero")
  ref <- list(c("associations:6", -0.19485945, -1.472743), c("prob:MAS", 0.15768545, 1.819940),
              c("associations:3", -0.13450148, -2.133642), c("impacts:2", -0.11738233, -2.444814),
              c("impacts:1", -0.10920483, -2.627887), c("associations:5", -0.07874080, -3.644591))
  for (x in ref) {
    e <- effect_of(p, x[1], "SCK")
    expect_equal(e$effect_per_sd, as.numeric(x[2]), tolerance = 1e-5)
    expect_equal(e$sd_needed, as.numeric(x[3]), tolerance = 1e-5)
  }
  # Listed: the four largest effects, and every input that alone removes the
  # flag within 3 SD (here the same five inputs).
  expect_setequal(unique(p$inputs$key[p$inputs$listed]),
                  c("associations:6", "prob:MAS", "associations:3", "impacts:2", "impacts:1"))
  # Without a distribution: the SCK row itself (exact; linear; the adjusted
  # impact reaches 10% of the raw impact as entered). The RP row would need
  # a value far outside a quarter to three times its own.
  expect_equal(threshold_of(p, "impacts:4", "SCK"), 1.27457697, tolerance = 1e-7)
  expect_equal(nrow(p$thresholds), 1L)
  expect_true(p$thresholds$in_range)
  expect_length(p$no_single_fix, 0L)
  expect_null(p$disagreeing)
  out <- capture.output(print(r))
  expect_true(any(grepl("SCK: 1.00 -> -0.187 (other side of zero)", out, fixed = TRUE)))
  expect_true(any(grepl("MAS-SCK prevalence", out, fixed = TRUE)))
  expect_true(any(grepl("(1.5 SD down)", out, fixed = TRUE)))
})

test_that("B: event impacts, a sign change and a shrink (five diseases, culling, SCK OR 1.2)", {
  m <- read_five("culling.csv")
  m$impacts$value[4] <- 1.2
  r <- suppressWarnings(deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0, diagnose = TRUE))
  p <- r$problems
  f <- p$flagged
  expect_equal(f$disease[f$flagged], c("LAM", "SCK"))
  expect_equal(f$reason[f$flagged], c("shrunk by more than 90%", "other side of 1"))
  expect_equal(f$adjusted[f$disease == "LAM"], 1.039582872, tolerance = 1e-8)
  expect_equal(f$criterion[f$flagged], c(-0.0165689632, -0.034317179), tolerance = 1e-7)
  ref <- list(c("impacts:6", "LAM", -0.58011059, -0.029692), c("impacts:1", "LAM", 0.45401032, 0.037939),
              c("prob:MAS", "LAM", -0.18553125, -0.092840), c("associations:1", "LAM", -0.08894553, -0.193656),
              c("associations:2", "LAM", -0.07676512, -0.224383),
              c("associations:6", "SCK", -0.04384415, -0.770219), c("prob:MAS", "SCK", 0.03291829, 1.025861),
              c("associations:12", "SCK", -0.01790559, -1.885981), c("impacts:6", "SCK", -0.01746960, -1.933050))
  for (x in ref) {
    e <- effect_of(p, x[1], x[2])
    expect_equal(e$effect_per_sd, as.numeric(x[3]), tolerance = 1e-4)
    expect_equal(e$sd_needed, as.numeric(x[4]), tolerance = 1e-4)
  }
  # Pooled LAM-MAS rows (associations 1 and 2) are moved with the pooling
  # weights of the central estimate.
  th <- list(c("associations:9", "SCK", 0.37067257), c("associations:10", "LAM", 0.44689494),
             c("impacts:2", "LAM", 1.61819061), c("impacts:2", "SCK", 0.79389011),
             c("impacts:4", "LAM", 1.40706670), c("impacts:4", "SCK", 1.24438782))
  for (x in th) expect_equal(threshold_of(p, x[1], x[2]), as.numeric(x[3]), tolerance = 1e-7)
  expect_equal(nrow(p$thresholds), length(th))
  out <- capture.output(print(p))
  expect_true(any(grepl("Flagged adjusted hazard ratios", out, fixed = TRUE)))
})

test_that("E: rows that disagree with other rows (impacts and pooled associations)", {
  m <- read_five("impacts_calving_interval.csv")
  m$impacts$value[3] <- 20
  m$associations$value[2] <- 6
  r <- suppressWarnings(deconflate(m, n_draws = 0, diagnose = TRUE))
  p <- r$problems
  expect_false(any(p$flagged$flagged))
  expect_null(p$inputs)
  dg <- p$disagreeing
  expect_equal(paste0(dg$table, ":", dg$row), c("impacts:2", "impacts:3", "associations:1", "associations:2"))
  expect_equal(dg$std_residual, c(-2.0021579842, 3.0179367530, -3.5848464736, 3.1348891416), tolerance = 1e-7)
  expect_equal(dg$changes_most, rep("MAS", 4))
  expect_equal(dg$before, rep(7.2529895371, 4), tolerance = 1e-8)
  expect_equal(dg$after, c(16.9511095432, 3.0457115821, 6.2784837983, 8.3555657975), tolerance = 1e-7)
  expect_equal(dg$aggregate_after, c(9.4597846429, 6.7984125964, 7.2669557850, 8.1316326692), tolerance = 1e-7)
  out <- capture.output(print(p))
  expect_true(any(grepl("impacts row 3 (MAS, Study 3): 20.0, fitted 10.9, residual +3.0", out, fixed = TRUE)))
})

test_that("nothing to report, and the argument is checked", {
  m <- read_five("impacts_yield.csv")
  r <- deconflate(m, n_draws = 0, diagnose = TRUE)
  expect_false(any(r$problems$flagged$flagged))
  expect_null(r$problems$disagreeing)
  expect_equal(capture.output(print(r$problems)),
               "Diagnosis: no adjusted impact is flagged and no row disagrees with other rows.")
  expect_null(deconflate(m, n_draws = 0)$problems)
  expect_error(deconflate(m, n_draws = 0, diagnose = NA), class = "deconflate_invalid_input")
  expect_error(deconflate(m, n_draws = 0, diagnose = "yes"), class = "deconflate_invalid_input")
})

test_that("the sign-change warning points to the check", {
  m <- read_five("impacts_yield.csv")
  m$impacts$value[4] <- 1.0
  w <- tryCatch(deconflate(m, n_draws = 0), deconflate_sign_change = function(w) w)
  expect_match(conditionMessage(w), "diagnose = TRUE", fixed = TRUE)
})

test_that("labels of inputs", {
  m <- read_five("impacts_calving_interval.csv")
  lab <- input_labels(m, c("prob:MAS", "associations:1", "associations:2", "associations:12", "impacts:2",
                           "impacts:6", "prob:MAS"))
  expect_equal(lab, c("MAS incidence rate", "LAM-MAS OR [row 1]", "LAM-MAS OR [row 2]",
                      "LAM-MAS OR_ratio (among SCK)", "MAS impact (adjusted for SCK)",
                      "RP impact (among no others)", "MAS incidence rate"))
})

test_that("old distribution keys, the pairwise shortcut, and a fixed overall risk", {
  m <- supp_model(c(0.2, 5, 7.5))
  m$distributions <- list("impact:d1" = dist_normal(0.2, 0.1), "assoc:d1:d2" = dist_lognormal_ci(2, 1.4, 2.9))
  r <- suppressWarnings(deconflate(m, n_draws = 0, diagnose = TRUE))
  expect_equal(r$diagnostics$method, "pairwise shortcut")
  p <- r$problems
  expect_equal(p$flagged$disease[p$flagged$flagged], "d1")
  expect_setequal(p$inputs$key, c("impacts:1", "associations:1"))
  expect_true(all(is.finite(p$inputs$effect_per_sd)))
  expect_true(all(c("prob:d2", "associations:2") %in% p$thresholds$key))
  # Event impacts with an overall risk given as a fixed distribution (an
  # input without a distribution) and no other distributions.
  h <- supp_hr_model(c(1.05, 2.0, 1.3))
  r2 <- suppressWarnings(deconflate(h, event_model = TRUE, overall_risk = dist_fixed(0.2), n_draws = 0,
                                    diagnose = TRUE))
  p2 <- r2$problems
  expect_equal(p2$flagged$reason[p2$flagged$flagged], "other side of 1")
  expect_null(p2$inputs)
  expect_true(all(c("associations:1", "impacts:1") %in% p2$thresholds$key))
  expect_true(any(grepl("d1: HR 1.05 -> 0.964 (other side of 1)", capture.output(print(p2)), fixed = TRUE)))
})

test_that("C: the 2024 global dairy yield inputs (reference case C)", {
  skip_on_cran()  # about 230 reruns of a twelve-disease model
  f <- function(x) extdata("global_dairy_2024", x)
  m <- suppressMessages(cm_read_inputs(diseases = f("diseases.csv"), associations = f("associations.csv"),
                                       impacts = f("impacts_yield.csv")))
  p <- suppressWarnings(deconflate(m, n_draws = 0, diagnose = TRUE))$problems
  fl <- p$flagged
  expect_equal(fl$disease[fl$flagged], c("CK", "CM", "DA", "MF"))
  expect_equal(fl$criterion[fl$flagged], c(-5.0739629504, -0.4645487011, -2.8989530754, -2.2638124942),
               tolerance = 1e-8)
  ref <- list(c("impacts:11", "CK", -0.49463560, -10.257982), c("impacts:2", "CM", 0.85490153, 0.543394),
              c("associations:12", "CM", -0.47630539, -0.975317), c("impacts:3", "DA", 2.27913310, 1.271954),
              c("associations:17", "DA", -0.45374369, -6.388966), c("impacts:4", "MF", -0.28673283, -7.895198),
              c("associations:34", "MF", -0.11628324, -19.468089))
  for (x in ref) {
    e <- effect_of(p, x[1], x[2])
    expect_equal(e$effect_per_sd, as.numeric(x[3]), tolerance = 1e-5)
    expect_equal(e$sd_needed, as.numeric(x[4]), tolerance = 1e-5)
  }
  expect_equal(threshold_of(p, "prob:SCM", "CM"), 0.19784578, tolerance = 1e-7)
  expect_equal(threshold_of(p, "associations:8", "CM"), 0.83945405, tolerance = 1e-7)
  expect_equal(threshold_of(p, "impacts:7", "MF"), 2.69413822, tolerance = 1e-7)
  expect_false(p$thresholds$in_range[p$thresholds$key == "impacts:7"])
  expect_equal(p$no_single_fix, c("CK", "MF"))
})
