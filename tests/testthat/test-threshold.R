# Threshold searches (cm_threshold). Reference values:
# inst/validation/reference_v030.py.

test_that("a relative change in the aggregate is found for an association", {
  m <- example_supplement()
  # d1:d3 is independent in the example (baseline OR 1).
  th <- cm_threshold(m, "pair:d1:d3", c(1, 50), conclusion = "total", target = 0.9, relative_to = "adjusted")
  expect_s3_class(th, "cm_threshold")
  expect_equal(th$baseline, 1)
  expect_true(th$log_scale)
  expect_equal(th$method, "pairwise shortcut")
  t1 <- th$thresholds[th$thresholds$status == "threshold", ]
  expect_equal(nrow(t1), 1)
  expect_equal(t1$threshold, 5.263907242734966, tolerance = 1e-6)
  expect_lt(t1$lower, t1$threshold)
  expect_gt(t1$upper, t1$threshold)
  expect_equal(th$summary$result, "threshold(s) found")
  expect_equal(nrow(th$regions), 0)
  # Results just below and above the threshold are returned.
  expect_equal(length(th$details), 1)
  expect_true(all(c("below", "above") %in% names(th$details[[1]])))
  expect_output(print(th), "Crossings")
  th2 <- cm_threshold(m, "pair:d2:d3", c(0.2, 20), conclusion = "total", target = 1.1, relative_to = "adjusted")
  expect_equal(th2$thresholds$threshold[th2$thresholds$status == "threshold"],
               1.4333305414665463, tolerance = 1e-6)
})

test_that("a model-implied sign change of an adjusted impact is found", {
  th <- cm_threshold(example_supplement(), "impacts:1", c(-5, 5), conclusion = "sign",
                     diseases = "d1")
  expect_false(th$log_scale)
  expect_equal(th$thresholds$item, "d1")
  expect_equal(th$thresholds$status, "threshold")
  expect_equal(th$thresholds$threshold, 0.37414199213982113, tolerance = 1e-6)
  expect_match(th$thresholds$description, "model-implied")
  expect_match(th$thresholds$description, "moves above no effect")
})

test_that("varying a probability or an impact draws no random numbers", {
  set.seed(1)
  r1 <- stats::runif(1)
  set.seed(1)
  cm_threshold(example_supplement(), "impacts:1", c(-5, 5), conclusion = "sign", diseases = "d1",
               n_grid = 5)
  cm_threshold(example_supplement(), "prob:d1", c(0.05, 0.5), conclusion = "rank", n_grid = 5)
  expect_equal(stats::runif(1), r1)
})

# p = 0.5 for a, b and c with phi(a, b) = phi(a, c) = 0.8: the conflation
# matrix is the correlation matrix, singular where phi(b, c) = 2 * 0.8^2 - 1
# = 0.28, i.e. at OR(b, c) = (1.28 / 0.72)^2 = 256 / 81. With raw impacts
# 1, 2, 3 the adjusted impact of a is (s - 3) / (s - 0.28) (s = phi(b, c)),
# which changes sign only at the pole. Below OR(b, c) = 16 (phi 0.6) the
# pairs are jointly infeasible, so the pole lies in the infeasible region.
# Reference values: python emulation of cm_threshold() with R's LINPACK
# rank check (qr tolerance 1e-7).
singular_pole_model <- function() {
  cm_model(cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                         cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(0.8, 0.8, 20),
                                         measure = c("phi", "phi", "OR"))),
           cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
}

test_that("an infeasible part of the range is a region, not a crossing", {
  m <- singular_pole_model()
  # The joint fit reports the region below the pole as infeasible.
  th0 <- cm_threshold(m, "pair:b:c", c(0.05, 100), conclusion = "sign", diseases = "a",
                      n_grid = 3)
  expect_equal(th0$scan$status, c("infeasible", "infeasible", "ok"))
  expect_equal(nrow(th0$thresholds), 0)
  expect_equal(th0$regions$n_points, 2L)
  expect_true(th0$log_scale)
  expect_equal(th0$baseline, 20)
  # The feasibility argument is gone.
  expect_error(cm_threshold(m, "pair:b:c", c(0.05, 100), feasibility = "none"),
               class = "deconflate_unsupported")
})

# Crude rows for a, for b and for the combination a;b (all vs none) do not
# tell b_a, b_b and the interaction apart when P(a) + P(b) = 1, whatever the
# association: with P(b) = 0.5, the adjusted impacts have a pole at
# P(a) = 0.5 and change sign across it.
pole_rows_model <- function(pa = 0.3) {
  cm_model(cm_population(cm_diseases(c("a", "b"), c(pa, 0.5)), cm_associations("a", "b", 2)),
           cm_impacts(c("a", "b", "a;b"), c(1, 2, 4)))
}

test_that("a pole of a singular system is a discontinuity, not a threshold", {
  expect_error(deconflate(pole_rows_model(0.5), n_draws = 0), class = "deconflate_singular")
  near <- suppressWarnings(deconflate(pole_rows_model(0.499), n_draws = 0))$adjusted$adjusted
  expect_lt(near[1], -1000)
  # The grid misses the pole; the bisection lands on it.
  th <- cm_threshold(pole_rows_model(), "prob:a", c(0.3, 0.7), conclusion = "sign", diseases = "a",
                     n_grid = 10)
  expect_equal(th$method, "full model")
  expect_false(th$log_scale)
  expect_equal(th$scan$status, rep("ok", 10))
  st <- th$thresholds
  expect_equal(nrow(st), 1)
  expect_equal(st$status, "discontinuity")
  expect_match(st$description, "nearly singular")
  expect_true(is.na(st$threshold))
  expect_lt(st$lower, 0.5)
  expect_gt(st$upper, 0.5)
  expect_lt(st$below, -10)
  expect_gt(st$above, 10)
  expect_equal(th$summary$n_thresholds, 0L)
  expect_equal(th$summary$result, "no threshold; see discontinuities or unresolved crossings")
  # A grid point on the pole is unusable: a region, and no crossing.
  th2 <- cm_threshold(pole_rows_model(), "prob:a", c(0.3, 0.7), conclusion = "sign", diseases = "a",
                      n_grid = 11)
  expect_equal(th2$scan$status[6], "singular")
  expect_equal(th2$regions$status, "singular")
  expect_equal(nrow(th2$thresholds), 0)
  expect_equal(th2$summary$result, "none found where the range could be evaluated")
})

test_that("cm_threshold rejects the retired `method` argument", {
  m <- example_supplement()
  expect_error(cm_threshold(m, "impacts:1", c(-3, 3), conclusion = "sign", diseases = "d1",
                            method = "global"),
               "`method` argument has been removed", class = "deconflate_unsupported")
  # The model stays finite at impact:d1 = -0.5266329181546854, where eq. 16
  # of the 2022 paper divides by zero.
  pole <- -0.5266329181546854
  expect_true(all(is.finite(deconflate(supp_model(c(pole, 5, 7.5)), warn = FALSE)$adjusted$adjusted)))
})

test_that("exact zeros on the grid count only when the sign differs on both sides", {
  f <- c(1, 0, 1, 0, -1, 0, 0, 2, NA, -1, 0)
  ok <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  # 1 -> 0 -> 1 touches zero; 1 -> 0 -> -1 and -1 -> 0 -> 0 -> 2 cross; the
  # sign change across the unusable point 9 and the final zero are not
  # crossings.
  expect_equal(threshold_brackets(f, ok), rbind(c(3L, 5L), c(5L, 8L)))
  expect_equal(threshold_brackets(c(2, -1), c(TRUE, TRUE)), rbind(c(1L, 2L)))
  expect_equal(nrow(threshold_brackets(c(1, 0, NA, -1), c(TRUE, TRUE, FALSE, TRUE))), 0L)
  expect_equal(nrow(threshold_brackets(c(0, 0, 0), rep(TRUE, 3))), 0L)
})

test_that("the aggregate crossing a target is found", {
  th <- cm_threshold(example_supplement(), "impacts:2", c(0, 20), conclusion = "total", target = 3)
  expect_equal(th$thresholds$threshold, 13.03899186228618, tolerance = 1e-6)
  expect_match(th$thresholds$description, "rises above")
  expect_error(cm_threshold(example_supplement(), "impacts:2", c(0, 20), conclusion = "total"),
               "target")
})

test_that("a change in ranking is found, by contribution", {
  th <- cm_threshold(example_supplement(), "pair:d2:d3", c(0.1, 100), conclusion = "rank",
                     n_grid = 201)
  expect_equal(th$summary$item, c("d1 vs d2", "d1 vs d3", "d2 vs d3"))
  expect_equal(th$summary$n_thresholds, c(1L, 0L, 0L))
  expect_equal(nrow(th$thresholds), 1)
  t1 <- th$thresholds[th$thresholds$item == "d1 vs d2", ]
  expect_equal(t1$threshold, 10.20365117522034, tolerance = 1e-6)
  expect_match(t1$description, "d1 moves above d2")
  expect_equal(th$summary$result, c("threshold(s) found", "none found in the range",
                                    "none found in the range"))
})

test_that("unusable parts of the range are reported, not mistaken for thresholds", {
  # A probability of 1 or more is invalid. Below 1 every grid point passes
  # the triple screen (reference_v030.py port), so the only unusable region
  # is P(d1) >= 1.
  th <- cm_threshold(example_supplement(), "prob:d1", c(0.05, 1.5), conclusion = "rank")
  expect_equal(nrow(th$regions), 1)
  expect_equal(th$regions$status, "infeasible")
  expect_equal(th$regions$from, 1.007, tolerance = 1e-12)
  expect_equal(th$regions$to, 1.5)
  expect_equal(th$regions$n_points, 35L)
  expect_equal(sort(unique(th$scan$status)), c("infeasible", "ok"))
  expect_true(all(th$scan$status[th$scan$value >= 1] == "infeasible"))
  expect_true(all(th$scan$status[th$scan$value < 1] == "ok"))
  expect_equal(th$summary$n_thresholds, c(1L, 1L, 0L))
  expect_equal(th$thresholds$item, c("d1 vs d2", "d1 vs d3"))
  expect_equal(th$thresholds$status, c("threshold", "threshold"))
  expect_equal(th$thresholds$threshold, c(0.22362189686298367, 0.6220804451704026), tolerance = 1e-6)
  expect_equal(th$summary$result, c("threshold(s) found", "threshold(s) found",
                                    "none found where the range could be evaluated"))
  expect_output(print(th), "could not be used")
})

test_that("an invalid value of the input itself is an infeasible point", {
  # Odds ratios <= 0 on a linear grid.
  th <- cm_threshold(example_supplement(), "pair:d1:d2", c(-1, 5), conclusion = "rank",
                     log_scale = FALSE)
  expect_equal(th$regions$status, "infeasible")
  expect_equal(th$regions$n_points, 17L)
  expect_equal(nrow(th$thresholds), 0)
  expect_true(all(th$summary$result == "none found where the range could be evaluated"))
})

test_that("interaction inputs fit the full model", {
  skip_on_cran()
  th <- cm_threshold(example_supplement(), "interaction:d1;d2", c(-1, 1), conclusion = "total",
                     target = 2.12, n_grid = 11)
  # Every point (with or without the interaction row) uses one joint fit.
  expect_equal(th$method, "full model")
  expect_equal(th$baseline, 0)
  # The aggregate falls as the interaction grows (crude impacts absorb it).
  expect_equal(th$thresholds$status, "threshold")
  expect_equal(th$thresholds$threshold, -0.5347852259874342, tolerance = 1e-5)
  expect_match(th$thresholds$description, "falls below")
  # Ranked by contribution (interaction shares included).
  thr <- cm_threshold(example_supplement(), "interaction:d1;d2", c(-1, 1), conclusion = "rank", n_grid = 5)
  expect_equal(thr$summary$item, c("d1 vs d2", "d1 vs d3", "d2 vs d3"))
  # An interaction of event impacts is a hazard ratio (an HR row adjusted
  # for all), so the range must be positive.
  expect_error(cm_threshold(supp_hr_model(), "interaction:d1;d2", c(-1, 1), conclusion = "total",
                            target = 0.1, event_model = TRUE, overall_risk = 0.25),
               class = "deconflate_invalid_input")
})

test_that("rows are varied by number and by alias", {
  m <- example_supplement()
  a <- cm_threshold(m, "pair:d2:d3", c(0.2, 20), conclusion = "total", target = 1.1, relative_to = "adjusted", n_grid = 21)
  b <- cm_threshold(m, "associations:3", c(0.2, 20), conclusion = "total", target = 1.1, relative_to = "adjusted", n_grid = 21)
  expect_equal(a$thresholds$threshold, b$thresholds$threshold)
  expect_equal(a$baseline, b$baseline)
  i1 <- cm_threshold(m, "impacts:2", c(0, 20), conclusion = "total", target = 3, n_grid = 11)
  i2 <- cm_threshold(m, "impacts:2", c(0, 20), conclusion = "total", target = 3, n_grid = 11)
  expect_equal(i1$thresholds$threshold, 13.03899186228618, tolerance = 1e-6)
  expect_equal(i2$thresholds$threshold, i1$thresholds$threshold)
  # A three-way ratio changes additive results only through combinations:
  # with an interaction it moves the total.
  mi <- set_interaction(m, "d1", "d2", 1)
  t3 <- cm_threshold(mi, "or_ratio:d1:d2|d3", c(0.5, 2), conclusion = "total", target = 2.0890878497857703,
                     n_grid = 5)
  expect_true(is.na(t3$baseline))   # no OR_ratio row in the model
  expect_equal(t3$scan$total[c(1, 5)], c(2.0913976662985534, 2.086745618267496), tolerance = 1e-7)
  # Two crude rows of a pair: assoc: is ambiguous; vary a row by number.
  m5 <- read_five()
  expect_error(cm_threshold(m5, "pair:LAM:MAS", c(1, 3)), "associations:")
  expect_error(cm_threshold(m, "impacts:9", c(0, 1)), "row")
})

test_that("descriptions of relative changes give the level", {
  th <- cm_threshold(example_supplement(), "pair:d2:d3", c(0.2, 20), conclusion = "total", target = 1.1, relative_to = "adjusted")
  expect_match(th$thresholds$description, "the aggregate (rises above|falls below) 1.1 times its adjusted value \\(")
})

test_that("inputs are validated", {
  m <- example_supplement()
  expect_error(cm_threshold(m, "pair:d1", c(1, 2)), "wrong number of parts")
  expect_error(cm_threshold(m, "pair:d1:zz", c(1, 2)), "Unknown diseases")
  expect_error(cm_threshold(m, "foo:d1", c(1, 2)), "Unknown input")
  expect_error(cm_threshold(m, "impacts:1", c(2, 1)), "increasing")
  expect_error(cm_threshold(m, "impacts:1", c(-1, 1), log_scale = TRUE), "positive range")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), conclusion = "rank", diseases = "d1"),
               "at least two")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), conclusion = "rank", diseases = c("d1", "d1")),
               "at least two")
  expect_error(cm_threshold(m, "pair:d1:d1", c(1, 2)), "same disease twice")
  expect_error(cm_threshold(m, c("impacts:1", "impacts:2"), c(1, 2)), "`range` is for one input")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), conclusion = "change", target = 0.1),
               "relative_to = 'adjusted'", class = "deconflate_unsupported")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), conclusion = "rank", target = 1), "only with conclusion")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), conclusion = "sign", relative_to = "raw"), "only with conclusion")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), conclusion = "total", target = 1, relative_to = "x"),
               class = "deconflate_invalid_input")
  expect_error(cm_threshold(m, conclusion = "rank", top = 0), "top")
  expect_error(cm_threshold(m, conclusion = "rank", sensitivity = list()), "cm_sensitivity")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), n_grid = NA), "n_grid")
  expect_error(cm_threshold(m, "impacts:1", c(0, 1), log_scale = NA), "log_scale")
  expect_error(cm_threshold(supp_model(c(0, 0, 0)), "impacts:1", c(0, 1), conclusion = "total",
                            target = 1.1, relative_to = "adjusted"), "zero")
  expect_error(cm_threshold(m, "prob:d1", c(0.05, 0.5), joint = fit_joint(m)), "fixed `joint`")
})

test_that("threshold plots draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  th <- cm_threshold(example_supplement(), "pair:d2:d3", c(0.1, 100), conclusion = "rank")
  expect_invisible(plot(th))
})

test_that("the printout fits 80 columns and says when the input was not entered", {
  m <- cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5)))
  th <- cm_threshold(m, "pair:d1:d3", c(0.2, 5), conclusion = "rank")
  out <- capture.output(print(th))
  expect_true(all(nchar(out) <= 80))
  expect_match(printed(th), "(now not entered;", fixed = TRUE)
})

# ---- v0.6 item 12: targets for sets of diseases, relative targets ---------
# Reference values: inst/validation/reference_v06_thresholds.py.

test_that("the loss of a set of diseases crosses a target", {
  m <- read_five("impacts_yield.csv")
  set <- c("LAM", "MAS", "SCK")
  a1 <- cm_threshold(m, "associations:3", c(1, 4), conclusion = "total", diseases = set, target = 1.95, n_grid = 31)
  expect_equal(a1$thresholds$threshold, 3.5001209402, tolerance = 1e-6)
  expect_equal(a1$thresholds$item, "total of LAM, MAS, SCK")
  expect_match(a1$thresholds$description, "the loss of LAM, MAS and SCK falls below 1.95", fixed = TRUE)
  a2 <- cm_threshold(m, "associations:3", c(1, 4), conclusion = "total", diseases = set, target = 0.97,
                     relative_to = "adjusted", n_grid = 31)
  expect_equal(a2$reference, 2.1011883335, tolerance = 1e-7)
  expect_equal(a2$level, 0.97 * 2.1011883335, tolerance = 1e-7)
  expect_equal(a2$thresholds$threshold, 2.4933147198, tolerance = 1e-6)
  # Relative to the raw value (the diseases unrelated: p * raw here).
  a3 <- cm_threshold(m, "prob:SCK", c(0.2, 0.6), conclusion = "total", diseases = set, target = 0.75,
                     relative_to = "raw", n_grid = 31)
  expect_equal(a3$reference, 2.9302998718, tolerance = 1e-8)
  expect_equal(a3$thresholds$threshold, 0.4347663073, tolerance = 1e-6)
  a4 <- cm_threshold(m, "prob:SCK", c(0.2, 0.6), conclusion = "total", target = 0.75, relative_to = "raw",
                     n_grid = 31)
  expect_equal(a4$reference, 3.7422998718, tolerance = 1e-8)
  expect_equal(a4$thresholds$threshold, 0.4319971303, tolerance = 1e-6)
  # All diseases named is the aggregate itself.
  a5 <- cm_threshold(m, "prob:SCK", c(0.2, 0.6), conclusion = "total", target = 0.75, relative_to = "raw",
                     diseases = c("LAM", "MAS", "MET", "SCK", "RP"), n_grid = 31)
  expect_equal(a5$thresholds$threshold, a4$thresholds$threshold)
})

test_that("raw contributions with interactions come from the refit with the diseases unrelated", {
  w <- read_five("impacts_welfare.csv")
  res <- deconflate(w, n_draws = 0)
  rc <- raw_total(w, res$settings$plan, NULL, res$row_se, by_disease = TRUE)
  expect_equal(unname(rc), c(2.5992670552, 1.4902952311, 0.4955, 1.2140625, 0.2355), tolerance = 1e-8)
  expect_equal(sum(rc), res$totals$raw_total, tolerance = 1e-10)
  # The result keeps them (and contribution_table() shows them).
  expect_equal(res$contributions$raw_contribution, unname(rc[res$contributions$disease]), tolerance = 1e-12)
  expect_equal(contribution_table(res)$raw_contribution, unname(rc[contribution_table(res)$disease]),
               tolerance = 1e-12)
  b1 <- cm_threshold(w, "impacts:6", c(0, 6), conclusion = "total", diseases = c("LAM", "MAS"), target = 0.8,
                     relative_to = "raw", n_grid = 31)
  expect_equal(b1$reference, 4.0895622863, tolerance = 1e-8)
  expect_equal(b1$thresholds$threshold, 2.8999084221, tolerance = 1e-6)
  b1b <- cm_threshold(w, "impacts:1", c(5, 15), conclusion = "total", diseases = c("LAM", "MAS"), target = 0.85,
                      relative_to = "raw", n_grid = 31)
  expect_equal(b1b$thresholds$threshold, 10.3887986429, tolerance = 1e-6)
  b2 <- cm_threshold(w, "prob:LAM", c(0.1, 0.4), conclusion = "total", diseases = c("LAM", "MAS"), target = 3,
                     n_grid = 31)
  expect_equal(b2$thresholds$threshold, 0.2082453761, tolerance = 1e-6)
})

test_that("event impacts: the attributable risk of a set of diseases", {
  m <- read_five("culling.csv")
  r <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  rc <- raw_total(m, r$settings$plan, 0.25, r$row_se, by_disease = TRUE)
  expect_equal(unname(rc), c(0.0245138239, 0.0232707068, 0.0049379599, 0.0225485554, 0.0053104297),
               tolerance = 1e-8)
  expect_equal(sum(rc), r$attributable$summary$raw_attributable, tolerance = 1e-8)
  expect_equal(r$contributions$raw_contribution, unname(rc[r$contributions$disease]), tolerance = 1e-12)
  set <- c("LAM", "MAS")
  c1 <- cm_threshold(r, "impacts:2", c(0.8, 2.4), conclusion = "total", diseases = set, target = 0.03, n_grid = 31)
  expect_equal(c1$thresholds$threshold, 1.5163936339, tolerance = 1e-6)
  expect_match(c1$thresholds$description, "the attributable risk of LAM and MAS rises above 0.0300", fixed = TRUE)
  c2 <- cm_threshold(r, "impacts:6", c(2, 5), conclusion = "total", diseases = set, target = 0.9,
                     relative_to = "adjusted", n_grid = 31)
  expect_equal(c2$reference, 0.0355830961, tolerance = 1e-8)
  expect_equal(c2$thresholds$threshold, 3.3069022137, tolerance = 1e-6)
  c3 <- cm_threshold(r, "risk", c(0.15, 0.35), conclusion = "total", diseases = set, target = 0.6,
                     relative_to = "raw", n_grid = 31)
  expect_equal(c3$reference, 0.0477845307, tolerance = 1e-8)
  expect_equal(c3$thresholds$threshold, 0.1880800223, tolerance = 1e-6)
})

test_that("default ranges follow the rules for each kind of input", {
  m <- read_five("impacts_yield.csv")
  res <- deconflate(m, n_draws = 0)
  th <- threshold_setup(res, FALSE, NULL, NULL, list())
  get <- function() res
  # With a distribution: its 0.5% and 99.5% quantiles (LAM-SCK OR: normal(2, 0.3) truncated at 0).
  r3 <- threshold_default_range("associations:3", th, get)
  expect_equal(c(r3$low, r3$high), stats::qnorm(c(0.005, 0.995), 2, 0.3), tolerance = 1e-8)
  # The pair key of an entered crude row is that row.
  expect_equal(threshold_default_range("pair:LAM:SCK", th, get)[c("low", "high")], r3[c("low", "high")])
  # Without one: -50% to +50% (an OR of 1.5 for MAS-RP).
  r12 <- threshold_default_range("associations:11", th, get)
  expect_equal(c(r12$low, r12$high, r12$value), c(0.75, 2.25, 1.5))
  # Probabilities stay inside (0, 1).
  pr <- threshold_default_range("prob:RP", th, get)
  expect_equal(c(pr$low, pr$high), c(0.03, 0.09))
  # Unknown pairs and interactions.
  mu <- m
  mu$associations <- mu$associations[-11, ]
  thu <- threshold_setup(mu, FALSE, NULL, NULL, list())
  expect_equal(unlist(threshold_default_range("pair:MAS:RP", thu, get)[c("low", "high", "how")]),
               c(low = "0.2", high = "5", how = "OR 0.2-5"))
  ri <- threshold_default_range("interaction:LAM;MAS", th, get)
  w <- min(abs(res$adjusted$adjusted[1:2]))
  expect_equal(c(ri$low, ri$high), c(-w, w))
  # A named input without a range uses its default range.
  t1 <- cm_threshold(res, "pair:LAM:SCK", conclusion = "rank", diseases = c("MAS", "SCK"), n_grid = 11)
  expect_equal(t1$range, c(r3$low, r3$high))
})

test_that("inputs are chosen by their effect on the compared quantity", {
  m <- read_five("impacts_yield.csv")
  res <- deconflate(m, n_draws = 0)
  s <- cm_sensitivity(res)
  # total of a set: the swing in the sum of its contributions
  set <- c("LAM", "SCK")
  bd <- s$by_disease
  eff <- vapply(s$inputs$input, function(k) {
    b <- bd[bd$input == k & bd$disease %in% set, ]
    abs(sum(b$contribution_high) - sum(b$contribution_low))
  }, numeric(1))
  eff <- eff[is.finite(eff) & eff > 0]
  top <- names(sort(eff, decreasing = TRUE))[1:2]
  ts <- cm_threshold(res, conclusion = "total", diseases = set, target = 1.3, top = 2, sensitivity = s, n_grid = 11)
  expect_s3_class(ts, "cm_threshold_set")
  expect_equal(ts$inputs$input, top)
  expect_equal(ts$inputs$effect, unname(eff[top]), tolerance = 1e-12)
  expect_equal(names(ts$searches), top)
  expect_true(all(vapply(ts$searches, inherits, logical(1), "cm_threshold")))
  expect_equal(names(ts$thresholds)[1:2], c("input", "label"))
  # rank of two diseases: the swing in the difference of their contributions
  rk <- vapply(s$inputs$input, function(k) {
    b <- bd[bd$input == k, ]
    d <- function(col) b[[col]][b$disease == "MAS"] - b[[col]][b$disease == "SCK"]
    abs(d("contribution_high") - d("contribution_low"))
  }, numeric(1))
  rk <- rk[is.finite(rk) & rk > 0]
  tr <- cm_threshold(res, conclusion = "rank", diseases = c("MAS", "SCK"), top = 3, sensitivity = s, n_grid = 11)
  expect_equal(tr$inputs$input, names(sort(rk, decreasing = TRUE))[1:3])
  out <- capture.output(print(tr))
  expect_true(all(nchar(out) <= 80))
  expect_match(out[1], "where MAS and SCK swap places", fixed = TRUE)
  # Several named inputs.
  tn <- cm_threshold(res, c("prob:SCK", "pair:LAM:SCK"), conclusion = "total", target = 2.8, n_grid = 11)
  expect_equal(tn$inputs$input, c("prob:SCK", "pair:LAM:SCK"))
  expect_true(all(is.na(tn$inputs$effect)))
  expect_false(tn$auto)
  # Only unknown inputs.
  tu <- cm_threshold(res, conclusion = "total", target = 2.8, inputs = "unknown", top = 1, n_grid = 11)
  expect_equal(tu$inputs$kind, "unknown")
})
