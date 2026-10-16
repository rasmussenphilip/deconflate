# Shapley attribution (R/shapley.R) and sensitivity / scenario tools
# (R/sensitivity.R): cm_sensitivity() and the set_*() helpers.
# Supplement impacts are in percent (2.5, 5, 7.5);
# reference values from inst/validation/reference_v020_tests.py,
# reference_2022_shapley.py (multiplicative Shapley values) and
# reference_v040.py (sections 6-8: unknown pairs, models without
# associations, event models).

t2_ids <- c("d1", "d2", "d3")
t2_pop <- function() {
  cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)))
}
t2_model <- function(values) cm_model(t2_pop(), cm_impacts(t2_ids, values))
# Supplement population with culling hazard ratios 1.5, 2.0 and 1.3.
t2_event <- function() {
  cm_model(t2_pop(), cm_impacts(t2_ids, c(1.5, 2.0, 1.3), measure = "HR"))
}
# The supplement without any association estimates.
t2_no_assoc <- function(m = example_supplement()) {
  m$associations <- NULL
  m
}

supp_truth <- c(d1 = 2, d2 = 4, d3 = 6)
supp_base_total <- 2.109228876453694       # adjusted aggregate (pairwise shortcut), percent
supp_int_total <- 2.0890878497857703       # global, interaction d1:d2 = 1
supp_int_tw_total <- c(`0.5` = 2.0913976662985534, `2` = 2.086745618267496)
supp_ar <- 0.03664898509534775             # attributable risk, event model at 0.25

test_that("additive contributions add up to the aggregate", {
  m <- example_supplement()
  res <- deconflate(m)
  expect_equal(res$diagnostics$method, "pairwise shortcut")
  expect_equal(sum(res$contributions$contribution), res$totals$adjusted_total)
  expect_equal(res$totals$adjusted_total, supp_base_total, tolerance = 1e-8)
  res_i <- deconflate(set_interaction(m, "d1", "d2", 1))
  expect_equal(res_i$diagnostics$method, "full model")
  expect_equal(res_i$adjusted$kind, c(rep("disease", 3), "interaction"))
  expect_equal(res_i$adjusted$adjusted[4], 1)
  expect_equal(sum(res_i$contributions$contribution), res_i$totals$adjusted_total)
  expect_equal(res_i$totals$adjusted_total, supp_int_total, tolerance = 1e-8)
  P12 <- res_i$joint_pairs["d1", "d2"]
  expect_equal(res_i$totals$interaction_total, P12)
  expect_equal(res_i$contributions$interaction, c(0.5 * P12, 0.5 * P12, 0))
  expect_equal(res_i$totals$adjusted_total,
               sum(c(0.10, 0.15, 0.20) * res_i$adjusted$adjusted[1:3]) + P12, tolerance = 1e-12)
})

test_that("cell-wise Shapley equals the closed form for pairwise interactions", {
  pop <- t2_pop()
  m2 <- cm_model(pop, simulate_raw_impacts(pop, supp_truth, interactions = c("d1;d2" = 1)))
  res <- deconflate(m2)
  expect_equal(res$adjusted$adjusted, c(unname(supp_truth), 1), tolerance = 1e-8)
  closed <- res$contributions
  adj <- stats::setNames(res$adjusted$adjusted[1:3], res$adjusted$disease[1:3])
  # Additive loss: main effects plus the d1:d2 interaction (1).
  loss <- function(x) sum(adj[names(x)] * x) + 1 * x[["d1"]] * x[["d2"]]
  expect_silent(sh <- shapley_by_cell(res$joint, loss))
  expect_equal(sh$disease, t2_ids)
  expect_equal(sh$shapley, closed$contribution, tolerance = 1e-7)
  expect_equal(sum(sh$shapley), attr(sh, "total"), tolerance = 1e-12)
  expect_equal(attr(sh, "total"), res$totals$adjusted_total, tolerance = 1e-8)
  expect_equal(attr(sh, "skipped_mass"), 0)
  expect_equal(attr(sh, "skipped_loss"), 0)
  expect_equal(sum(sh$share), 1)
})

test_that("higher-order additive terms are split equally", {
  j <- fit_joint(t2_pop())
  main <- c(d1 = 2, d2 = 4, d3 = 6)
  loss <- function(x) sum(main[names(x)] * x) + 3 * prod(x[c("d1", "d2", "d3")])
  sh <- shapley_by_cell(j, loss)
  p123 <- sum(j$prob[rowSums(j$cells) == 3])
  expect_equal(sh$shapley, unname(main * c(0.10, 0.15, 0.20) + 3 * p123 / 3), tolerance = 1e-8)
  expect_error(shapley_by_cell(t2_pop(), loss), "fit_joint")
})

# Multiplicative loss: 1 - prod(1 - m_i) over the diseases present.
t2_mult_loss <- function(m = c(d1 = 0.02, d2 = 0.04, d3 = 0.06)) {
  function(x) 1 - prod((1 - m[names(x)])^x)
}

test_that("multiplicative losses match the Python reference", {
  j <- fit_joint(t2_pop())
  sh <- shapley_by_cell(j, t2_mult_loss())
  expect_equal(sh$shapley, c(0.0019783456, 0.0059222732, 0.0119200649), tolerance = 1e-7)
  expect_equal(attr(sh, "total"), 0.01982068370187362, tolerance = 1e-9)
  expect_equal(attr(sh, "allocated"), sum(sh$shapley))
})

test_that("skipped combinations are reported, not silently dropped", {
  j <- fit_joint(t2_pop())
  f <- t2_mult_loss()
  expect_warning(sh <- shapley_by_cell(j, f, max_present = 2), class = "deconflate_incomplete_allocation")
  p123 <- j$prob[rowSums(j$cells) == 3]
  expect_equal(attr(sh, "skipped_mass"), p123)
  expect_equal(attr(sh, "skipped_loss"), p123 * f(c(d1 = 1, d2 = 1, d3 = 1)))
  expect_gt(attr(sh, "skipped_loss"), 0)
  expect_equal(attr(sh, "allocated"), sum(sh$shapley))
  expect_equal(attr(sh, "allocated") + attr(sh, "skipped_loss"), attr(sh, "total"), tolerance = 1e-12)
  # The total still covers every combination.
  expect_equal(attr(sh, "total"), 0.01982068370187362, tolerance = 1e-9)
  expect_equal(sum(sh$share), 1)
})

test_that("set_three_way, set_interaction and set_association edit models", {
  base <- example_supplement()
  tw <- set_three_way(base, "d1", "d2", "d3", 2)
  a <- tw$associations
  expect_equal(nrow(a), 4)
  expect_equal(a$measure[4], "OR_ratio")
  expect_equal(c(a$disease[4], a$associated_with[4], a$among[4], a$note[4]), c("d1", "d2", "d3", "scenario"))
  # The ratio is symmetric in the three diseases: the same row is replaced.
  tw2 <- set_three_way(tw, "d3", "d1", "d2", 0.5)
  expect_equal(nrow(tw2$associations), 4)
  expect_equal(tw2$associations$value[4], 0.5)
  expect_equal(tw2$associations$row, 1:4)
  expect_error(set_three_way(base, "d1", "d2", "zz", 2), "Unknown diseases")
  expect_error(set_three_way(base, "d1", "d2", "d3", -1), "positive")
  expect_error(set_three_way(base, "d1", "d2", "d1", 2), "three different")
  pop_tw <- set_three_way(t2_pop(), "d1", "d2", "d3", 2)
  expect_s3_class(pop_tw, "cm_population")

  wi <- set_interaction(base, "d1", "d2", 1)
  expect_equal(nrow(wi$impacts), 4)
  expect_equal(c(wi$impacts$disease[4], wi$impacts$adjusted_for[4]), c("d1;d2", "all"))
  wi2 <- set_interaction(wi, "d2", "d1", 2)
  expect_equal(nrow(wi2$impacts), 4)
  expect_equal(wi2$impacts$value[4], 2)
  wi3 <- set_interaction(wi2, "d2", "d3", 0.5)
  expect_equal(nrow(wi3$impacts), 5)
  expect_equal(wi3$impacts$disease[5], "d2;d3")
  expect_equal(deconflate(wi3, n_draws = 0)$adjusted$adjusted[4:5], c(2, 0.5), tolerance = 1e-8)
  expect_error(set_interaction(t2_pop(), "d1", "d2", 1), "cm_model")
  expect_error(set_interaction(base, "d1", "zz", 1), "Unknown diseases")
  expect_error(set_interaction(base, "d1", "d1", 1), "two different")
  expect_error(set_association(base, "zz", "d1", 2), "Unknown diseases")

  alt <- set_association(base, "d3", "d1", 2)
  a <- alt$associations
  expect_equal(nrow(a), 3)
  hit <- pair_row_hits(a, "d1", "d3")
  expect_equal(sum(hit), 1)
  expect_equal(a$value[hit], 2)
  expect_equal(a$note[hit], "scenario")
  expect_equal(a$row, 1:3)
  new_pop <- set_association(cm_population(cm_diseases(t2_ids, c(0.1, 0.15, 0.2))), "d1", "d2", 2)
  expect_equal(nrow(new_pop$associations), 1)
  expect_equal(pair_tables(new_pop)$status, c("specified", "unknown", "unknown"))
  # Several crude rows of a pair are replaced by one; rows within subgroups
  # are kept, and the distributions of replaced rows are dropped.
  m <- read_five()
  lm <- m$associations$row[pair_row_hits(m$associations, "LAM", "MAS")]
  expect_length(lm, 2)
  m2 <- set_association(m, "MAS", "LAM", 3)
  hit <- pair_row_hits(m2$associations, "LAM", "MAS")
  expect_equal(sum(hit), 1)
  expect_equal(m2$associations$value[hit], 3)
  expect_equal(nrow(m2$associations), nrow(m$associations) - 1)
  expect_true("OR_ratio" %in% m2$associations$measure)
  expect_false(any(paste0("associations:", lm) %in% names(m2$distributions)))
  # Retired measures cannot be set.
  expect_error(set_association(base, "d1", "d3", 1, measure = "independent"), class = "deconflate_unsupported")
})

# The sensitivity tool (cm_sensitivity(); more in test-sensitivity.R), with
# reference values of the earlier screens: each scenario is one end of an
# input's range.
sens_total <- function(s, key, end = c("low", "high")) {
  end <- match.arg(end)
  s$inputs[[paste0("total_", end)]][s$inputs$input == key]
}

test_that("entered associations: the supplement's scenarios", {
  res <- deconflate(example_supplement(), n_draws = 0)
  s <- cm_sensitivity(res, inputs = "entered", keys = c("associations:2", "associations:3"),
                      ranges = data.frame(input = c("associations:2", "associations:3"), low = c(0.5, 1.5),
                                          high = c(2, 6)))
  expect_s3_class(s, "cm_sensitivity")
  expect_equal(s$central$total, supp_base_total, tolerance = 1e-10)
  expect_equal(s$metric, "aggregate")
  # d1:d3 (OR 1) at 2; d2:d3 (OR 3) at 1.5.
  expect_equal(sens_total(s, "associations:2", "high"), 2.0128810813034628, tolerance = 1e-10)
  expect_equal(sens_total(s, "associations:3", "low"), 2.306649129926649, tolerance = 1e-10)
  expect_equal(s$inputs$range, c("given", "given"))
  expect_true(all(is.na(s$inputs$failed)))
  expect_output(print(s), "Sensitivity of the aggregate")
  # The same scenario adjusted directly.
  d13 <- deconflate(set_association(example_supplement(), "d1", "d3", 2), n_draws = 0)
  expect_equal(d13$totals$adjusted_total, 2.0128810813034628, tolerance = 1e-10)
  # The aggregate is linear in the impacts: impact rows give symmetric totals.
  si <- cm_sensitivity(res, keys = c("impacts:1", "impacts:2", "impacts:3"))
  expect_equal((si$inputs$total_low + si$inputs$total_high) / 2, rep(supp_base_total, 3), tolerance = 1e-10)
  expect_equal(si$inputs$range, rep("20%", 3))
})

test_that("runs that cannot be adjusted are kept with the reason", {
  m <- cm_model(cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                              cm_associations(c("d1", "d2"), c("d2", "d3"), c(0.3, 3),
                                              measure = c("prevalence", "OR"))),
                cm_impacts(t2_ids, c(2.5, 5, 7.5)))
  # P(d1 | d2) = 0.9 is a valid value but impossible with P(d1) = 0.10 and
  # P(d2) = 0.15 (0.9 x 0.15 > 0.10); a negative odds ratio is invalid.
  s <- cm_sensitivity(deconflate(m, n_draws = 0), keys = c("associations:1", "associations:2"),
                      ranges = data.frame(input = c("associations:1", "associations:2"), low = c(0.15, -3),
                                          high = c(0.9, 1.5)))
  f <- stats::setNames(s$inputs$failed, s$inputs$input)
  expect_equal(unname(f["associations:1"]), "high: infeasible")
  expect_equal(unname(f["associations:2"]), "low: infeasible")
  expect_true(is.finite(sens_total(s, "associations:1", "low")))
  expect_true(is.na(s$inputs$swing[s$inputs$input == "associations:1"]))
  # Inputs without a swing come last.
  expect_true(all(is.na(utils::tail(s$inputs$swing, 2))))
})

test_that("the aggregate is reported whatever its size", {
  # Impacts 45 times the supplement's (percent): 103.8% when the d2:d3 odds
  # ratio is halved (reference_v020_tests.py screen values, scaled).
  m <- t2_model(45 * c(2.5, 5, 7.5))
  s <- cm_sensitivity(deconflate(m, n_draws = 0), keys = "associations:3",
                      ranges = data.frame(input = "associations:3", low = 1.5, high = 6))
  expect_equal(s$central$total, 45 * supp_base_total, tolerance = 1e-9)
  expect_equal(sens_total(s, "associations:3", "low"), 45 * 2.306649129926649, tolerance = 1e-9)
})

test_that("interactions: every pair without a combination row, at the default range", {
  res <- deconflate(example_supplement(), n_draws = 0)
  s <- cm_sensitivity(res, inputs = "unknown")
  expect_setequal(s$inputs$input, c("interaction:d1;d2", "interaction:d1;d3", "interaction:d2;d3"))
  half <- 0.5 * min(abs(res$adjusted$adjusted[1:2]))
  expect_equal(s$inputs$high[s$inputs$input == "interaction:d1;d2"], half)
  expect_equal(s$inputs$low[s$inputs$input == "interaction:d1;d2"], -half)
  expect_true(all(is.finite(s$inputs$swing)))
  s1 <- cm_sensitivity(res, keys = "interaction:d1;d2",
                       ranges = data.frame(input = "interaction:d1;d2", low = 0.5, high = 1))
  expect_equal(sens_total(s1, "interaction:d1;d2", "high"), supp_int_total, tolerance = 1e-8)
  expect_error(cm_sensitivity(res, keys = "interaction:d1;zz"), "Unknown input keys")
  expect_error(cm_sensitivity(res, ranges = data.frame(input = "pair:d1:d2", low = 1, high = 2)),
               "No input 'pair:d1:d2'")
})

test_that("three-way terms change global results with interactions only", {
  base <- example_supplement()
  st <- cm_sensitivity(deconflate(base, n_draws = 0), inputs = "unknown", three_way = list(c("d1", "d2", "d3")),
                       keys = "or_ratio:d1:d2|d3")
  expect_equal(st$inputs$range, "ratio 0.5-2")
  expect_lt(abs(st$inputs$swing), 1e-7 * supp_base_total)
  tw <- set_three_way(base, "d1", "d2", "d3", 2)
  res_tw0 <- deconflate(tw)
  expect_equal(res_tw0$diagnostics$method, "full model")
  expect_equal(res_tw0$adjusted$adjusted, deconflate(base)$adjusted$adjusted, tolerance = 1e-7)
  # With an interaction, the probability of the triple matters.
  int <- set_interaction(base, "d1", "d2", 1)
  sti <- cm_sensitivity(deconflate(int, n_draws = 0), inputs = "unknown", three_way = "all")
  expect_equal(sti$central$total, supp_int_total, tolerance = 1e-8)
  expect_equal(c(sens_total(sti, "or_ratio:d1:d2|d3", "low"), sens_total(sti, "or_ratio:d1:d2|d3", "high")),
               unname(supp_int_tw_total), tolerance = 1e-7)
  expect_error(cm_sensitivity(deconflate(base, n_draws = 0), three_way = list(c("d1", "d1", "d2"))),
               "three different diseases")
})

test_that("unknown pairs fit the full model", {
  m <- cm_model(cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25)),
                cm_impacts(c("a", "b", "c"), c(3, 4, 5)),
                associations = cm_associations(c("a", "b"), c("b", "c"), c(3, 2)))
  s <- cm_sensitivity(deconflate(m, n_draws = 0), ranges = data.frame(input = "pair:a:c", low = 0.5, high = 2))
  expect_equal(s$central$total, 2.4364966316833225, tolerance = 1e-9)
  expect_true("pair:a:c" %in% s$inputs$input)
  expect_false(any(c("pair:a:b", "pair:b:c") %in% s$inputs$input))
  ac2 <- deconflate(set_association(m, "a", "c", 2))
  expect_equal(ac2$diagnostics$method, "pairwise shortcut")
  expect_equal(sens_total(s, "pair:a:c", "high"), ac2$totals$adjusted_total, tolerance = 1e-8)
  expect_setequal(s$inputs$input[s$inputs$kind == "entered"],
                  c("prob:a", "prob:b", "prob:c", "associations:1", "associations:2", "impacts:1", "impacts:2",
                    "impacts:3"))
})

test_that("without any association, the baseline is the raw impacts", {
  m0 <- t2_no_assoc()
  expect_error(deconflate(m0), class = "deconflate_unsupported")
  s <- cm_sensitivity(m0, inputs = "unknown",
                      ranges = data.frame(input = c("pair:d1:d2", "pair:d1:d3", "pair:d2:d3"), low = 0.5, high = 2))
  # Everything independent: the aggregate is sum(p * raw) = 2.5.
  expect_equal(s$central$total, 2.5, tolerance = 1e-9)
  # One pair at a time, the others unknown (reference_v040.py, section 7).
  got <- c(sens_total(s, "pair:d1:d2", "low"), sens_total(s, "pair:d1:d2", "high"),
           sens_total(s, "pair:d1:d3", "low"), sens_total(s, "pair:d1:d3", "high"),
           sens_total(s, "pair:d2:d3", "low"), sens_total(s, "pair:d2:d3", "high"))
  expect_equal(got, c(2.558837937037114, 2.4267130832998585, 2.6041868147508, 2.382923828238541,
                      2.698420923951448, 2.2841996245852565), tolerance = 1e-8)
  expect_equal(s$inputs$input[1], "pair:d2:d3")
  # Entered inputs: equal swings per disease (the aggregate is sum(p * raw)).
  se <- cm_sensitivity(m0, inputs = "entered")
  sw <- stats::setNames(se$inputs$swing, se$inputs$input)
  expect_equal(unname(sw[c("impacts:1", "impacts:2", "impacts:3")]), c(0.1, 0.3, 0.6), tolerance = 1e-8)
  expect_equal(unname(sw[c("prob:d1", "prob:d2", "prob:d3")]), c(0.1, 0.3, 0.6), tolerance = 1e-8)
})

test_that("event models: the attributable risk", {
  ev <- t2_event()
  res <- deconflate(ev, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  s <- cm_sensitivity(res, three_way = "all",
                      ranges = data.frame(input = "associations:3", low = 1.5, high = 6))
  expect_equal(s$metric, "attributable risk")
  expect_equal(s$central$total, supp_ar, tolerance = 1e-8)
  expect_equal(c(sens_total(s, "associations:3", "low"), sens_total(s, "associations:3", "high")),
               c(0.04039585852943817, 0.03337529500836167), tolerance = 1e-8)
  # Three-way terms change the event model even with every pair given.
  expect_equal(c(sens_total(s, "or_ratio:d1:d2|d3", "low"), sens_total(s, "or_ratio:d1:d2|d3", "high")),
               c(0.03685486412032779, 0.036447342123939075), tolerance = 1e-8)
  # The overall risk is an input (20% down and up).
  r <- s$inputs[s$inputs$input == "risk", ]
  expect_equal(r$value, 0.25)
  expect_equal(c(r$total_low, r$total_high), c(0.030686662060154746, 0.041864593971753294), tolerance = 1e-8)
  # Interactions as hazard ratios from 0.67 to 1.5.
  expect_equal(unique(s$inputs$range[startsWith(s$inputs$input, "interaction:")]), "HR 0.67-1.5")
  bm <- burden_metric(res)
  expect_equal(sum(bm$by_disease), bm$total, tolerance = 1e-10)
  # A model needs the settings; a result brings them.
  expect_error(cm_sensitivity(ev, keys = "risk"), class = "deconflate_unsupported")
  expect_error(cm_sensitivity(ev, event_model = TRUE, keys = "risk"), "overall_risk")
  expect_error(cm_sensitivity(res, overall_risk = 0.3), "settings come from the result")
})

test_that("event models without associations: the baseline hazard ratios are the raw ones", {
  ev0 <- t2_no_assoc(t2_event())
  expect_error(deconflate(ev0, event_model = TRUE, overall_risk = 0.25), class = "deconflate_unsupported")
  s <- cm_sensitivity(ev0, event_model = TRUE, overall_risk = 0.25, inputs = "unknown",
                      keys = c("pair:d1:d2", "pair:d2:d3"),
                      ranges = data.frame(input = c("pair:d1:d2", "pair:d2:d3"), low = 0.5, high = 2))
  expect_equal(unname(s$central$adjusted), c(1.5, 2.0, 1.3), tolerance = 1e-9)
  expect_equal(s$central$total, 0.04539536116776455, tolerance = 1e-8)
  got <- c(sens_total(s, "pair:d1:d2", "low"), sens_total(s, "pair:d1:d2", "high"),
           sens_total(s, "pair:d2:d3", "low"), sens_total(s, "pair:d2:d3", "high"))
  expect_equal(got, c(0.04781684106024564, 0.04271371201235663, 0.049357621033920274, 0.04136372148403489),
               tolerance = 1e-8)
})
