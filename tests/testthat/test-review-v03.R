# Regression tests for the second external review (deconflate 0.3.0),
# adapted to the 0.4 interface (importance sampling, scenario reweighting
# and the stand-alone Monte Carlo functions were removed in 0.4).
# Reference values: inst/validation/reference_v030.py.

t2_ids <- c("d1", "d2", "d3")
t2_pop <- function(three_way = NULL) {
  cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)),
                three_way = three_way)
}
t2_model <- function(values) cm_model(t2_pop(), cm_impacts(t2_ids, values))
# A model where the published approximation divides by zero for d1:
# b1 = m1^2 / (m1 + A12 m2 + A13 m3) with m = (-A12, 1, 0).
t2_pole_model <- function() {
  e <- deconflate(t2_model(c(1, 1, 1)))$conflation$A["d1", "d2"]
  t2_model(c(-e, 1, 0))
}

# ---- 1. Distribution support ---------------------------------------------------

test_that("a mixture's support is the union of its positive-weight components", {
  gap <- dist_mixture(dist_uniform(0, 0.25), dist_uniform(0.5, 1))
  expect_equal(gap$support, matrix(c(0, 0.5, 0.25, 1), 2))
  expect_equal(c(gap$lower, gap$upper), c(0, 1))
  zero_w <- dist_mixture(dist_uniform(0, 0.5), dist_uniform(0, 1), weights = c(1, 0))
  expect_equal(zero_w$upper, 0.5)
  expect_false(support_covers(gap, dist_uniform(0, 1)))
  expect_false(support_covers(zero_w, dist_uniform(0, 1)))
  expect_true(support_covers(dist_mixture(dist_uniform(0, 1), dist_normal(0.5, 1)), dist_uniform(0, 1)))
  # Nested mixtures: a gap is closed, or left open.
  closed <- dist_mixture(gap, dist_uniform(0.2, 0.6))
  open <- dist_mixture(gap, dist_uniform(0.3, 0.6))
  expect_true(support_covers(closed, dist_uniform(0, 1)))
  expect_false(support_covers(open, dist_uniform(0, 1)))
  expect_equal(format_support(gap), "[0, 0.25] u [0.5, 1]")
})

test_that("the stand-alone Monte Carlo and reweighting functions are gone", {
  ex <- getNamespaceExports("deconflate")
  gone <- c("cm_monte_carlo", "cm_sampler", "cm_batch_sampler", "sampler_global_dairy", "deconflate_hr",
            "cm_hr_model", "cm_hazard_ratios", "attributable_risk", "cm_analyses", "cm_reweight",
            "cm_scenario", "cm_diagnose", "cm_suggest_proposal", "compare_scenarios",
            "example_global_dairy_hr")
  for (f in gone) expect_false(f %in% ex, info = f)
  expect_false("proposal" %in% names(formals(deconflate.cm_model)))
  expect_false("missing_associations" %in% names(formals(cm_population)))
  expect_false("missing_associations" %in% names(formals(cm_model)))
})

# ---- 2. Undefined results are never successful estimates ----------------------

test_that("compare_methods reports an undefined published result as failed", {
  pm <- t2_pole_model()
  expect_true(is.infinite(adjust_impacts(pm, method = "published", warn = FALSE)$adjusted$adjusted[1]))
  cmp <- compare_methods(pm, methods = c("published", "simultaneous"))
  expect_equal(cmp$methods, "simultaneous")
  expect_match(cmp$failed[["published"]], "^undefined")
  expect_s3_class(cmp$undefined$published, "cm_result")
  expect_true(all(is.finite(cmp$impacts$simultaneous)))
  expect_false("published" %in% names(cmp$impacts))
})

test_that("sensitivity tools report undefined scenarios with a reason", {
  # The published approximation is the only way to get a non-finite result,
  # and it is not a method of deconflate(), so the scenario helpers are
  # tested directly on its result.
  pm <- t2_pole_model()
  pub <- adjust_impacts(pm, method = "published", warn = FALSE)
  expect_true(failed_run(pub))
  sm <- scenario_metric(pub)
  expect_null(sm$met)
  expect_match(sm$reason, "^undefined: non-finite adjusted values \\(d1\\)")
  expect_true(is.na(scenario_metric(deconflate(pm, warn = FALSE))$reason))
  base <- burden_metric(deconflate(t2_model(c(1, 1, 1))))
  row <- screen_row(data.frame(scenario = "pole", stringsAsFactors = FALSE), pub, base)
  expect_true(is.na(row$total))
  expect_true(is.na(row$rel_change))
  expect_match(row$failed, "^undefined")
  # A failed run carries its reason.
  sm2 <- scenario_metric(structure(list(), reason = "infeasible"))
  expect_equal(sm2$reason, "infeasible")
  # A screen's baseline must be a finite estimate.
  expect_error(check_baseline(pub), class = "deconflate_nonfinite")
  expect_error(check_baseline(structure(list(), reason = "singular")), "singular")
  # The screens do not accept the published approximation.
  expect_error(screen_associations(pm, method = "published"), class = "deconflate_unsupported")
  expect_error(sensitivity_oat(pm, method = "published"), class = "deconflate_unsupported")
  # With the exact method, the pole model is an ordinary scenario.
  sc <- screen_associations(pm)
  expect_true(all(is.na(sc$failed)))
  expect_true(all(is.finite(sc$total)))
  expect_true(all(is.finite(sensitivity_oat(pm, inputs = "impact")$swing)))
})

# ---- 3. Unavailable precision is not "ok" -------------------------------------

test_that("one usable replicate block gives insufficient_info, not ok", {
  g <- data.frame(x = c(1, 2, 3, 4), w = 1, disease = "q", method = "", block = 1L,
                  raw = NA_real_, adjusted = NA_real_)
  st <- mc_stats(g, c(0.025, 0.5, 0.975), trim = 0.05, check_pole = FALSE, n_blocks = 1L)$row
  expect_true(is.na(st$mcse))
  expect_equal(st$stability, "insufficient_info")
  # Blocks without any accepted draw still count as replicates.
  st4 <- mc_stats(g, c(0.025, 0.5, 0.975), trim = 0.05, check_pole = FALSE, n_blocks = 4L)$row
  expect_true(is.na(st4$mcse))
  expect_equal(st4$stability, "insufficient_info")
})

test_that("Latin hypercube blocks give the ratio-estimator standard error", {
  ratio_se <- function(x, block, R) {
    w <- rep(1 / length(x), length(x))
    blk <- factor(block, levels = seq_len(R))
    W <- tapply(w, blk, sum); W[is.na(W)] <- 0
    S <- tapply(w * x, blk, sum); S[is.na(S)] <- 0
    mu <- sum(w * x)
    sqrt(sum((S - mu * W)^2) / (R * (R - 1))) / (sum(W) / R)
  }
  m <- cm_model(t2_pop(), cm_impacts(t2_ids, c(2.5, 5, 7.5)),
                distributions = list("impact:d1" = dist_normal(2.5, 0.5)))
  res <- deconflate(m, n_draws = 80, sampling = "lhs", lhs_replicates = 4, seed = 5)
  dr <- res$draws
  expect_equal(dr$n_blocks, 4L)
  expect_equal(dr$sampling, "lhs")
  expect_equal(as.vector(table(dr$block)), rep(20L, 4))
  s <- dr$summary[dr$summary$quantity == "adjusted:d1", ]
  expect_equal(s$mcse, ratio_se(dr$values[, "adjusted:d1"], dr$block, 4), tolerance = 1e-10)
  # A block whose draws were all rejected (removed) still counts.
  keep <- dr$block != 3
  g <- data.frame(x = dr$values[keep, "adjusted:d1"], w = 1, disease = "q", method = "",
                  block = dr$block[keep], raw = NA_real_, adjusted = NA_real_)
  st <- mc_stats(g, c(0.025, 0.5, 0.975), trim = 0.05, check_pole = FALSE, n_blocks = 4L)$row
  expect_equal(st$mcse, ratio_se(g$x, g$block, 4), tolerance = 1e-10)
  expect_true(st$stability %in% c("ok", "imprecise", "insufficient_info", "heavy_tail"))
  expect_error(deconflate(m, n_draws = 3, sampling = "lhs"), "at least 4 draws")
})

test_that("a zero mean with positive uncertainty gives insufficient_info", {
  g <- data.frame(x = rep(c(1, -1), 25), w = 1, disease = "a", method = "", block = NA_integer_,
                  raw = NA_real_, adjusted = NA_real_)
  st <- mc_stats(g, c(0.025, 0.5, 0.975), trim = 0.05, check_pole = FALSE)$row
  expect_equal(st$mean, 0)
  expect_gt(st$mcse, 0)
  expect_equal(st$stability, "insufficient_info")
  g$x <- 2 + g$x / 100
  expect_equal(mc_stats(g, c(0.025, 0.5, 0.975), trim = 0.05, check_pole = FALSE)$row$stability, "ok")
})

# ---- 4. Three-way terms and unknown pairs -------------------------------------

test_that("with unknown pairs, a three-way term can change additive results (Bob)", {
  pop <- function(ratio) {
    cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                  three_way = if (ratio != 1) cm_three_way("a", "b", "c", ratio))
  }
  m1 <- cm_model(pop(1), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  m4 <- cm_model(pop(4), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  # Every pair is unknown (no associations): deconflate() refuses, the
  # global engine and the screens do not.
  expect_error(deconflate(m4), class = "deconflate_unsupported")
  expect_equal(global_reasons(m4), c("three-way terms", "3 pairs without an association (unknown)"))
  r1 <- adjust_impacts(m1, method = "global")
  r4 <- adjust_impacts(m4, method = "global")
  expect_equal(r1$totals$adjusted_total, 3, tolerance = 1e-8)
  expect_equal(r1$adjusted$adjusted, c(1, 2, 3), tolerance = 1e-8)
  expect_equal(unname(r4$joint_pairs["a", "b"]), 0.2988400309685928, tolerance = 1e-7)
  expect_equal(r4$totals$adjusted_total, 2.157155621237494, tolerance = 1e-7)
  expect_equal(r4$adjusted$adjusted, c(0.195311748647, 1.438103747492, 2.680895746336),
               tolerance = 1e-7)
  st <- screen_three_way(m1, ratios = 4)
  expect_equal(attr(st, "baseline_total"), 3, tolerance = 1e-8)
  expect_equal(st$total, 2.157155621237494, tolerance = 1e-7)
  expect_equal(st$rel_change, (2.157155621237494 - 3) / 3, tolerance = 1e-7)
})

test_that("with constrained pairs, a three-way term leaves additive results unchanged", {
  m <- cm_model(t2_pop(three_way = cm_three_way("d1", "d2", "d3", 4)),
                cm_impacts(t2_ids, c(2.5, 5, 7.5)))
  res <- deconflate(m)
  expect_equal(res$method, "global")
  expect_equal(res$adjusted$adjusted, deconflate(example_supplement())$adjusted$adjusted, tolerance = 1e-8)
})

# ---- 5. Simulation needs a converged joint ------------------------------------

test_that("simulate_raw_impacts() stops on a joint that did not converge", {
  bad <- cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                       cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
  truth <- c(a = 1, b = 2, c = 3)
  expect_error(simulate_raw_impacts(bad, truth), class = "deconflate_nonconvergence")
  j <- suppressWarnings(fit_joint(bad, max_iter = 50))
  expect_false(j$converged)
  expect_error(simulate_raw_impacts(bad, truth, joint = j), class = "deconflate_nonconvergence")
})

# ---- 7. Event estimands are explicit -------------------------------------------

test_that("event estimands must be stated with their snapshot names", {
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR"), "State the estimand")
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR", estimand = "crude"), "estimand")
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR", estimand = "adjusted_linear",
                          adjusted_for = "all"),
               "estimand")
  h <- cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR", estimand = "snapshot_crude")
  expect_equal(h$estimand, rep("snapshot_crude", 3))
})

test_that("historical conversions are not in the general examples", {
  uk <- example_uk_dairy_2022()
  expect_s3_class(uk, "cm_model")
  expect_equal(attr(uk$impacts, "label"), "milk yield loss")
  expect_null(attr(uk, "hazard_ratios"))
  expect_equal(names(formals(example_uk_dairy_2022)), c("outcome", "yield_sck"))
  # Culling is an event impact table of hazard ratios, not converted.
  ukc <- example_uk_dairy_2022("culling")
  expect_equal(impact_kind(ukc$impacts), "event")
  expect_equal(ukc$impacts$value[ukc$impacts$disease == "DA"], 3.83)
  # The (internal) economic inputs value the yield and fertility analyses.
  eco <- uk_dairy_2022_economics()
  expect_equal(names(eco$valuation), c("yield", "fertility"))
  # Unlisted pairs are unknown: the pairwise methods cannot run.
  cmp <- compare_methods(uk, methods = c("simultaneous", "global"))
  expect_equal(cmp$methods, "global")
  expect_match(cmp$failed[["simultaneous"]], "No association for")
  # Valuation is not exported.
  ex <- getNamespaceExports("deconflate")
  for (f in c("uk_dairy_2022_economics", "productivity_gap", "value_losses", "cm_mc_gap")) {
    expect_false(f %in% ex)
  }
  gd <- example_global_dairy()
  expect_s3_class(gd, "cm_model")
  expect_equal(attr(gd$impacts, "label"), "milk yield loss")
  expect_equal(names(formals(example_global_dairy)), c("outcome", "inputs"))
  gc <- example_global_dairy("culling")
  expect_equal(impact_kind(gc$impacts), "event")
  expect_true(all(paste0("impact:", gc$impacts$disease) %in% names(gc$distributions)))
})
