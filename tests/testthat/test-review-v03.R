# Regression tests for the second external review (deconflate 0.3.0),
# adapted to the 0.4 interface (importance sampling, scenario reweighting
# and the stand-alone Monte Carlo functions were removed in 0.4).
# Reference values: inst/validation/reference_v030.py.

t2_ids <- c("d1", "d2", "d3")
t2_pop <- function(extra = NULL) supp_population(extra)
t2_model <- function(values) cm_model(t2_pop(), cm_impacts(t2_ids, values))

# ---- 1. Distribution support ---------------------------------------------------

test_that("a mixture's support is the union of its positive-weight components", {
  gap <- dist_mixture(dist_uniform(0, 0.25), dist_uniform(0.5, 1))
  expect_equal(gap$support, matrix(c(0, 0.5, 0.25, 1), 2))
  expect_equal(c(gap$lower, gap$upper), c(0, 1))
  zero_w <- dist_mixture(dist_uniform(0, 0.5), dist_uniform(0, 1), weights = c(1, 0))
  expect_equal(zero_w$upper, 0.5)
  expect_equal(format_support(zero_w), "[0, 0.5]")
  # Nested mixtures: a gap is closed, or left open.
  closed <- dist_mixture(gap, dist_uniform(0.2, 0.6))
  open <- dist_mixture(gap, dist_uniform(0.3, 0.6))
  expect_equal(format_support(closed), "[0, 1]")
  expect_equal(format_support(open), "[0, 0.25] u [0.3, 1]")
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
  # Retired in 0.4 (item 4): the three-way and interaction constructors.
  for (f in c("cm_three_way", "cm_interactions")) expect_false(f %in% ex, info = f)
  expect_false("method" %in% names(formals(deconflate.cm_model)))
})

# ---- 2. Undefined results are never successful estimates ----------------------

test_that("cm_sensitivity() needs a finite central result and reports failed runs", {
  bad <- deconflate(t2_model(c(1, 1, 1)), n_draws = 0)
  bad$adjusted$adjusted[1] <- Inf
  expect_error(cm_sensitivity(bad), "not finite", class = "deconflate_nonfinite")
  ok <- deconflate(t2_model(c(1, 1, 1)), n_draws = 0)
  s <- cm_sensitivity(ok, keys = "associations:1", ranges = data.frame(input = "associations:1", low = -1, high = 2))
  expect_match(s$inputs$failed, "^low: infeasible")
  expect_true(is.na(s$inputs$total_low))
  expect_true(is.finite(s$inputs$total_high))
})

# ---- 3. Unavailable precision is not "ok" -------------------------------------

test_that("one usable replicate block gives insufficient_info, not ok", {
  probs <- c(0.025, 0.5, 0.975)
  st <- draw_stats(c(1, 2, 3, 4), rep(1L, 4), probs, n_blocks = 1L)
  expect_true(is.na(st$mcse))
  expect_equal(st$stability, "insufficient_info")
  # Blocks without any accepted draw still count as replicates.
  st4 <- draw_stats(c(1, 2, 3, 4), rep(1L, 4), probs, n_blocks = 4L)
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
  x <- dr$values[keep, "adjusted:d1"]
  st <- draw_stats(x, dr$block[keep], c(0.025, 0.5, 0.975), n_blocks = 4L)
  expect_equal(st$mcse, ratio_se(x, dr$block[keep], 4), tolerance = 1e-10)
  expect_true(st$stability %in% c("ok", "imprecise", "insufficient_info", "heavy_tail"))
  expect_error(deconflate(m, n_draws = 3, sampling = "lhs"), "at least 4 draws")
})

test_that("a zero mean with positive uncertainty gives insufficient_info", {
  probs <- c(0.025, 0.5, 0.975)
  x <- rep(c(1, -1), 25)
  st <- draw_stats(x, NA, probs)
  expect_equal(st$mean, 0)
  expect_gt(st$mcse, 0)
  expect_equal(st$stability, "insufficient_info")
  expect_equal(draw_stats(2 + x / 100, NA, probs)$stability, "ok")
})

# ---- 4. Three-way associations and unknown pairs ------------------------------

test_that("with unknown pairs, a three-way association can change additive results (Bob)", {
  # Three diseases, every pair unknown, and one OR_ratio row: the a-b odds
  # ratio among c is 4 times that among not-c. (With three diseases this is
  # the log-linear three-way term of earlier versions.)
  pop <- function(ratio) {
    cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                  if (ratio != 1) cm_associations("a", "b", ratio, measure = "OR_ratio", among = "c"))
  }
  m1 <- cm_model(pop(1), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  m4 <- cm_model(pop(4), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(deconflate(m1), class = "deconflate_unsupported")
  r1 <- adjust_additive(m1, joint = fit_joint(m1))
  r4 <- deconflate(m4, n_draws = 0)
  expect_equal(r1$totals$adjusted_total, 3, tolerance = 1e-8)
  expect_equal(r1$adjusted$adjusted, c(1, 2, 3), tolerance = 1e-8)
  # Python reference (reference_v041.py, bob_three_way()).
  expect_equal(unname(r4$joint_pairs["a", "b"]), 0.298840030968, tolerance = 1e-7)
  expect_equal(r4$totals$adjusted_total, 2.157155621241, tolerance = 1e-7)
  expect_equal(r4$adjusted$adjusted, c(0.195311748651, 1.438103747494, 2.680895746336),
               tolerance = 1e-7)
  # The OR_ratio row is about a and b; a:c and b:c are unknown.
  expect_equal(paste(r4$unknown_pairs$disease1, r4$unknown_pairs$disease2), c("a c", "b c"))
  st <- cm_sensitivity(m1, inputs = "unknown", keys = "or_ratio:a:b|c",
                       ranges = data.frame(input = "or_ratio:a:b|c", low = 4, high = 5))
  expect_equal(st$central$total, 3, tolerance = 1e-8)
  expect_equal(st$inputs$total_low, 2.157155621241, tolerance = 1e-7)
})

test_that("with crude pair rows, a three-way association leaves additive results unchanged", {
  m <- cm_model(t2_pop(cm_associations("d1", "d2", 4, measure = "OR_ratio", among = "d3")),
                cm_impacts(t2_ids, c(2.5, 5, 7.5)))
  res <- deconflate(m)
  expect_equal(res$method, "full model")
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

# ---- 7. Event estimands follow from adjusted_for --------------------------------

test_that("event impacts need a measure; adjusted_for decides crude or stratified", {
  h <- cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR")
  expect_equal(attr(h, "kind"), "event")
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR", estimand = "snapshot_crude"),
               "no longer used", class = "deconflate_unsupported")
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR2"), "measure")
})

test_that("the general examples keep the impacts on their own scale", {
  # Valuation is not exported.
  ex <- getNamespaceExports("deconflate")
  for (f in c("productivity_gap", "value_losses", "cm_mc_gap")) {
    expect_false(f %in% ex)
  }
  gd <- example_global_dairy()
  expect_s3_class(gd, "cm_model")
  expect_equal(attr(gd$impacts, "label"), "milk yield loss")
  expect_null(attr(gd, "hazard_ratios"))
  expect_equal(names(formals(example_global_dairy)), c("outcome", "inputs"))
  # Culling is an event impact table of hazard ratios, not converted.
  gc <- example_global_dairy("culling")
  expect_equal(impact_kind(gc$impacts), "event")
  expect_true(all(paste0("impacts:", gc$impacts$row) %in% names(gc$distributions)))
})
