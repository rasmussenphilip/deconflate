# Regression tests for the second external review (deconflate 0.3.0).
# Reference values: inst/validation/reference_v030.py.

# ---- 1. Importance-sampling support -------------------------------------------

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
})

test_that("importance sampling rejects proposals with holes in their support", {
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d1 = dist_uniform(0, 1)))
  gap <- dist_mixture(dist_uniform(0, 0.25), dist_uniform(0.5, 1))
  zero_w <- dist_mixture(dist_uniform(0, 0.5), dist_uniform(0, 1), weights = c(1, 0))
  expect_error(cm_monte_carlo(s, 10, proposal = list("impact:d1" = gap), seed = 1),
               class = "deconflate_unsupported")
  expect_error(cm_monte_carlo(s, 10, proposal = list("impact:d1" = zero_w), seed = 1),
               class = "deconflate_unsupported")
  ok <- dist_mixture(dist_uniform(0, 1), dist_uniform(0.4, 0.6), weights = c(0.5, 0.5))
  mc <- cm_monte_carlo(s, 50, proposal = list("impact:d1" = ok), seed = 1)
  expect_s3_class(mc, "cm_mc")
})

test_that("scenario reweighting checks the support of the scenario", {
  # Sampled from a distribution with a gap: a scenario that puts weight in
  # the gap cannot be represented.
  gap <- dist_mixture(dist_uniform(0, 0.25), dist_uniform(0.5, 1))
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d1 = gap))
  mc <- cm_monte_carlo(s, 40, seed = 2)
  expect_error(cm_scenario(mc, list("impact:d1" = dist_uniform(0, 1))),
               class = "deconflate_unsupported")
  expect_s3_class(cm_scenario(mc, list("impact:d1" = dist_uniform(0.6, 0.9))), "cm_mc")
})

# ---- 2. Undefined results are never successful estimates ----------------------

test_that("compare_methods reports an undefined published result as failed", {
  pm <- pole_model()
  expect_true(is.infinite(deconflate(pm, method = "published", warn = FALSE)$adjusted$adjusted[1]))
  cmp <- compare_methods(pm, methods = c("published", "simultaneous"))
  expect_equal(cmp$methods, "simultaneous")
  expect_match(cmp$failed[["published"]], "^undefined")
  expect_s3_class(cmp$undefined$published, "cm_result")
  expect_true(all(is.finite(cmp$impacts$simultaneous)))
  expect_false("published" %in% names(cmp$impacts))
})

test_that("sensitivity tools report undefined scenarios with a reason", {
  pm <- pole_model()
  cs <- compare_scenarios(base = supp_model(c(1, 1, 1)), pole = pm, method = "published")
  expect_true(is.na(cs$totals$failed[1]))
  expect_true(is.na(cs$totals$total[2]))
  expect_match(cs$totals$failed[2], "^undefined")
  expect_error(screen_associations(pm, method = "published"), class = "deconflate_nonfinite")
  expect_true(failed_run(deconflate(pm, method = "published", warn = FALSE)))
})

# ---- 3. Unavailable precision is not "ok" -------------------------------------

lhs_run <- function(n = 80, R = 4, seed = 5) {
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d1 = dist_normal(2.5, 0.5)))
  cm_monte_carlo(s, n, sampling = "lhs", lhs_replicates = R, seed = seed)
}

test_that("one usable replicate block gives insufficient_info, not ok", {
  mc <- lhs_run()
  only1 <- cm_reweight(mc, function(p) ifelse(mc$block == 1, 0, -Inf))
  sm <- summary(only1, diagnose = FALSE)
  expect_true(all(is.na(sm$mcse)))
  expect_true(all(sm$stability == "insufficient_info"))
  expect_true(all(cm_diagnose(only1)$stability == "insufficient_info"))
})

test_that("blocks with zero weight or no accepted draws count as replicates", {
  mc <- lhs_run()
  expect_equal(mc$n_blocks, 4)
  ratio_se <- function(mc, dis) {
    d <- mc$draws[mc$draws$disease == dis & mc$draws$method == "simultaneous", ]
    idx <- match(d$draw, mc$params$draw)
    w <- mc$weights[idx] / sum(mc$weights[idx])
    blk <- factor(mc$block[idx], levels = seq_len(mc$n_blocks))
    W <- tapply(w, blk, sum); W[is.na(W)] <- 0
    S <- tapply(w * d$adjusted, blk, sum); S[is.na(S)] <- 0
    mu <- sum(w * d$adjusted)
    R <- mc$n_blocks
    sqrt(sum((S - mu * W)^2) / (R * (R - 1))) / (sum(W) / R)
  }
  # A block with zero weight.
  zero <- cm_reweight(mc, function(p) ifelse(mc$block == 2, -Inf, 0))
  sz <- summary(zero, diagnose = FALSE)
  expect_equal(sz$mcse[sz$disease == "d1"], ratio_se(zero, "d1"), tolerance = 1e-10)
  # A block whose draws were all rejected (removed from the run).
  keep <- mc$block != 3
  gone <- mc
  gone$params <- mc$params[keep, , drop = FALSE]
  gone$block <- mc$block[keep]
  gone$weights <- mc$weights[keep] / sum(mc$weights[keep])
  gone$log_weights <- mc$log_weights[keep]
  gone$draws <- mc$draws[mc$draws$draw %in% gone$params$draw, , drop = FALSE]
  gone$totals <- mc$totals[mc$totals$draw %in% gone$params$draw, , drop = FALSE]
  sg <- summary(gone, diagnose = FALSE)
  expect_equal(sg$mcse[sg$disease == "d1"], ratio_se(gone, "d1"), tolerance = 1e-10)
  expect_true(all(sg$stability %in% c("ok", "imprecise", "insufficient_info", "heavy_tail")))
})

test_that("a zero mean with positive uncertainty gives insufficient_info", {
  pop <- cm_population(cm_diseases(c("a", "b"), c(0.2, 0.3)))
  s <- function(i) cm_model(pop, cm_impacts(c("a", "b"), c(if (i %% 2) 1 else -1, 2)))
  mc <- cm_monte_carlo(s, 50, seed = 1)
  sm <- summary(mc, diagnose = FALSE)
  a <- sm[sm$disease == "a", ]
  expect_equal(a$mean, 0)
  expect_gt(a$mcse, 0)
  expect_equal(a$stability, "insufficient_info")
  expect_equal(sm$stability[sm$disease == "b"], "ok")
})

# ---- 4. Three-way terms and unknown pairs -------------------------------------

test_that("with unknown pairs, a three-way term can change additive results (Bob)", {
  pop <- function(ratio) {
    cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                  three_way = if (ratio != 1) cm_three_way("a", "b", "c", ratio),
                  missing_associations = "unknown")
  }
  m1 <- cm_model(pop(1), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  m4 <- cm_model(pop(4), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  r1 <- deconflate(m1, method = "global")
  r4 <- deconflate(m4, method = "global")
  expect_equal(r1$totals$adjusted_total, 3, tolerance = 1e-8)
  expect_equal(r1$adjusted$adjusted, c(1, 2, 3), tolerance = 1e-8)
  expect_equal(unname(r4$joint_pairs["a", "b"]), 0.2988400309685928, tolerance = 1e-7)
  expect_equal(r4$totals$adjusted_total, 2.157155621237494, tolerance = 1e-7)
  expect_equal(r4$adjusted$adjusted, c(0.195311748647, 1.438103747492, 2.680895746336),
               tolerance = 1e-7)
})

test_that("with constrained pairs, a three-way term leaves additive results unchanged", {
  m <- cm_model(supp_population(three_way = cm_three_way("d1", "d2", "d3", 4)),
                cm_impacts(ids3, c(2.5, 5, 7.5)))
  expect_equal(deconflate(m, method = "global")$adjusted$adjusted,
               deconflate(example_supplement())$adjusted$adjusted, tolerance = 1e-8)
})

# ---- 5. Simulation needs a converged joint ------------------------------------

test_that("simulate_raw_impacts() stops on a joint that did not converge", {
  bad <- bad_population()
  truth <- c(a = 1, b = 2, c = 3)
  expect_error(simulate_raw_impacts(bad, truth), class = "deconflate_nonconvergence")
  j <- suppressWarnings(fit_joint(bad, max_iter = 50))
  expect_false(j$converged)
  expect_error(simulate_raw_impacts(bad, truth, joint = j), class = "deconflate_nonconvergence")
})

# ---- 7. Hazard-ratio estimands are explicit -----------------------------------

test_that("hazard-ratio estimands must be stated with their snapshot names", {
  expect_error(cm_hazard_ratios(ids3, c(1.5, 2, 1.3)), "State the estimand")
  expect_error(cm_hazard_ratios(ids3, c(1.5, 2, 1.3), estimand = "crude"),
               class = "deconflate_unsupported")
  expect_error(cm_hazard_ratios(ids3, c(1.5, 2, 1.3), estimand = "adjusted", adjusted_for = "all"),
               class = "deconflate_unsupported")
  h <- cm_hazard_ratios(ids3, c(1.5, 2, 1.3), estimand = "snapshot_crude")
  expect_equal(h$estimand, rep("snapshot_crude", 3))
})

test_that("historical conversions are not in the general examples", {
  uk <- example_uk_dairy_2022()
  expect_equal(names(uk$models), c("yield", "fertility"))
  expect_equal(attr(uk, "hazard_ratios")$estimand[1], "snapshot_crude")
  # The economic inputs value exactly the analyses of the example.
  eco <- uk_dairy_2022_economics()
  expect_equal(names(eco$valuation), names(uk$models))
  expect_equal(compare_methods(uk, methods = "simultaneous", valuation = eco$valuation)$methods,
               "simultaneous")
  gd <- example_global_dairy()
  expect_equal(names(gd$models), c("yield", "fertility"))
  expect_equal(names(sampler_global_dairy()$samplers), c("yield", "fertility"))
  expect_false("culling" %in% names(formals(example_global_dairy)))
})
