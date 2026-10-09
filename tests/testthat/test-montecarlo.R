# Monte Carlo propagation (R/sampler.R, R/montecarlo.R).

# A smooth sampler on the supplement model: two uncertain impacts.
mc_smooth_sampler <- function() {
  cm_sampler(supp_model(c(2.5, 5, 7.5)),
             impacts = list(d1 = dist_normal(2.5, 0.5), d3 = dist_normal(7.5, 1)))
}

# Rows of a summary for one disease (or quantity) and method.
mc_row <- function(s, item, method) {
  key <- if (!is.null(s$quantity)) s$quantity else s$disease
  s[key == item & s$method == method, , drop = FALSE]
}

# Draws of one disease and method, with their normalised weights and blocks.
mc_group <- function(mc, disease, method) {
  d <- mc$draws[mc$draws$disease == disease & mc$draws$method == method, , drop = FALSE]
  idx <- match(d$draw, mc$params$draw)
  d$w <- mc$weights[idx] / sum(mc$weights[idx])
  if (!is.null(mc$block)) d$block <- mc$block[idx]
  d
}

stability_levels <- c("ok", "imprecise", "insufficient_info", "heavy_tail", "possible_pole")

# ---- cm_sampler -------------------------------------------------------------

test_that("cm_sampler maps prob:, assoc: (either order), impact: and inter: keys", {
  m <- supp_model(c(2.5, 5, 7.5), interactions = cm_interactions("d1", "d2", 0.5))
  s <- cm_sampler(m,
                  diseases = list(d1 = dist_fixed(0.12)),
                  associations = list("d2:d1" = dist_fixed(2.5), "d3:d2" = dist_fixed(1.5)),
                  impacts = list(d3 = dist_fixed(8)),
                  interactions = list("d2:d1" = dist_fixed(1)))
  expect_s3_class(s, "cm_sampler")
  specs <- attr(s, "specs")
  expect_setequal(names(specs),
                  c("prob:d1", "assoc:d1:d2", "assoc:d2:d3", "impact:d3", "inter:d1:d2"))
  expect_true(all(vapply(specs, inherits, logical(1), "cm_dist")))
  expect_equal(specs[["impact:d3"]]$mean, 8)

  d <- s(1)
  expect_s3_class(d, "cm_model")
  expect_equal(d$diseases$value, c(0.12, 0.15, 0.20))
  expect_equal(d$diseases$prob, c(0.12, 0.15, 0.20))
  # Association rows d1:d2, d1:d3, d2:d3; d1:d3 is not varied.
  expect_equal(d$associations$value, c(2.5, 1, 1.5))
  expect_equal(d$impacts$value, c(2.5, 5, 8))
  expect_equal(d$interactions$value, 1)
  # The original model is not modified.
  expect_equal(m$impacts$value, c(2.5, 5, 7.5))
  expect_output(print(s), "<cm_sampler> 5 uncertain inputs")

  # Incidence rates are drawn on their own scale and converted.
  pop_ir <- cm_population(
    cm_diseases(ids3, c(0.10, 0.15, 0.20), type = "incidence_rate"),
    cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)))
  s_ir <- cm_sampler(cm_model(pop_ir, cm_impacts(ids3, c(2.5, 5, 7.5))),
                     diseases = list(d2 = dist_fixed(0.5)))
  d_ir <- s_ir(1)
  expect_equal(d_ir$diseases$value[2], 0.5)
  expect_equal(d_ir$diseases$prob[2], 1 - exp(-0.5))
})

test_that("cm_sampler accepts stratified uniforms and fixed values", {
  s <- mc_smooth_sampler()
  expect_equal(s(1, u = c("impact:d1" = 0.5))$impacts$value[1], 2.5)
  d <- s(1, values = c("impact:d1" = 3, "impact:d3" = 6))
  expect_equal(d$impacts$value, c(3, 5, 6))
})

test_that("cm_sampler rejects keys that do not match the model", {
  m <- supp_model(c(2.5, 5, 7.5))
  expect_error(cm_sampler(supp_population()), "cm_model")
  expect_error(cm_sampler(m, impacts = list(d1 = 3)), "named list of cm_dist")
  expect_error(cm_sampler(m, impacts = list(dist_fixed(3))), "named list of cm_dist")
  expect_error(cm_sampler(m, diseases = list(zz = dist_fixed(0.1))), "Unknown disease")
  expect_error(cm_sampler(m, associations = list("d1-d2" = dist_fixed(2))), "does not match")
  expect_error(cm_sampler(m, associations = list("d1:d2:d3" = dist_fixed(2))), "does not match")
  expect_error(cm_sampler(m, associations = list("d1:d2" = dist_fixed(2), "d2:d1" = dist_fixed(3))),
               "given twice")
  expect_error(cm_sampler(m, impacts = list(zz = dist_fixed(1))), "does not match an impact")
  expect_error(cm_sampler(m, impacts = list(d1 = dist_fixed(1), d1 = dist_fixed(2))),
               "more than once")
  expect_error(cm_sampler(m, interactions = list("d1:d2" = dist_fixed(1))), "no interactions")
  mi <- supp_model(c(2.5, 5, 7.5), interactions = cm_interactions("d1", "d2", 0.5))
  expect_error(cm_sampler(mi, interactions = list("d1:d3" = dist_fixed(1))), "does not match an interaction")
  m_noassoc <- cm_model(cm_population(cm_diseases(ids3, c(0.1, 0.15, 0.2))),
                        cm_impacts(ids3, c(2.5, 5, 7.5)))
  expect_error(cm_sampler(m_noassoc, associations = list("d1:d2" = dist_fixed(2))),
               "no associations")
  m_indep <- cm_model(
    cm_population(cm_diseases(ids3, c(0.1, 0.15, 0.2)),
                  cm_associations(c("d1", "d2"), c("d2", "d3"), c(NA, 3),
                                  measure = c("independent", "OR"))),
    cm_impacts(ids3, c(2.5, 5, 7.5)))
  expect_error(cm_sampler(m_indep, associations = list("d1:d2" = dist_fixed(2))),
               "numeric measure")
})

test_that("impossible sampled inputs raise deconflate_infeasible", {
  m <- supp_model(c(2.5, 5, 7.5))
  s <- cm_sampler(m, diseases = list(d1 = dist_normal(0.1, 0.05)),
                  associations = list("d1:d2" = dist_normal(2, 1)),
                  impacts = list(d1 = dist_normal(2.5, 0.5)))
  v <- c("prob:d1" = 0.1, "assoc:d1:d2" = 2, "impact:d1" = 2.5)
  expect_s3_class(s(1, values = v), "cm_model")
  bad <- function(k, x) { v[[k]] <- x; v }
  expect_error(s(1, values = bad("prob:d1", -0.01)), "invalid probability",
               class = "deconflate_infeasible")
  expect_error(s(1, values = bad("prob:d1", 1.2)), "invalid probability",
               class = "deconflate_infeasible")
  expect_error(s(1, values = bad("assoc:d1:d2", -0.5)), "not a valid OR",
               class = "deconflate_infeasible")
  expect_error(s(1, values = bad("impact:d1", Inf)), "not finite",
               class = "deconflate_infeasible")
})

test_that("three-way ratios can be sampled, in any order of the diseases", {
  pop <- supp_population(three_way = cm_three_way("d2", "d1", "d3", 2))
  m <- cm_model(pop, cm_impacts(ids3, c(2.5, 5, 7.5)), cm_interactions("d1", "d2", 0.5))
  s <- cm_sampler(m, three_way = list("d3:d1:d2" = dist_lognormal_ci(2, 1.2, 3.5)))
  expect_equal(names(attr(s, "specs")), "three:d2:d1:d3")
  expect_equal(s(1, values = c("three:d2:d1:d3" = 1.7))$three_way$ratio, 1.7)
  expect_error(s(1, values = c("three:d2:d1:d3" = 0)), "positive", class = "deconflate_infeasible")
  expect_error(s(1, values = c("three:d2:d1:d3" = -1)), class = "deconflate_infeasible")
  expect_error(cm_sampler(m, three_way = list("d1:d2" = dist_fixed(2))), "does not match")
  expect_error(cm_sampler(m, three_way = list("d1:d2:d4" = dist_fixed(2))), "does not match")
  expect_error(cm_sampler(m, three_way = list("d1:d2:d3" = dist_fixed(2), "d3:d2:d1" = dist_fixed(3))),
               "given twice")
  expect_error(cm_sampler(supp_model(c(2.5, 5, 7.5)), three_way = list("d1:d2:d3" = dist_fixed(2))),
               "no three-way terms")
  mc <- cm_monte_carlo(s, 20, method = "global", seed = 1)
  expect_equal(mc$n_rejected, 0)
  expect_gt(length(unique(mc$params[["three:d2:d1:d3"]])), 1)
  # The ratio changes the global result when interactions are present.
  lo <- deconflate(set_three_way(m, "d1", "d2", "d3", 0.5), method = "global")$totals$adjusted_total
  hi <- deconflate(set_three_way(m, "d1", "d2", "d3", 4), method = "global")$totals$adjusted_total
  expect_false(isTRUE(all.equal(lo, hi)))

  # In a batch, the three-way ratio is a shared population input.
  an <- cm_analyses(pop, a = cm_impacts(ids3, c(2.5, 5, 7.5)), b = cm_impacts(ids3, c(1, 2, 3)),
                    interactions = list(a = cm_interactions("d1", "d2", 0.5)))
  bs <- cm_batch_sampler(an, three_way = list("d1:d2:d3" = dist_uniform(1, 3)),
                         impacts = list(b = list(d1 = dist_normal(1, 0.2))))
  expect_equal(bs$population_keys, "three:d2:d1:d3")
  mcb <- cm_monte_carlo(bs, 15, method = "global", seed = 2)
  pa <- mcb$analyses$a$params
  pb <- mcb$analyses$b$params
  expect_equal(pa[["three:d2:d1:d3"]], pb[["three:d2:d1:d3"]][match(pa$draw, pb$draw)])
})

# ---- Rejections -------------------------------------------------------------

test_that("invalid draws are rejected, counted and tabulated", {
  m <- supp_model(c(2.5, 5, 7.5))
  # A normal prevalence and a normal odds ratio both go negative in ~5% of draws.
  s <- cm_sampler(m, diseases = list(d1 = dist_normal(0.1, 0.06)),
                  associations = list("d1:d2" = dist_normal(2, 1.2)),
                  impacts = list(d3 = dist_normal(7.5, 1)))
  mc <- cm_monte_carlo(s, 300, seed = 5)
  expect_s3_class(mc, "cm_mc")
  expect_gt(mc$n_rejected, 0)
  expect_equal(mc$n_rejected + nrow(mc$params), 300)
  expect_equal(nrow(mc$rejections), mc$n_rejected)
  expect_named(mc$rejections, c("draw", "type", "reason"))
  expect_true(all(mc$rejections$type == "infeasible"))
  expect_true(any(grepl("invalid probability", mc$rejections$reason)))
  expect_true(any(grepl("is not a valid OR", mc$rejections$reason)))
  # Rejected draws leave no trace in the accepted results.
  expect_length(intersect(mc$rejections$draw, mc$params$draw), 0)
  expect_setequal(c(mc$rejections$draw, mc$params$draw), 1:300)
  expect_false(any(mc$draws$draw %in% mc$rejections$draw))
  expect_false(any(mc$totals$draw %in% mc$rejections$draw))
  expect_length(mc$weights, nrow(mc$params))
  # Every accepted draw is valid.
  expect_true(all(mc$params[["prob:d1"]] > 0))
  expect_true(all(mc$params[["assoc:d1:d2"]] > 0))
  expect_true(all(c("prob:d1", "prob:d2", "prob:d3", "assoc:d1:d2", "assoc:d1:d3",
                    "assoc:d2:d3", "impact:d1", "impact:d2", "impact:d3") %in% names(mc$params)))
  expect_equal(unique(mc$params[["impact:d2"]]), 5)

  rt <- summary(mc, what = "rejections")
  expect_named(rt, c("type", "reason", "n"))
  expect_equal(sum(rt$n), mc$n_rejected)
  expect_true(all(rt$type == "infeasible"))
  expect_false(is.unsorted(rev(rt$n)))
  expect_output(print(mc), sprintf("rejected: %d", mc$n_rejected))

  # No rejections: an empty table with the same columns.
  mc0 <- cm_monte_carlo(mc_smooth_sampler(), 20, seed = 1)
  expect_equal(mc0$n_rejected, 0)
  rt0 <- summary(mc0, what = "rejections")
  expect_equal(nrow(rt0), 0)
  expect_named(rt0, c("type", "reason", "n"))
})

test_that("unsupported combinations are rejected by type", {
  m <- supp_model(c(2.5, 5, 7.5), interactions = cm_interactions("d1", "d2", 0.5))
  s <- cm_sampler(m, interactions = list("d2:d1" = dist_normal(0.5, 0.2)))
  # Interactions need the global method: every draw is rejected.
  mc <- cm_monte_carlo(s, 10, method = "simultaneous", seed = 1)
  expect_equal(mc$n_rejected, 10)
  expect_null(mc$draws)
  expect_error(summary(mc), "All draws were rejected")
  rt <- summary(mc, what = "rejections")
  expect_equal(nrow(rt), 1)
  expect_equal(rt$type, "unsupported")
  expect_equal(rt$n, 10)
  # With the global method the interaction is propagated.
  mg <- cm_monte_carlo(s, 10, method = "global", seed = 1)
  expect_equal(mg$n_rejected, 0)
  expect_gt(length(unique(mg$params[["inter:d1:d2"]])), 1)
})

test_that("a non-finite result rejects the whole draw for every method", {
  # Published method: b1 = m1^2 / (m1 + A12 m2 + A13 m3). With m2 = 1 and
  # m3 = 0 the denominator is -A12 + A12 = 0 exactly when m1 = -A12.
  e <- deconflate(supp_model(c(1, 1, 1)))$conflation$A["d1", "d2"]
  expect_gt(e, 0)
  pole <- supp_model(c(-e, 1, 0))
  expect_true(is.infinite(adjust_impacts(pole, method = "published", warn = FALSE)$adjusted$adjusted[1]))
  expect_true(all(is.finite(deconflate(pole, method = "simultaneous", warn = FALSE)$adjusted$adjusted)))

  good <- supp_model(c(2.5, 5, 7.5))
  sampler <- function(i) if (i == 3) pole else good
  mc <- cm_monte_carlo(sampler, 5, method = c("simultaneous", "published"))
  expect_equal(mc$n_rejected, 1)
  expect_equal(mc$rejections$draw, 3)
  expect_equal(mc$rejections$type, "nonfinite")
  expect_match(mc$rejections$reason, "published")
  # Draw 3 is dropped for both methods, including the finite simultaneous one.
  expect_setequal(unique(mc$draws$draw), c(1, 2, 4, 5))
  expect_equal(nrow(mc$draws), 4 * 2 * 3)
  expect_false(3 %in% mc$totals$draw)
  expect_false(3 %in% mc$params$draw)
  expect_true(all(is.finite(mc$draws$adjusted)))
  diag <- attr(summary(mc, diagnose = FALSE), "diagnosis")
  nf <- diag[diag$stability == "non_finite", , drop = FALSE]
  expect_equal(nrow(nf), 1)
  expect_equal(nf$item, "all")
  expect_equal(nf$method, "simultaneous, published")
  expect_match(nf$detail, "1 draw")
  expect_true("non_finite" %in% cm_diagnose(mc)$stability)

  # The simultaneous method alone accepts every draw.
  expect_equal(cm_monte_carlo(sampler, 5, method = "simultaneous")$n_rejected, 0)

  # A sampled value can also make every draw non-finite.
  s <- cm_sampler(supp_model(c(2.5, 1, 0)), impacts = list(d1 = dist_fixed(-e)))
  mc_all <- cm_monte_carlo(s, 4, method = "published", seed = 1)
  expect_equal(mc_all$n_rejected, 4)
  expect_true(all(mc_all$rejections$type == "nonfinite"))
  expect_null(mc_all$draws)
  expect_equal(mc_all$ess, 0)
  expect_error(summary(mc_all), "All draws were rejected")
  expect_equal(summary(mc_all, what = "rejections")$n, 4)
  expect_output(print(mc_all), "rejected: 4")
})

# ---- Summaries ----------------------------------------------------------------

test_that("summaries report means, MCSE, quantiles and consistent totals", {
  s <- mc_smooth_sampler()
  mc <- cm_monte_carlo(s, 100, method = c("published", "simultaneous"), seed = 3)
  expect_equal(mc$method, c("published", "simultaneous"))
  expect_equal(mc$sampling, "random")
  expect_null(mc$block)
  expect_equal(mc$weights, rep(1 / 100, 100))
  expect_equal(mc$ess, 100)
  sm <- summary(mc, diagnose = FALSE)
  expect_equal(nrow(sm), 6)
  expect_true(all(c("disease", "method", "mean", "sd", "mcse", "q0.025", "q0.5", "q0.975",
                    "trimmed_mean", "rel_mcse", "tail_share", "stability") %in% names(sm)))
  expect_true(all(sm$stability %in% stability_levels))
  expect_equal(attr(sm, "n_draws"), 100)
  expect_equal(attr(sm, "n_rejected"), 0)

  # Equal weights: mean, SD (divisor n) and MCSE = sd / sqrt(n).
  g <- mc_group(mc, "d1", "simultaneous")
  r <- mc_row(sm, "d1", "simultaneous")
  mu <- mean(g$adjusted)
  expect_equal(r$mean, mu)
  expect_equal(r$sd, sqrt(mean((g$adjusted - mu)^2)))
  expect_equal(r$mcse, sqrt(sum((g$adjusted - mu)^2)) / 100)
  expect_equal(r$rel_mcse, r$mcse / abs(mu))
  expect_true(r$`q0.025` <= r$`q0.5` && r$`q0.5` <= r$`q0.975`)

  # Other quantiles; no trimming gives the mean.
  sq <- summary(mc, probs = c(0.1, 0.9), trim = 0, diagnose = FALSE)
  expect_true(all(c("q0.1", "q0.9") %in% names(sq)))
  expect_equal(sq$trimmed_mean, sq$mean)

  # Contributions add up to the adjusted total, per method.
  sc <- summary(mc, what = "contribution", diagnose = FALSE)
  st <- summary(mc, what = "total", diagnose = FALSE)
  expect_true("quantity" %in% names(st))
  expect_setequal(st$quantity, c("adjusted_total", "raw_sum"))
  for (m in mc$method) {
    expect_equal(sum(sc$mean[sc$method == m]), mc_row(st, "adjusted_total", m)$mean)
  }
  expect_equal(mc_row(st, "raw_sum", "raw")$mean,
               mean(mc$totals$raw_sum[mc$totals$method == "published"]))
})

test_that("progress messages are printed on request", {
  msgs <- capture_messages(cm_monte_carlo(mc_smooth_sampler(), 10, progress = TRUE, seed = 1))
  expect_length(msgs, 10)
  expect_true(any(grepl("Draw 10 of 10", msgs)))
})

test_that("cm_monte_carlo checks its sampler", {
  expect_error(cm_monte_carlo(supp_model(c(2.5, 5, 7.5)), 5), "must be a function")
})

# ---- Importance sampling ----------------------------------------------------

test_that("proposals must cover the support of the input distribution", {
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)),
                  impacts = list(d1 = dist_normal(2.5, 0.5), d2 = dist_fixed(5),
                                 d3 = dist_pert(6, 7.5, 9)))
  ok <- function(p) cm_monte_carlo(s, 3, proposal = p, seed = 1)
  # Bounded proposal for an unbounded input.
  expect_error(ok(list("impact:d1" = dist_uniform(1, 4))), "does not cover",
               class = "deconflate_unsupported")
  expect_error(ok(list("impact:d1" = dist_normal(2.5, 1, lower = 0))), "does not cover",
               class = "deconflate_unsupported")
  # Lower and upper bounds of a bounded input.
  expect_error(ok(list("impact:d3" = dist_uniform(6.5, 10))), "does not cover",
               class = "deconflate_unsupported")
  expect_error(ok(list("impact:d3" = dist_uniform(5, 8.5))), "does not cover",
               class = "deconflate_unsupported")
  # Point masses cannot be importance-sampled, in either role.
  expect_error(ok(list("impact:d2" = dist_normal(5, 1))), "point masses",
               class = "deconflate_unsupported")
  expect_error(ok(list("impact:d3" = dist_fixed(7.5))), "point masses",
               class = "deconflate_unsupported")
  expect_error(ok(list("impact:d1" = dist_mixture(dist_fixed(2.5), dist_normal(2.5, 1)))),
               "point masses", class = "deconflate_unsupported")
  # Unknown keys and malformed proposals.
  expect_error(ok(list("impact:zz" = dist_normal(1, 1))), "No sampled input")
  expect_error(ok(list(dist_normal(1, 1))), "named list")
  expect_error(cm_monte_carlo(function(i) supp_model(c(2.5, 5, 7.5)), 3,
                              proposal = list("impact:d1" = dist_normal(2.5, 1))),
               "cm_sampler")
  # Covering proposals are accepted.
  expect_s3_class(ok(list("impact:d1" = dist_normal(2.5, 1), "impact:d3" = dist_uniform(6, 9))),
                  "cm_mc")
  expect_s3_class(ok(list("impact:d3" = dist_pert(6, 7.5, 9))), "cm_mc")
})

test_that("importance sampling gives weights, ESS and the self-normalised MCSE", {
  s <- mc_smooth_sampler()
  prop <- list("impact:d1" = dist_normal(2.5, 1))
  mc <- cm_monte_carlo(s, 200, proposal = prop, seed = 11)
  expect_equal(mc$n_rejected, 0)
  expect_identical(mc$proposal, prop)
  x <- mc$params[["impact:d1"]]
  lw <- stats::dnorm(x, 2.5, 0.5, log = TRUE) - stats::dnorm(x, 2.5, 1, log = TRUE)
  expect_equal(mc$log_weights, lw)
  w <- exp(lw - max(lw))
  expect_equal(mc$weights, w / sum(w))
  expect_equal(sum(mc$weights), 1)
  expect_equal(mc$ess, 1 / sum(mc$weights^2))
  expect_lt(mc$ess, 200)
  expect_gt(mc$ess, 50)
  expect_output(print(mc), "importance sampling of impact:d1")

  sm <- summary(mc, diagnose = FALSE)
  for (dis in ids3) {
    g <- mc_group(mc, dis, "simultaneous")
    mu <- sum(g$w * g$adjusted)
    r <- mc_row(sm, dis, "simultaneous")
    expect_equal(r$mean, mu)
    expect_equal(r$sd, sqrt(sum(g$w * (g$adjusted - mu)^2)))
    expect_equal(r$mcse, sqrt(sum(g$w^2 * (g$adjusted - mu)^2)))
  }
  # The weighted raw impact recovers the input's own mean (2.5).
  g <- mc_group(mc, "d1", "simultaneous")
  expect_lt(abs(sum(g$w * g$raw) - 2.5), 0.2)

  # Draws outside the input's support get zero weight.
  sb <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d3 = dist_pert(6, 7.5, 9)))
  mb <- cm_monte_carlo(sb, 60, proposal = list("impact:d3" = dist_uniform(5, 10)), seed = 3)
  x3 <- mb$params[["impact:d3"]]
  expect_true(any(x3 < 6 | x3 > 9))
  expect_true(all(mb$weights[x3 < 6 | x3 > 9] == 0))
  expect_true(all(mb$weights[x3 > 6 & x3 < 9] > 0))
})

# ---- Latin hypercube sampling -------------------------------------------------

test_that("LHS uses replicate blocks and block-mean standard errors", {
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)),
                  impacts = list(d1 = dist_uniform(2, 3), d3 = dist_normal(7.5, 1)))
  mc <- cm_monte_carlo(s, 100, sampling = "lhs", lhs_replicates = 5, seed = 7)
  expect_equal(mc$sampling, "lhs")
  expect_equal(mc$n_rejected, 0)
  expect_length(mc$block, nrow(mc$params))
  expect_equal(sort(unique(mc$block)), 1:5)
  expect_equal(as.vector(table(mc$block)), rep(20, 5))
  expect_output(print(mc), "Latin hypercube, 5 replicate blocks")
  # Each block is stratified: one draw in each of 20 equal-probability strata.
  for (b in 1:5) {
    sel <- mc$block == b
    expect_equal(sort(floor((mc$params[["impact:d1"]][sel] - 2) * 20)), 0:19)
    expect_equal(sort(floor(stats::pnorm(mc$params[["impact:d3"]][sel], 7.5, 1) * 20)), 0:19)
  }
  sm <- summary(mc, diagnose = FALSE)
  for (dis in c("d1", "d3")) {
    g <- mc_group(mc, dis, "simultaneous")
    bm <- vapply(split(seq_len(nrow(g)), g$block),
                 function(ix) sum(g$w[ix] * g$adjusted[ix]) / sum(g$w[ix]), numeric(1))
    expect_length(bm, 5)
    expect_equal(bm, as.vector(tapply(g$adjusted, g$block, mean)), ignore_attr = TRUE)
    expect_equal(mc_row(sm, dis, "simultaneous")$mcse, stats::sd(bm) / sqrt(5))
    expect_equal(mc_row(sm, dis, "simultaneous")$mean, mean(g$adjusted))
  }

  # The number of blocks is capped at n_draws / 2.
  m10 <- cm_monte_carlo(s, 10, sampling = "lhs", lhs_replicates = 10, seed = 1)
  expect_equal(as.vector(table(m10$block)), rep(2, 5))
  expect_error(cm_monte_carlo(s, 3, sampling = "lhs"), "at least 4 draws")
  expect_error(cm_monte_carlo(function(i) supp_model(c(2.5, 5, 7.5)), 10, sampling = "lhs"),
               "cm_sampler")
})

test_that("LHS with importance sampling uses weighted block means", {
  s <- mc_smooth_sampler()
  mc <- cm_monte_carlo(s, 80, sampling = "lhs", lhs_replicates = 4,
                       proposal = list("impact:d1" = dist_normal(2.5, 1)), seed = 8)
  expect_equal(as.vector(table(mc$block)), rep(20, 4))
  # The proposal is stratified within each block.
  for (b in 1:4) {
    x <- mc$params[["impact:d1"]][mc$block == b]
    expect_equal(sort(floor(stats::pnorm(x, 2.5, 1) * 20)), 0:19)
  }
  g <- mc_group(mc, "d1", "simultaneous")
  expect_false(isTRUE(all.equal(g$w, rep(1 / 80, 80))))
  bm <- vapply(split(seq_len(nrow(g)), g$block),
               function(ix) sum(g$w[ix] * g$adjusted[ix]) / sum(g$w[ix]), numeric(1))
  Wb <- vapply(split(g$w, g$block), sum, numeric(1))
  mu <- sum(g$w * g$adjusted)
  expect_equal(sum(Wb * bm), mu)
  sm <- summary(mc, diagnose = FALSE)
  expect_equal(mc_row(sm, "d1", "simultaneous")$mcse, sqrt(4 / 3 * sum(Wb^2 * (bm - mu)^2)),
               ignore_attr = TRUE)
  expect_equal(mc_row(sm, "d1", "simultaneous")$mean, mu)
})

# ---- Stability statuses -------------------------------------------------------

test_that("stability statuses follow the documented rules", {
  check_rules <- function(sm, n) {
    expect_true(all(sm$stability %in% stability_levels))
    np <- sm$stability != "possible_pole"
    expected <- ifelse(n >= 50 & sm$tail_share > 0.6, "heavy_tail",
                ifelse(is.na(sm$mcse) | (is.na(sm$rel_mcse) & sm$mcse > 0), "insufficient_info",
                       ifelse(!is.na(sm$rel_mcse) & sm$rel_mcse > 0.05, "imprecise", "ok")))
    expect_equal(sm$stability[np], expected[np])
  }
  # ok: a smooth, precise estimate.
  s_ok <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d1 = dist_normal(2.5, 0.05)))
  mc_ok <- cm_monte_carlo(s_ok, 50, seed = 2)
  sm_ok <- summary(mc_ok, diagnose = FALSE)
  check_rules(sm_ok, 50)
  expect_true(all(sm_ok$stability == "ok"))
  dg <- cm_diagnose(mc_ok)
  expect_s3_class(dg, "cm_diagnosis")
  expect_named(dg, c("item", "method", "stability", "detail", "suggestion"))
  expect_equal(nrow(dg), 0)
  expect_output(print(dg), "No unstable Monte Carlo estimates")

  # imprecise: d1 alternates between 1 and 9 (deterministic).
  imp <- function(i) supp_model(c(if (i %% 2) 1 else 9, 5, 7.5))
  mc_imp <- cm_monte_carlo(imp, 50)
  sm_imp <- summary(mc_imp, diagnose = FALSE)
  check_rules(sm_imp, 50)
  r <- mc_row(sm_imp, "d1", "simultaneous")
  expect_equal(r$stability, "imprecise")
  dg <- cm_diagnose(mc_imp)
  d1 <- dg[dg$item == "d1", , drop = FALSE]
  expect_equal(d1$stability, "imprecise")
  need <- ceiling(50 * (r$rel_mcse / 0.02)^2)
  expect_match(d1$suggestion, format(need, big.mark = ",", scientific = FALSE), fixed = TRUE)
  expect_match(d1$suggestion, "lhs", fixed = TRUE)

  # heavy_tail: one extreme draw dominates the variance (deterministic).
  ht <- function(i) supp_model(c(if (i == 1) 1000 else 2.5 + i / 1000, 5, 7.5))
  mc_ht <- cm_monte_carlo(ht, 60)
  sm_ht <- summary(mc_ht, diagnose = FALSE)
  check_rules(sm_ht, 60)
  r <- mc_row(sm_ht, "d1", "simultaneous")
  expect_gt(r$tail_share, 0.6)
  expect_equal(r$stability, "heavy_tail")
  dg <- cm_diagnose(mc_ht)
  d1 <- dg[dg$item == "d1", , drop = FALSE]
  expect_equal(d1$stability, "heavy_tail")
  expect_match(d1$suggestion, 'cm_suggest_proposal(mc, "d1", method = "simultaneous")', fixed = TRUE)
  expect_message(summary(mc_ht), "may be unstable")
  expect_output(print(dg), "d1 \\(simultaneous\\)")
})

test_that("possible_pole is flagged when m + c changes sign against m", {
  # c1 = A12 * m2 (A13 = 0 because d1:d3 has odds ratio 1). With m2 = 5,
  # c1 is about 0.53, so m1 ~ N(-0.5, 0.5) puts m1 + c1 on both sides of 0,
  # and m1 + c1 > 0 > m1 in a share of the draws.
  s <- cm_sampler(supp_model(c(-0.5, 5, 7.5)), impacts = list(d1 = dist_normal(-0.5, 0.5)))
  mc <- cm_monte_carlo(s, 100, method = c("published", "simultaneous"), seed = 5)
  expect_equal(mc$n_rejected, 0)
  g <- mc_group(mc, "d1", "published")
  ok <- abs(g$raw) > 1e-12 & abs(g$adjusted) > 1e-12
  den <- g$raw[ok]^2 / g$adjusted[ok]
  expect_true(any(den > 0) && any(den < 0))
  flip <- mean(sign(den) != sign(g$raw[ok]))
  expect_gt(flip, 0)

  sm <- summary(mc, diagnose = FALSE)
  expect_equal(mc_row(sm, "d1", "published")$stability, "possible_pole")
  # The exact method has no division and is never flagged as a pole.
  expect_false("possible_pole" %in% sm$stability[sm$method == "simultaneous"])
  expect_false("possible_pole" %in% sm$stability[sm$disease != "d1"])
  dg <- cm_diagnose(mc)
  p <- dg[dg$stability == "possible_pole", , drop = FALSE]
  expect_equal(p$item, "d1")
  expect_equal(p$method, "published")
  expect_match(p$detail, sprintf("%.1f%%", 100 * flip), fixed = TRUE)
  expect_match(p$suggestion, "simultaneous", fixed = TRUE)
  expect_message(summary(mc), "possible|pole")
  # Contributions inherit the flag; totals are not checked for poles.
  sc <- summary(mc, what = "contribution", diagnose = FALSE)
  expect_equal(mc_row(sc, "d1", "published")$stability, "possible_pole")
  st <- summary(mc, what = "total", diagnose = FALSE)
  expect_false("possible_pole" %in% st$stability)
  expect_true(all(st$stability %in% stability_levels))
})

test_that("possible_pole is not flagged for removable or sign-stable denominators", {
  # Removable: the other impacts are 0, so c1 = 0 and m^2 / (m + c) = m.
  # m1 takes both signs and is exactly 0 in every fourth draw, where the
  # published method returns 0 instead of 0 / 0.
  rem <- function(i) supp_model(c(if (i %% 4 == 0) 0 else stats::rnorm(1, 0.2, 1), 0, 0))
  mc <- cm_monte_carlo(rem, 100, method = "published", seed = 6)
  expect_equal(mc$n_rejected, 0)
  g <- mc_group(mc, "d1", "published")
  expect_equal(sum(g$raw == 0), 25)
  expect_true(all(g$adjusted[g$raw == 0] == 0))
  expect_equal(g$adjusted, g$raw)
  expect_true(any(g$raw > 0) && any(g$raw < 0))
  sm <- summary(mc, diagnose = FALSE)
  expect_true(all(sm$stability %in% stability_levels))
  expect_false("possible_pole" %in% sm$stability)
  expect_equal(mc_row(sm, "d2", "published")$stability, "ok")

  # The denominator never changes sign (m1 + c1 > 0 in every draw).
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d1 = dist_normal(2.5, 0.5)))
  mc2 <- cm_monte_carlo(s, 100, method = "published", seed = 7)
  g2 <- mc_group(mc2, "d1", "published")
  expect_true(all(g2$raw^2 / g2$adjusted > 0))
  sm2 <- summary(mc2, diagnose = FALSE)
  expect_false("possible_pole" %in% sm2$stability)

  # A consistent sign change (m1 < 0 < m1 + c1 in every draw) is not a pole.
  s3 <- cm_sampler(supp_model(c(-0.3, 5, 7.5)), impacts = list(d1 = dist_normal(-0.3, 0.02)))
  mc3 <- cm_monte_carlo(s3, 60, method = "published", seed = 8)
  g3 <- mc_group(mc3, "d1", "published")
  expect_true(all(g3$raw < 0) && all(g3$raw^2 / g3$adjusted > 0))
  expect_false("possible_pole" %in% summary(mc3, diagnose = FALSE)$stability)
})

test_that("reported MCSEs are calibrated over repeated runs", {
  skip_on_cran()
  s <- cm_sampler(supp_model(c(2.5, 5, 7.5)),
                  associations = list("d1:d2" = dist_lognormal_ci(2, 1.5, 2.7)),
                  impacts = list(d1 = dist_normal(2.5, 0.5), d2 = dist_normal(5, 1)))
  for (smp in c("random", "lhs")) {
    est <- vapply(1:30, function(r) {
      mc <- cm_monte_carlo(s, 80, sampling = smp, lhs_replicates = 8, seed = 1000 + r)
      row <- mc_row(summary(mc, diagnose = FALSE), "d1", "simultaneous")
      c(mean = row$mean, mcse = row$mcse)
    }, c(mean = 0, mcse = 0))
    ratio <- stats::sd(est["mean", ]) / mean(est["mcse", ])
    expect_gt(ratio, 0.5, label = paste("sd(run means) / mean(mcse),", smp))
    expect_lt(ratio, 2, label = paste("sd(run means) / mean(mcse),", smp))
  }
})

# ---- Batch runs ---------------------------------------------------------------

test_that("batch runs share population draws across analyses", {
  an <- cm_analyses(supp_population(),
                    a = cm_impacts(ids3, c(2.5, 5, 7.5)),
                    b = cm_impacts(ids3, c(1, 2, 3)))
  bs <- cm_batch_sampler(an, diseases = list(d1 = dist_normal(0.1, 0.1)),
                         associations = list("d3:d2" = dist_lognormal_ci(3, 2, 4.5)),
                         impacts = list(a = list(d1 = dist_normal(2.5, 0.5))))
  expect_s3_class(bs, "cm_batch_sampler")
  expect_equal(bs$population_keys, c("prob:d1", "assoc:d2:d3"))
  expect_output(print(bs), "2 analyses \\(a, b\\); 2 shared population inputs")
  mc <- cm_monte_carlo(bs, 60, seed = 9)
  expect_s3_class(mc, "cm_mc_batch")
  ra <- mc$analyses$a
  rb <- mc$analyses$b
  # Negative prevalence draws are rejected in both analyses, for the same draws.
  expect_gt(ra$n_rejected, 0)
  expect_equal(ra$rejections$draw, rb$rejections$draw)
  keys <- c("draw", "prob:d1", "assoc:d2:d3")
  expect_equal(ra$params[, keys], rb$params[, keys], ignore_attr = TRUE)
  expect_gt(length(unique(ra$params[["impact:d1"]])), 1)
  expect_equal(unique(rb$params[["impact:d1"]]), 1)

  expect_error(cm_batch_sampler(supp_model(c(2.5, 5, 7.5))), "cm_analyses")
  expect_error(cm_batch_sampler(an, impacts = list(zz = list(d1 = dist_fixed(1)))),
               "Unknown analyses")
  expect_error(cm_monte_carlo(bs, 10, sampling = "lhs"), class = "deconflate_unsupported")
  expect_error(cm_monte_carlo(bs, 10, proposal = list("prob:d1" = dist_uniform(0, 1))),
               class = "deconflate_unsupported")
})

test_that("the global dairy batch sampler runs, summarises and prints", {
  bs <- sampler_global_dairy()
  expect_s3_class(bs, "cm_batch_sampler")
  expect_equal(names(bs$samplers), c("yield", "fertility"))
  expect_length(bs$population_keys, 35)
  expect_true(all(startsWith(bs$population_keys, "assoc:")))
  mc <- cm_monte_carlo(bs, 20, method = "published", seed = 1)
  expect_s3_class(mc, "cm_mc_batch")
  expect_equal(mc$n_draws, 20)
  expect_equal(mc$population_keys, bs$population_keys)
  py <- mc$analyses$yield$params
  pf <- mc$analyses$fertility$params
  expect_gt(nrow(py), 0)
  expect_setequal(py$draw, pf$draw)
  pop_cols <- grep("^(prob|assoc):", names(py), value = TRUE)
  expect_true(all(bs$population_keys %in% pop_cols))
  expect_equal(py[match(pf$draw, py$draw), pop_cols], pf[, pop_cols], ignore_attr = TRUE)
  # Disease probabilities are fixed in the analysis inputs; odds ratios vary.
  expect_length(unique(py[["prob:SCK"]]), 1)
  expect_gt(length(unique(py[["assoc:CK:SCK"]])), 1)
  # Each analysis draws its own impacts.
  expect_false(isTRUE(all.equal(py[["impact:CM"]], pf[["impact:CM"]])))

  s <- summary(mc, diagnose = FALSE)
  expect_equal(names(s)[1], "analysis")
  expect_equal(unique(s$analysis), c("yield", "fertility"))
  expect_equal(nrow(s), 2 * 12)
  expect_true(all(s$stability %in% stability_levels))
  st <- summary(mc, what = "total", diagnose = FALSE)
  expect_true(all(c("analysis", "quantity") %in% names(st)))
  expect_output(print(mc), "<cm_mc_batch> 20 draws, 2 analyses")
  expect_output(print(mc), "yield: \\d+ accepted")
})

# ---- Method comparison, diagnosis, proposals, reweighting, scenarios, gaps ----

test_that("compare_methods summarises a two-method run", {
  mc <- cm_monte_carlo(mc_smooth_sampler(), 100, method = c("published", "simultaneous"), seed = 3)
  cmp <- compare_methods(mc)
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$methods, c("published", "simultaneous"))
  expect_named(cmp$impacts, c("disease", "raw_mean", "published", "simultaneous"))
  expect_equal(cmp$impacts$disease, ids3)
  sm <- summary(mc, diagnose = FALSE)
  for (m in mc$method) {
    expect_equal(cmp$impacts[[m]], sm$mean[sm$method == m][match(ids3, sm$disease[sm$method == m])])
    expect_true(all(cmp$stability[[m]] %in% stability_levels))
  }
  g <- mc_group(mc, "d1", "published")
  expect_equal(cmp$impacts$raw_mean[1], mean(g$raw))
  expect_named(cmp$totals, c("quantity", "method", "mean", "q0.5", "trimmed_mean", "mcse", "stability"))
  med <- compare_methods(mc, stat = "median")
  expect_equal(med$impacts$simultaneous,
               sm$`q0.5`[sm$method == "simultaneous"][match(ids3, sm$disease[sm$method == "simultaneous"])])
  expect_output(print(cmp), "Monte Carlo: 100 draws \\(0 rejected\\)")
  expect_error(compare_methods(cm_monte_carlo(mc_smooth_sampler(), 10, seed = 1)), "one method")
})

test_that("cm_suggest_proposal returns a usable defensive mixture", {
  s <- cm_sampler(supp_model(c(0.2, 5, 7.5)),
                  impacts = list(d1 = dist_normal(0.2, 1.5), d3 = dist_fixed(7.5)))
  mc <- cm_monte_carlo(s, 100, seed = 4)
  expect_message(cm_suggest_proposal(mc, "d1"), "Proposal for impact:d1")
  prop <- suppressMessages(cm_suggest_proposal(mc, "d1"))
  # Point masses are never proposed.
  expect_named(prop, "impact:d1")
  expect_s3_class(prop[[1]], "cm_dist")
  expect_equal(prop[[1]]$type, "mixture")
  expect_equal(prop[[1]]$params$weights, c(0.5, 0.5))
  expect_equal(c(prop[[1]]$lower, prop[[1]]$upper), c(-Inf, Inf))
  expect_match(attr(prop, "explanation"), "50% its own distribution")

  mc_is <- cm_monte_carlo(s, 100, proposal = prop, seed = 5)
  expect_equal(mc_is$n_rejected, 0)
  # The defensive mixture bounds the weights by 1 / (1 - 0.5).
  expect_true(all(mc_is$log_weights <= log(2) + 1e-9))
  expect_gt(mc_is$ess, 0)
  expect_lte(mc_is$ess, 100)
  r <- mc_row(summary(mc_is, diagnose = FALSE), "d1", "simultaneous")
  expect_true(is.finite(r$mcse))

  expect_error(cm_suggest_proposal(cm_monte_carlo(s, 10, seed = 1), "d1"), "Too few")
  expect_error(cm_suggest_proposal(mc, "d1", weight = 1), "between 0 and 1")
  mc_fn <- cm_monte_carlo(function(i) supp_model(c(i, 5, 7.5)), 25)
  expect_error(cm_suggest_proposal(mc_fn, "d1"), "no recorded distributions")
})

test_that("cm_reweight keeps the weights with a zero log ratio", {
  s <- mc_smooth_sampler()
  mc <- cm_monte_carlo(s, 100, seed = 12)
  zero <- function(p) rep(0, nrow(p))
  mc2 <- cm_reweight(mc, zero)
  expect_equal(mc2$weights, mc$weights)
  expect_equal(mc2$ess, 100)
  expect_equal(summary(mc2, diagnose = FALSE)$mean, summary(mc, diagnose = FALSE)$mean)
  # Importance weights from the run are kept.
  mc_is <- cm_monte_carlo(s, 100, proposal = list("impact:d1" = dist_normal(2.5, 1)), seed = 13)
  mc_is2 <- cm_reweight(mc_is, zero)
  expect_equal(mc_is2$weights, mc_is$weights)
  expect_equal(mc_is2$ess, mc_is$ess)

  # Shifting impact:d1 upwards raises the weighted mean of d1.
  mc3 <- cm_reweight(mc, function(p) 2 * (p[["impact:d1"]] - 2.5))
  d1 <- function(x) mc_row(summary(x, diagnose = FALSE), "d1", "simultaneous")$mean
  expect_gt(d1(mc3), d1(mc))
  expect_lt(mc3$ess, 100)
  expect_warning(cm_reweight(mc, function(p) 50 * p[["impact:d1"]]), class = "deconflate_low_ess")

  expect_error(cm_reweight(mc, function(p) 0), "one value per accepted draw")
  expect_error(cm_reweight(mc, function(p) rep(Inf, nrow(p))), class = "deconflate_nonfinite")
  expect_error(cm_reweight(mc, function(p) rep(-Inf, nrow(p))), class = "deconflate_nonfinite")
  expect_error(cm_reweight(list(), zero), "cm_monte_carlo")
})

test_that("cm_scenario reweights draws to scenario distributions", {
  s <- mc_smooth_sampler()
  mc <- cm_monte_carlo(s, 150, seed = 14)
  sc <- cm_scenario(mc, list("impact:d1" = dist_normal(2.7, 0.5)))
  x <- mc$params[["impact:d1"]]
  lw <- stats::dnorm(x, 2.7, 0.5, log = TRUE) - stats::dnorm(x, 2.5, 0.5, log = TRUE)
  w <- exp(lw - max(lw))
  expect_equal(sc$weights, w / sum(w))
  expect_lt(sc$ess, 150)
  g <- mc_group(sc, "d1", "simultaneous")
  expect_gt(sum(g$w * g$raw), mean(g$raw))
  # The unchanged input distribution gives uniform weights.
  same <- cm_scenario(mc, list("impact:d1" = dist_normal(2.5, 0.5)))
  expect_equal(same$weights, mc$weights)

  # On an importance-sampled run, the draws were sampled from the proposal.
  mc_is <- cm_monte_carlo(s, 100, proposal = list("impact:d1" = dist_normal(2.5, 1)), seed = 15)
  back <- cm_scenario(mc_is, list("impact:d1" = dist_normal(2.5, 0.5)))
  expect_equal(back$weights, mc_is$weights)

  expect_error(cm_scenario(mc, list("impact:zz" = dist_normal(1, 1))), "No sampled input")
  expect_error(cm_scenario(mc, list("impact:d1" = dist_fixed(2.5))), class = "deconflate_unsupported")
  su <- cm_sampler(supp_model(c(2.5, 5, 7.5)), impacts = list(d1 = dist_uniform(2, 3)))
  mu <- cm_monte_carlo(su, 30, seed = 16)
  expect_error(cm_scenario(mu, list("impact:d1" = dist_normal(2.5, 0.1))), "is not covered",
               class = "deconflate_unsupported")
  expect_s3_class(cm_scenario(mu, list("impact:d1" = dist_uniform(2.2, 2.8))), "cm_mc")
  mc_fn <- cm_monte_carlo(function(i) supp_model(c(2.5, 5, 7.5)), 5)
  expect_error(cm_scenario(mc_fn, list("impact:d1" = dist_normal(2.5, 0.5))),
               "no recorded distributions")
})

test_that("each draw's contributions and totals match deconflate() on its own model", {
  s <- mc_smooth_sampler()
  mc <- cm_monte_carlo(s, 50, seed = 17)
  # Contributions add up to the adjusted total of each draw.
  sums <- tapply(mc$draws$contribution, mc$draws$draw, sum)
  expect_equal(as.vector(sums[as.character(mc$totals$draw)]), mc$totals$adjusted_total,
               ignore_attr = TRUE)
  # Same as deconflate() on the draw's own model.
  keys <- names(attr(s, "specs"))
  v1 <- unlist(mc$params[1, keys])
  r1 <- deconflate(s(1, values = v1))
  expect_equal(mc$totals$adjusted_total[1], r1$totals$adjusted_total, ignore_attr = TRUE)
  expect_equal(mc$draws$contribution[mc$draws$draw == mc$params$draw[1]], r1$contributions$total,
               ignore_attr = TRUE)
  # The Monte Carlo gap helper was removed with valuation.
  expect_false(exists("cm_mc_gap", envir = asNamespace("deconflate"), inherits = FALSE))
})

# ---- Reproducibility ----------------------------------------------------------

test_that("seeds reproduce runs exactly", {
  s <- mc_smooth_sampler()
  same <- function(a, b) {
    expect_identical(a$draws, b$draws)
    expect_identical(a$params, b$params)
    expect_identical(a$weights, b$weights)
    expect_identical(a$rejections, b$rejections)
  }
  same(cm_monte_carlo(s, 30, seed = 42), cm_monte_carlo(s, 30, seed = 42))
  same(cm_monte_carlo(s, 30, sampling = "lhs", lhs_replicates = 3, seed = 42),
       cm_monte_carlo(s, 30, sampling = "lhs", lhs_replicates = 3, seed = 42))
  prop <- list("impact:d1" = dist_normal(2.5, 1))
  same(cm_monte_carlo(s, 30, proposal = prop, seed = 42),
       cm_monte_carlo(s, 30, proposal = prop, seed = 42))
  expect_false(identical(cm_monte_carlo(s, 30, seed = 42)$params,
                         cm_monte_carlo(s, 30, seed = 43)$params))

  an <- cm_analyses(supp_population(), a = cm_impacts(ids3, c(2.5, 5, 7.5)),
                    b = cm_impacts(ids3, c(1, 2, 3)))
  bs <- cm_batch_sampler(an, associations = list("d1:d2" = dist_lognormal_ci(2, 1.5, 2.7)),
                         impacts = list(b = list(d2 = dist_normal(2, 0.5))))
  b1 <- cm_monte_carlo(bs, 15, seed = 42)
  b2 <- cm_monte_carlo(bs, 15, seed = 42)
  same(b1$analyses$a, b2$analyses$a)
  same(b1$analyses$b, b2$analyses$b)
})
