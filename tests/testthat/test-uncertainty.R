# Uncertainty: draws of the inputs with a distribution in deconflate(n_draws =)
# and their summary statistics (draw_stats(), R/draws.R).

# The supplement model with two uncertain impacts.
unc_model <- function(d1 = dist_normal(2.5, 0.5), d3 = dist_normal(7.5, 1)) {
  cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5), units = "%"),
           distributions = list("impact:d1" = d1, "impact:d3" = d3))
}

# The supplement population with event impacts (hazard ratios), d1 uncertain.
unc_event_model <- function() {
  cm_model(supp_population(),
           cm_impacts(ids3, c(1.5, 2.0, 1.3), measure = "HR", label = "culling"),
           distributions = list("impact:d1" = dist_lognormal_ci(1.5, 1.2, 1.9)))
}

# Summary row of one quantity.
qrow <- function(s, q) s[s$quantity == q, , drop = FALSE]

stability_levels <- c("ok", "imprecise", "insufficient_info", "heavy_tail")

# ---- Distributions in the model ------------------------------------------------

test_that("distributions are keyed by input and checked against the model", {
  m <- cm_model(supp_population(), cm_impacts(c(ids3, "d1;d2"), c(2.5, 5, 7.5, 0.5),
                                              adjusted_for = c(NA, NA, NA, "all")),
                distributions = list("prob:d1" = dist_fixed(0.12), "assoc:d2:d1" = dist_fixed(2.5),
                                     "assoc:d3:d2" = dist_fixed(1.5), "impact:d3" = dist_fixed(8),
                                     "inter:d2:d1" = dist_fixed(1)))
  # Keys in canonical form: the table and row.
  expect_equal(names(m$distributions),
               c("prob:d1", "associations:1", "associations:3", "impacts:3", "impacts:4"))
  expect_true(all(vapply(m$distributions, inherits, logical(1), "cm_dist")))
  expect_output(print(m), "Uncertain inputs \\(with a distribution\\): 5")
  # set_inputs() puts drawn values into a copy of the model, by canonical
  # key or alias.
  d <- set_inputs(m, c("prob:d1" = 0.12, "associations:1" = 2.5, "assoc:d2:d3" = 1.5,
                       "impacts:3" = 8, "inter:d1:d2" = 1))
  expect_s3_class(d, "cm_model")
  expect_equal(d$diseases$value, c(0.12, 0.15, 0.20))
  expect_equal(d$diseases$prob, c(0.12, 0.15, 0.20))
  expect_equal(d$associations$value, c(2.5, 1, 1.5))
  expect_equal(d$impacts$value, c(2.5, 5, 8, 1))
  expect_equal(m$impacts$value, c(2.5, 5, 7.5, 0.5))
  # Incidence rates are drawn on their own scale and converted.
  m_ir <- cm_model(cm_diseases(ids3, c(0.10, 0.15, 0.20), measure = "incidence_rate"),
                   cm_impacts(ids3, c(2.5, 5, 7.5)),
                   associations = cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)),
                   distributions = list("prob:d2" = dist_uniform(0.4, 0.6)))
  d_ir <- set_inputs(m_ir, c("prob:d2" = 0.5))
  expect_equal(d_ir$diseases$value[2], 0.5)
  expect_equal(d_ir$diseases$prob[2], 1 - exp(-0.5))

  mk <- function(dists, model = supp_model(c(2.5, 5, 7.5))) {
    cm_model(as_population(model), model$impacts, distributions = dists)
  }
  expect_error(mk(list("impact:d1" = 3)), "named list of distributions")
  expect_error(mk(list(dist_fixed(3))), "named list of distributions")
  expect_error(mk(list("prob:zz" = dist_fixed(0.1))), "unknown disease")
  expect_error(mk(list("impact:zz" = dist_fixed(1))), "no impact row matches")
  expect_error(mk(list("assoc:d1-d2" = dist_fixed(2))), "assoc:<d1>:<d2>")
  expect_error(mk(list("assoc:d1:d2:d3" = dist_fixed(2))), "assoc:<d1>:<d2>")
  expect_error(mk(list("assoc:d1:d2" = dist_fixed(2), "assoc:d2:d1" = dist_fixed(3))), "more than once")
  expect_error(mk(list("impact:d1" = dist_fixed(1), "impact:d1" = dist_fixed(2))), "more than once")
  expect_error(mk(list("inter:d1:d2" = dist_fixed(1))), "no impact row matches")
  expect_error(mk(list("assoc:d1:d3" = dist_fixed(1)),
                  cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5)))),
               "no association row")
  expect_error(mk(list("assoc:d1:d2" = dist_fixed(1)),
                  cm_model(cm_population(cm_diseases(ids3, c(0.1, 0.15, 0.2))), cm_impacts(ids3, c(2.5, 5, 7.5)))),
               "no association")
  expect_error(mk(list("foo:d1" = dist_fixed(1))), "unknown input type")
  expect_error(cm_model(supp_population(), distributions = list("impact:d1" = dist_fixed(1))),
               "no impacts")
})

test_that("impossible drawn values raise deconflate_infeasible", {
  m <- supp_model(c(2.5, 5, 7.5))
  expect_error(set_inputs(m, c("prob:d1" = -0.01)), "not a valid disease value", class = "deconflate_infeasible")
  expect_error(set_inputs(m, c("prob:d1" = 1.2)), "not a valid disease value", class = "deconflate_infeasible")
  expect_error(set_inputs(m, c("assoc:d1:d2" = -0.5)), "is not a valid OR", class = "deconflate_infeasible")
  expect_error(set_inputs(m, c("impact:d1" = Inf)), "not finite", class = "deconflate_infeasible")
  e <- supp_hr_model(c(1.5, 2, 0.05), measure = c("HR", "HR", "RD"))
  expect_error(set_inputs(e, c("impact:d1" = -1)), "not a valid HR", class = "deconflate_infeasible")
  expect_error(set_inputs(e, c("impacts:3" = 1.5)), "not a valid risk difference",
               class = "deconflate_infeasible")
  expect_equal(set_inputs(e, c("impact:d3" = -0.05))$impacts$value[3], -0.05)
  expect_error(set_inputs(m, c("foo:d1" = 1)), "unknown input type")
  p <- cm_model(supp_population(cm_associations("d1", "d2", 0.3, measure = "prevalence", among = "d3")),
                m$impacts)
  expect_error(set_inputs(p, c("associations:4" = 1)), "not a valid prevalence", class = "deconflate_infeasible")
})

test_that("three-way ratios (OR_ratio rows) can have distributions, in any order of the diseases", {
  pop <- supp_population(cm_associations("d2", "d1", 2, measure = "OR_ratio", among = "d3"))
  m <- cm_model(pop, cm_impacts(c(ids3, "d1;d2"), c(2.5, 5, 7.5, 0.5), adjusted_for = c(NA, NA, NA, "all")),
                distributions = list("three:d3:d1:d2" = dist_lognormal_ci(2, 1.2, 3.5)))
  expect_equal(names(m$distributions), "associations:4")
  expect_equal(set_inputs(m, c("associations:4" = 1.7))$associations$value[4], 1.7)
  expect_error(set_inputs(m, c("associations:4" = 0)), "not a valid OR_ratio", class = "deconflate_infeasible")
  mk <- function(d) cm_model(pop, m$impacts, distributions = d)
  expect_error(mk(list("three:d1:d2" = dist_fixed(2))), "three:<d1>:<d2>:<d3>")
  expect_error(mk(list("three:d1:d2:d4" = dist_fixed(2))), "no association row matches")
  expect_error(mk(list("three:d1:d2:d3" = dist_fixed(2), "three:d3:d2:d1" = dist_fixed(3))),
               "more than once")
  res <- deconflate(m, n_draws = 20, seed = 1)
  expect_equal(res$method, "full model")
  expect_equal(res$draws$n_rejected, 0)
  expect_gt(length(unique(res$draws$params[, "associations:4"])), 1)
  # The ratio changes the result when interactions are present.
  lo <- deconflate(set_inputs(m, c("associations:4" = 0.5)), n_draws = 0)$totals$adjusted_total
  hi <- deconflate(set_inputs(m, c("associations:4" = 4)), n_draws = 0)$totals$adjusted_total
  expect_false(isTRUE(all.equal(lo, hi)))
  expect_gt(stats::sd(res$draws$values[, "total"]), 0)
})

# ---- Draws in deconflate() --------------------------------------------------------

test_that("deconflate() draws the uncertain inputs and summarises the results", {
  m <- unc_model()
  res <- deconflate(m, n_draws = 100, seed = 3)
  expect_equal(res$method, "pairwise shortcut")
  dr <- res$draws
  expect_true(all(c("summary", "values", "params", "draw", "n_draws", "n_rejected", "rejections",
                    "seed", "sampling", "specs") %in% names(dr)))
  expect_equal(dr$n_draws, 100)
  expect_equal(dr$n_rejected, 0)
  expect_equal(dr$seed, 3)
  expect_equal(dr$sampling, "random")
  expect_null(dr$block)
  expect_equal(dr$draw, 1:100)
  expect_named(dr$rejections, c("draw", "type", "reason"))
  expect_equal(nrow(dr$rejections), 0)
  expect_equal(names(dr$specs), c("impacts:1", "impacts:3"))
  qn <- c(paste0("adjusted:", ids3), paste0("contribution:", ids3), "total", "raw_total")
  expect_equal(colnames(dr$values), qn)
  expect_equal(dim(dr$values), c(100L, 8L))
  expect_equal(colnames(dr$params), c("impacts:1", "impacts:3"))

  # Draws are quantiles of the distributions at seeded uniforms.
  set.seed(3)
  u <- matrix(stats::runif(200), 100, 2)
  expect_equal(unname(dr$params[, "impacts:1"]), stats::qnorm(u[, 1], 2.5, 0.5))
  expect_equal(unname(dr$params[, "impacts:3"]), stats::qnorm(u[, 2], 7.5, 1))

  # Each draw is the exact adjustment of its own inputs (linear here).
  A <- res$conflation$A
  raw <- cbind(dr$params[, "impacts:1"], 5, dr$params[, "impacts:3"])
  expect_equal(unname(dr$values[, paste0("adjusted:", ids3)]), unname(t(solve(A, t(raw)))), tolerance = 1e-10)
  expect_equal(unname(rowSums(dr$values[, paste0("contribution:", ids3)])), unname(dr$values[, "total"]),
               tolerance = 1e-12)
  expect_equal(unname(dr$values[, "raw_total"]), as.vector(raw %*% c(0.10, 0.15, 0.20)), tolerance = 1e-12)

  # Summary: mean, SD (divisor n), MCSE = SD / sqrt(n), quantiles, stability.
  s <- dr$summary
  expect_named(s, c("quantity", "mean", "sd", "mcse", "lower", "median", "upper", "stability"))
  expect_equal(s$quantity, qn)
  expect_true(all(s$stability %in% stability_levels))
  for (q in qn) {
    x <- dr$values[, q]
    mu <- mean(x)
    r <- qrow(s, q)
    expect_equal(r$mean, mu, info = q)
    expect_equal(r$sd, sqrt(mean((x - mu)^2)), info = q)
    expect_equal(r$mcse, sqrt(sum((x - mu)^2)) / 100, info = q)
    expect_equal(c(r$lower, r$median, r$upper), sort(x)[c(3, 50, 98)], info = q)
  }

  # The central estimate is the point estimate at the input values, not the
  # mean of the draws.
  pt <- deconflate(m, n_draws = 0)
  expect_equal(res$adjusted$adjusted, pt$adjusted$adjusted)
  expect_equal(res$totals$adjusted_total, pt$totals$adjusted_total)
  expect_false(isTRUE(all.equal(res$adjusted$adjusted, qrow(s, "adjusted:d1")$mean)))
  # Interval columns.
  expect_equal(res$adjusted$lower, s$lower[match(paste0("adjusted:", ids3), s$quantity)])
  expect_equal(res$adjusted$upper, s$upper[match(paste0("adjusted:", ids3), s$quantity)])
  expect_equal(res$contributions$lower, s$lower[match(paste0("contribution:", ids3), s$quantity)])
  expect_equal(res$contributions$upper, s$upper[match(paste0("contribution:", ids3), s$quantity)])
  expect_equal(res$totals$adjusted_total_lower, qrow(s, "total")$lower)
  expect_equal(res$totals$adjusted_total_upper, qrow(s, "total")$upper)
  expect_true(all(res$adjusted$lower <= res$adjusted$upper))
  # d2 has no distribution, but its adjusted impact depends on d1 and d3.
  expect_gt(qrow(s, "adjusted:d2")$sd, 0)
  expect_output(print(res), "95% intervals from 100 draws \\(0 rejected; random sampling; seed 3\\)")
  expect_output(print(res), "95% interval")

  # Point estimates only.
  expect_null(pt$draws)
  expect_null(pt$adjusted$lower)
  expect_null(pt$totals$adjusted_total_lower)
  expect_match(pt$notes, "a distribution, but n_draws = 0: the results are point estimates without intervals", fixed = TRUE)
  expect_error(deconflate(m, n_draws = -1), "whole number")
  expect_error(deconflate(m, n_draws = c(10, 20)), "whole number")
  expect_error(deconflate(m, n_draws = 10, seed = "a"), "`seed` must be a single whole number")
})

test_that("the same seed gives the same draws, and an unseeded run stores its seed", {
  m <- unc_model()
  a <- deconflate(m, n_draws = 30, seed = 42)$draws
  b <- deconflate(m, n_draws = 30, seed = 42)$draws
  expect_identical(a$values, b$values)
  expect_identical(a$params, b$params)
  expect_identical(a$summary, b$summary)
  expect_false(identical(a$params, deconflate(m, n_draws = 30, seed = 43)$draws$params))
  l1 <- deconflate(m, n_draws = 30, seed = 42, sampling = "lhs", lhs_replicates = 3)$draws
  l2 <- deconflate(m, n_draws = 30, seed = 42, sampling = "lhs", lhs_replicates = 3)$draws
  expect_identical(l1$values, l2$values)
  expect_identical(l1$block, l2$block)
  r0 <- deconflate(m, n_draws = 10)
  s0 <- r0$draws$seed
  expect_true(is.numeric(s0) && length(s0) == 1L && is.finite(s0))
  expect_identical(deconflate(m, n_draws = 10, seed = s0)$draws$values, r0$draws$values)
})

test_that("without distributions there are no draws, and a note says so", {
  res <- deconflate(example_supplement())
  expect_null(res$draws)
  expect_null(res$adjusted$lower)
  expect_match(res$notes, "No input has a distribution, so no draws were run", all = FALSE)
  expect_match(summary(res)$notes, "No input has a distribution", all = FALSE)
  expect_output(print(res), "No input has a distribution")
  expect_length(deconflate(example_supplement(), n_draws = 0)$notes, 0)
})

test_that("Latin hypercube sampling uses replicate blocks and block-mean standard errors", {
  m <- unc_model(d1 = dist_uniform(2, 3), d3 = dist_normal(7.5, 1))
  res <- deconflate(m, n_draws = 100, sampling = "lhs", lhs_replicates = 5, seed = 7)
  dr <- res$draws
  expect_equal(dr$sampling, "lhs")
  expect_equal(dr$n_rejected, 0)
  expect_equal(dr$n_blocks, 5)
  expect_length(dr$block, 100)
  expect_equal(as.vector(table(dr$block)), rep(20, 5))
  expect_output(print(res), "Latin hypercube sampling")
  # Each block is stratified: one draw in each of 20 equal-probability strata.
  for (b in 1:5) {
    sel <- dr$block == b
    expect_equal(sort(floor((dr$params[sel, "impacts:1"] - 2) * 20)), 0:19)
    expect_equal(sort(floor(stats::pnorm(dr$params[sel, "impacts:3"], 7.5, 1) * 20)), 0:19)
  }
  for (q in c("adjusted:d1", "adjusted:d3", "total")) {
    x <- dr$values[, q]
    bm <- as.vector(tapply(x, dr$block, mean))
    expect_equal(qrow(dr$summary, q)$mean, mean(x))
    expect_equal(qrow(dr$summary, q)$mcse, stats::sd(bm) / sqrt(5), info = q)
  }
  # The number of blocks is capped at n_draws / 2.
  r10 <- deconflate(m, n_draws = 10, sampling = "lhs", lhs_replicates = 10, seed = 1)
  expect_equal(as.vector(table(r10$draws$block)), rep(2, 5))
  expect_error(deconflate(m, n_draws = 3, sampling = "lhs"), "at least 4 draws")
  expect_error(deconflate(m, n_draws = 10, sampling = "sobol"))
})

test_that("rejected draws are counted, tabulated and not replaced", {
  m <- cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5)),
                distributions = list("prob:d1" = dist_normal(0.1, 0.06),
                                     "assoc:d1:d2" = dist_normal(2, 1.2),
                                     "impact:d3" = dist_normal(7.5, 1)))
  res <- deconflate(m, n_draws = 300, seed = 5)
  dr <- res$draws
  expect_gt(dr$n_rejected, 0)
  expect_equal(dr$n_rejected + nrow(dr$params), 300)
  expect_equal(nrow(dr$values), nrow(dr$params))
  expect_equal(nrow(dr$rejections), dr$n_rejected)
  expect_named(dr$rejections, c("draw", "type", "reason"))
  expect_true(all(dr$rejections$type == "infeasible"))
  expect_true(any(grepl("not a valid disease value", dr$rejections$reason)))
  expect_true(any(grepl("is not a valid OR", dr$rejections$reason)))
  # Rejected draws are not replaced and leave no trace in the accepted ones.
  expect_length(intersect(dr$rejections$draw, dr$draw), 0)
  expect_setequal(c(dr$rejections$draw, dr$draw), 1:300)
  expect_true(all(dr$params[, "prob:d1"] > 0))
  expect_true(all(dr$params[, "associations:1"] > 0))
  expect_output(print(res), sprintf("\\(%d rejected;", dr$n_rejected))

  # A high rejection share gets a note (print and summary).
  mh <- cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5)),
                 distributions = list("prob:d1" = dist_normal(0.05, 0.1)))
  rh <- deconflate(mh, n_draws = 100, seed = 1)
  expect_gt(rh$draws$n_rejected, 10)
  expect_true(any(grepl("of the draws were rejected \\(infeasible", summary(rh)$notes)))
  expect_output(print(rh), "of the draws were rejected")
  # No note for a low share.
  expect_false(any(grepl("were rejected", summary(res)$notes)) && dr$n_rejected <= 30)

  # Every draw rejected: no intervals.
  ma <- cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5)),
                 distributions = list("prob:d1" = dist_uniform(1.1, 1.5)))
  ra <- deconflate(ma, n_draws = 5, seed = 1)
  expect_equal(ra$draws$n_rejected, 5)
  expect_null(ra$draws$values)
  expect_null(ra$draws$summary)
  expect_equal(nrow(ra$draws$params), 0)
  expect_null(ra$adjusted$lower)
  expect_equal(ra$adjusted$adjusted, deconflate(example_supplement(), n_draws = 0)$adjusted$adjusted)
  expect_match(summary(ra)$notes, "All 5 draws were rejected", all = FALSE)
  expect_output(print(ra), "All 5 draws were rejected")
})

test_that("rejections are typed by the class of the failure", {
  m <- unc_model()
  ev <- function(mm, r) {
    v <- mm$impacts$value[1]
    if (v > 2.5) cm_abort("too large", class = "deconflate_unsupported")
    if (v < 2) return(c(a = Inf))
    c(a = v)
  }
  dr <- run_draws(m, m$distributions, 60, 1, "random", 10L, ev)
  p <- deconflate(m, n_draws = 60, seed = 1)$draws$params[, "impacts:1"]
  expect_equal(dr$n_rejected, sum(p > 2.5 | p < 2))
  expect_setequal(unique(dr$rejections$type), c("unsupported", "nonfinite"))
  expect_equal(dr$rejections$draw[dr$rejections$type == "unsupported"], which(p > 2.5))
  expect_equal(dr$rejections$draw[dr$rejections$type == "nonfinite"], which(p < 2))
  expect_match(dr$rejections$reason[dr$rejections$type == "nonfinite"], "Non-finite")
  expect_true(all(dr$values[, "a"] >= 2 & dr$values[, "a"] <= 2.5))
  # Errors that are not deconflate conditions are not swallowed.
  expect_error(run_draws(m, m$distributions, 5, 1, "random", 10L, function(mm, r) stop("boom")), "boom")
  expect_equal(condition_type(simpleError("x")), "error")
})

test_that("a precision note is given when the Monte Carlo error is large", {
  m <- unc_model(d1 = dist_normal(0.2, 3))
  res <- deconflate(m, n_draws = 20, seed = 1)
  expect_true("imprecise" %in% res$draws$summary$stability)
  expect_true(any(grepl("Too few draws for precise intervals of the adjusted impacts? of d1", summary(res)$notes)))
  expect_match(printed(res), "Too few draws for precise intervals", fixed = TRUE)
})

test_that("stability statuses follow the documented rules", {
  stab <- function(x, block = NULL) {
    summarise_draws(list(values = cbind(q = x), block = block,
                         n_blocks = if (!is.null(block)) length(unique(block))))
  }
  expect_equal(stab(2.5 + seq(-0.01, 0.01, length.out = 50))$stability, "ok")
  imp <- stab(rep(c(1, 9), 25))
  expect_equal(imp$stability, "imprecise")
  expect_equal(imp$mcse, 4 / sqrt(50))
  expect_equal(stab(c(1000, 2.5 + (2:60) / 1000))$stability, "heavy_tail")
  expect_equal(stab(rep(c(-1, 1), 5))$stability, "insufficient_info")    # mean 0
  expect_equal(stab(1:10, block = rep(1, 10))$stability, "insufficient_info")  # one block
  expect_true(is.na(stab(1:10, block = rep(1, 10))$mcse))
  x <- c(3.1, 2.7, 4.4, 3.9, 2.2, 3.3, 4.1, 2.9, 3.6, 3.0, 2.5, 4.8, 3.4, 3.8, 2.6, 4.0, 3.2, 3.7, 2.8, 3.5)
  b <- rep(1:5, each = 4)
  expect_equal(stab(x, b)$mcse, stats::sd(as.vector(tapply(x, b, mean))) / sqrt(5))
  s <- stab(rev(seq_len(100)) / 10)
  expect_equal(c(s$lower, s$median, s$upper), c(0.3, 5.0, 9.8))
  expect_null(summarise_draws(list(values = NULL)))
  # draw_stats() directly: equal weights, MCSE = SD / sqrt(n) without blocks.
  ds <- draw_stats(x, NA, c(0.025, 0.5, 0.975))
  expect_named(ds, c("mean", "sd", "mcse", "q", "stability"))
  expect_equal(ds$mean, mean(x))
  expect_equal(ds$sd, sqrt(mean((x - mean(x))^2)))
  expect_equal(ds$mcse, ds$sd / sqrt(20))
  expect_length(ds$q, 3)
  expect_true(ds$stability %in% stability_levels)
  # Blocks with no accepted draw still count as replicates (n_blocks).
  expect_equal(draw_stats(x, b, c(0.025, 0.5, 0.975), n_blocks = 5)$mcse, stab(x, b)$mcse)
  expect_false(isTRUE(all.equal(draw_stats(x, b, c(0.025, 0.5, 0.975), n_blocks = 6)$mcse,
                                stab(x, b)$mcse)))
  # One block, even with others counted, gives no precision information.
  one <- draw_stats(1:4, rep(1L, 4), c(0.025, 0.5, 0.975), n_blocks = 4L)
  expect_true(is.na(one$mcse))
  expect_equal(one$stability, "insufficient_info")
})

test_that("each draw matches deconflate() on its own inputs (full model, unknown pair)", {
  m <- cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5)),
                distributions = list("assoc:d1:d2" = dist_lognormal_ci(2, 1.5, 2.7),
                                     "impact:d2" = dist_normal(5, 1)))
  res <- deconflate(m, n_draws = 20, seed = 1)
  expect_equal(res$method, "full model")
  expect_equal(res$draws$n_rejected, 0)
  expect_equal(nrow(res$unknown_pairs), 1)
  expect_gt(stats::sd(res$draws$values[, "adjusted:d1"]), 0)
  for (d in c(1, 7)) {
    rd <- deconflate(set_inputs(m, res$draws$params[d, ]), n_draws = 0)
    expect_equal(res$draws$values[d, ], result_quantities(rd), tolerance = 1e-7)
  }
})

# ---- Event impacts ------------------------------------------------------------------

test_that("event-model draws, with the overall risk as a distribution", {
  m <- unc_event_model()
  res <- deconflate(m, event_model = TRUE, overall_risk = dist_beta(250, 750), n_draws = 50, seed = 2)
  expect_s3_class(res, "cm_event_result")
  expect_equal(res$method, "snapshot")
  # The central estimate uses the mean of the overall risk (0.25)
  # (reference_v040.py, section 5).
  expect_equal(res$overall_risk, 0.25)
  expect_equal(res$adjusted$adjusted, c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 1e-7)
  expect_equal(res$attributable$summary$attributable, 0.0366489851, tolerance = 1e-8)
  dr <- res$draws
  expect_equal(colnames(dr$params), c("impacts:1", "risk"))
  expect_true(all(dr$params[, "risk"] > 0 & dr$params[, "risk"] < 1))
  expect_gt(length(unique(dr$params[, "risk"])), 1)
  qn <- c(paste0("adjusted:", ids3), "total", "fraction", paste0("contribution:", ids3), "raw_total")
  expect_equal(colnames(dr$values), qn)
  expect_equal(dr$summary$quantity, qn)
  expect_equal(unname(dr$values[, "fraction"]), unname(dr$values[, "total"] / dr$params[, "risk"]),
               tolerance = 1e-12)
  expect_equal(unname(rowSums(dr$values[, paste0("contribution:", ids3)])), unname(dr$values[, "total"]),
               tolerance = 1e-10)
  # Interval columns.
  s <- dr$summary
  expect_equal(res$adjusted$lower, s$lower[match(paste0("adjusted:", ids3), s$quantity)])
  expect_equal(res$adjusted$upper, s$upper[match(paste0("adjusted:", ids3), s$quantity)])
  expect_equal(res$attributable$summary$attributable_lower, qrow(s, "total")$lower)
  expect_equal(res$attributable$summary$attributable_upper, qrow(s, "total")$upper)
  expect_equal(res$contributions$lower, s$lower[match(paste0("contribution:", ids3), s$quantity)])
  expect_equal(res$contributions$upper, s$upper[match(paste0("contribution:", ids3), s$quantity)])
  expect_output(print(res), "95% interval")

  # Hazard-ratio-only models: the adjusted hazard ratios do not depend on the
  # overall risk, only the attributable risk does.
  rr <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = dist_beta(250, 750),
                   n_draws = 20, seed = 1)
  expect_equal(names(rr$draws$specs), "risk")
  for (i in 1:3) {
    expect_equal(unname(rr$draws$values[, paste0("adjusted:", ids3[i])]), rep(rr$adjusted$adjusted[i], 20),
                 tolerance = 1e-8)
  }
  expect_gt(stats::sd(rr$draws$values[, "total"]), 0)
  r30 <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.30, n_draws = 0)
  expect_equal(r30$adjusted$adjusted, rr$adjusted$adjusted, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(r30$attributable$summary$attributable,
                                rr$attributable$summary$attributable)))
  # A fixed overall risk and no distributions: no draws.
  r0 <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.25)
  expect_null(r0$draws)
  expect_match(r0$notes, "No input has a distribution", all = FALSE)

  # The overall risk must be a proportion, or a distribution within (0, 1).
  expect_error(deconflate(m, event_model = TRUE, overall_risk = dist_normal(0.25, 0.05)),
               "must lie between 0 and 1")
  expect_error(deconflate(m, event_model = TRUE, overall_risk = 1.5), "strictly between 0 and 1")
  expect_error(deconflate(m, event_model = TRUE, overall_risk = c(0.2, 0.3)), "strictly between 0 and 1")
  expect_error(deconflate(m, event_model = TRUE), "needs `overall_risk`")
  expect_s3_class(deconflate(m, event_model = TRUE, overall_risk = dist_normal(0.25, 0.05, lower = 0, upper = 1),
                             n_draws = 0),
                  "cm_event_result")
})

test_that("event-model draws of the five-disease culling table", {
  skip_on_cran()
  cu <- read_five("culling.csv")
  res <- deconflate(cu, event_model = TRUE, overall_risk = dist_beta(250, 750), n_draws = 20, seed = 1)
  dr <- res$draws
  expect_equal(colnames(dr$params), c(names(cu$distributions), "risk"))
  expect_equal(dr$n_rejected + nrow(dr$params), 20)
  expect_gt(nrow(dr$params), 10)
  # Central values (reference_v041.py), at the mean overall risk 0.25.
  expect_equal(res$adjusted$adjusted, c(1.060298968, 1.2342866774, 1.080381173, 1.1753668593, 1.509273874,
                                        2.2779452712), tolerance = 1e-6)
  expect_true(all(res$adjusted$lower <= res$adjusted$upper))
  expect_false(is.null(res$attributable$summary$attributable_lower))
})

# ---- Monte Carlo error ----------------------------------------------------------------

test_that("reported MCSEs are calibrated over repeated runs", {
  skip_on_cran()
  m <- cm_model(supp_population(), cm_impacts(ids3, c(2.5, 5, 7.5)),
                distributions = list("assoc:d1:d2" = dist_lognormal_ci(2, 1.5, 2.7),
                                     "impact:d1" = dist_normal(2.5, 0.5),
                                     "impact:d2" = dist_normal(5, 1)))
  for (smp in c("random", "lhs")) {
    est <- vapply(1:30, function(r) {
      s <- deconflate(m, n_draws = 80, sampling = smp, lhs_replicates = 8, seed = 1000 + r)$draws$summary
      c(mean = qrow(s, "adjusted:d1")$mean, mcse = qrow(s, "adjusted:d1")$mcse)
    }, c(mean = 0, mcse = 0))
    ratio <- stats::sd(est["mean", ]) / mean(est["mcse", ])
    expect_gt(ratio, 0.5, label = paste("sd(run means) / mean(mcse),", smp))
    expect_lt(ratio, 2, label = paste("sd(run means) / mean(mcse),", smp))
  }
})

test_that("the Monte Carlo functions of 0.3 are gone", {
  ns <- asNamespace("deconflate")
  for (f in c("cm_monte_carlo", "cm_sampler", "cm_batch_sampler", "cm_reweight", "cm_scenario",
              "cm_diagnose", "cm_suggest_proposal", "sampler_global_dairy", "cm_mc_gap")) {
    expect_false(f %in% getNamespaceExports(ns), info = f)
  }
  expect_false(exists("cm_monte_carlo", envir = ns, inherits = FALSE))
  expect_false(exists("cm_reweight", envir = ns, inherits = FALSE))
})
