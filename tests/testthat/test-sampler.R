test_that("samplers draw on the input scale and convert", {
  m <- example_global_dairy(inputs = "tables")
  s <- cm_sampler(m, diseases = list(SCK = dist_fixed(0.5)),
                  impacts = list("yield:SCK" = dist_fixed(10)),
                  associations = list("SCK:CK" = dist_fixed(5)))
  d <- s(1)
  expect_equal(d$diseases$prob[d$diseases$id == "SCK"], 1 - exp(-0.5))
  expect_equal(d$impacts$value[d$impacts$outcome == "yield" & d$impacts$disease == "SCK"], 0.10)
  a <- d$associations
  expect_equal(a$value[pair_key(a$disease1, a$disease2) == pair_key("CK", "SCK")], 5)
  expect_true(all(c("prob:SCK", "assoc:CK:SCK", "impact:yield:SCK") %in% names(attr(s, "specs"))))
  expect_error(cm_sampler(m, impacts = list("yield:XX" = dist_fixed(1))), "does not match")
})

test_that("outcome correlation induces dependence but keeps marginals", {
  m <- example_supplement()
  m$impacts <- combine_impacts(m$impacts,
                               cm_impacts(c("d1", "d2", "d3"), c(1, 2, 3), outcome = "fertility",
                                          scale = "percent", direction = "increase"))
  R <- matrix(c(1, 0.9, 0.9, 1), 2, dimnames = list(c("yield", "fertility"), c("yield", "fertility")))
  s <- cm_sampler(m, impacts = list("yield:d1" = dist_normal(2.5, 0.5),
                                    "fertility:d1" = dist_lognormal(log(1), 0.3)),
                  outcome_correlation = R)
  set.seed(11)
  x <- t(vapply(1:2000, function(i) {
    d <- s(i)$impacts
    c(d$value[d$outcome == "yield" & d$disease == "d1"],
      d$value[d$outcome == "fertility" & d$disease == "d1"])
  }, numeric(2)))
  expect_gt(stats::cor(x[, 1], x[, 2], method = "spearman"), 0.85)
  expect_equal(mean(x[, 1]), 0.025, tolerance = 0.02)
  expect_equal(length(attr(s, "correlated")), 2)
})

test_that("Monte Carlo with economics, rejections and scenarios", {
  m <- example_supplement()
  s <- cm_sampler(m,
                  diseases = list(d1 = dist_normal(0.1, 0.06)),  # sometimes negative: rejected
                  associations = list("d1:d2" = dist_mixture(dist_lognormal_ci(2, 1.5, 2.7),
                                                             dist_lognormal(log(2), 0.8),
                                                             weights = c(0.7, 0.3))),
                  impacts = list("yield:d3" = dist_normal(7.5, 1)))
  mc <- cm_monte_carlo(s, 300, economics = list(observed = c(yield = 10000),
                                                unit_value = list(yield = dist_uniform(0.25, 0.35))),
                       seed = 5)
  expect_gt(mc$n_rejected, 0)
  expect_equal(mc$n_rejected + nrow(mc$params), 300)
  expect_true(all(c("prob:d1", "assoc:d1:d2", "impact:yield:d3", "unit_value:yield") %in% names(mc$params)))
  tot <- summary(mc, what = "total")
  expect_true("total" %in% tot$outcome)
  expect_gt(tot$mean[tot$outcome == "total"], 0)
  expect_gt(nrow(summary(mc, what = "rejections")), 0)

  sc <- suppressWarnings(cm_scenario(mc, list("assoc:d1:d2" = dist_lognormal_ci(3, 2.2, 4.1))))
  expect_lt(sc$ess, nrow(mc$params))
  d1 <- function(x) { s <- summary(x); s$mean[s$disease == "d1"] }
  expect_lt(d1(sc), d1(mc))
  expect_error(cm_scenario(mc, list("assoc:d1:zz" = dist_fixed(1))), "No sampled input")
})

test_that("Monte Carlo of the 2024 analysis inputs reproduces Table 5", {
  skip_on_cran()
  mc <- cm_monte_carlo(sampler_global_dairy(), 300, method = "published", seed = 2024)
  s <- summary(mc)
  # Disease probabilities are fixed in the analysis inputs.
  expect_length(unique(mc$params[["prob:SCK"]]), 1)
  expect_false(any(startsWith(names(attr(sampler_global_dairy(), "specs")), "prob:")))
  # Negative normal draws of odds ratios are truncated, not rejected.
  expect_equal(mc$n_rejected, 0)
  # DA (and DYS for culling) are excluded: their raw impacts can be near zero
  # or negative, which makes the published approximation heavy-tailed and its
  # mean unstable in small runs. Long-run means (200,000 draws,
  # inst/validation/reference_2024_analysis.py): yield OC 2.51, PTB 3.32;
  # culling OC 1.48; all others within 0.02 of Table 5.
  y <- s[s$outcome == "yield" & s$disease != "DA", ]
  table5_yield <- c(CK = 0.03, CM = 1.36, DYS = 3.48, LAM = 2.62, MET = 2.87, MF = 0.07,
                    OC = 2.59, PTB = 3.37, RP = 2.30, SCK = 7.11, SCM = 5.58)
  expect_lt(max(abs(100 * y$mean - table5_yield[y$disease])), 0.35)
  cu <- s[s$outcome == "culling" & !s$disease %in% c("DA", "DYS"), ]
  table5_cull <- c(CK = 1.18, CM = 1.90, LAM = 1.40, MET = 1.03, MF = 2.64, OC = 1.51,
                   PTB = 2.07, RP = 1.29, SCK = 1.67, SCM = 1.25)
  expect_lt(max(abs(1 + cu$mean - table5_cull[cu$disease])), 0.12)
})

test_that("the Table 2-4 sampler draws incidence and shifts hazard ratios", {
  s <- sampler_global_dairy(inputs = "tables")
  specs <- attr(s, "specs")
  expect_true("prob:SCK" %in% names(specs))
  expect_equal(specs[["impact:culling:CK"]]$params$sd, 0.30)
  expect_equal(specs[["impact:culling:MET"]]$params$mean, 0.05)
  s2 <- sampler_global_dairy(culling = FALSE)
  expect_false(any(startsWith(names(attr(s2, "specs")), "impact:culling")))
})
