# Method comparison tables and Monte Carlo stability tools.

test_that("compare_methods tabulates the methods side by side with totals", {
  skip_on_cran()
  cmp <- compare_methods(example_uk_dairy_2022(), economics = uk_dairy_2022_economics())
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_equal(names(cmp$impacts),
               c("outcome", "disease", "scale", "raw", "published", "simultaneous", "global"))
  expect_equal(cmp$impacts$global, cmp$impacts$simultaneous, tolerance = 1e-6)
  expect_equal(unname(cmp$total_value[c("published", "simultaneous")]),
               c(402.2475530293503, 393.1318185491835), tolerance = 1e-9)
  expect_true(all(c("outcome", "method", "raw_loss", "adjusted_loss", "gap", "value") %in%
                    names(cmp$totals)))
  expect_equal(nrow(cmp$long), 3 * 39)
  expect_output(print(cmp), "Adjusted impacts")
})

test_that("compare_methods values hazard-ratio outcomes and reports failures", {
  eco <- list(observed = c(yield = 10000, culling = 0.25),
              unit_value = c(yield = 0.3, culling = 1000))
  cmp <- compare_methods(supp_hr(), economics = eco)
  g <- cmp$totals[cmp$totals$outcome == "culling" & cmp$totals$method == "global", ]
  expect_equal(g$gap, 0.03664898509534775, tolerance = 1e-7)
  expect_equal(g$value, 1000 * g$gap)
  # Interactions need the global method: the others fail and are reported.
  m <- example_supplement()
  m_int <- cm_model(m$diseases, m$associations, m$impacts,
                    interactions = cm_interactions("d1", "d2", 0.01, outcome = "yield"))
  cmp2 <- compare_methods(m_int)
  expect_equal(cmp2$methods, "global")
  expect_equal(names(cmp2$failed), c("published", "simultaneous"))
  expect_output(print(cmp2), "Failed")
})

test_that("Monte Carlo with several methods compares them on the same draws", {
  s <- cm_sampler(example_supplement(),
                  impacts = list("yield:d1" = dist_normal(2.5, 0.5),
                                 "yield:d3" = dist_normal(7.5, 1)))
  mc <- cm_monte_carlo(s, 200, method = c("published", "simultaneous"), seed = 3)
  expect_equal(mc$method, c("published", "simultaneous"))
  expect_equal(nrow(mc$draws), 200 * 2 * 3)
  sm <- summary(mc, diagnose = FALSE)
  expect_equal(nrow(sm), 6)
  expect_true(all(c("method", "trimmed_mean", "rel_mcse", "tail_share", "stability") %in% names(sm)))
  cmp <- compare_methods(mc)
  expect_equal(names(cmp$impacts),
               c("outcome", "disease", "scale", "raw_mean", "published", "simultaneous"))
  # The exact method is linear in the impacts (associations fixed), so its
  # mean equals the adjustment of the mean raw impacts.
  A <- deconflate(example_supplement())$conflation$yield$A
  expect_equal(cmp$impacts$simultaneous, as.vector(solve(A, cmp$impacts$raw_mean)),
               tolerance = 1e-10)
  expect_output(print(cmp), "Monte Carlo")
  expect_error(compare_methods(cm_monte_carlo(s, 10, seed = 1)), "one method")
})

test_that("a mean that does not exist is flagged, with suggestions", {
  s <- cm_sampler(example_supplement(), impacts = list("yield:d2" = dist_normal(0.2, 1.5)))
  mc <- cm_monte_carlo(s, 400, method = c("published", "simultaneous"), seed = 4)
  expect_message(sm <- summary(mc), "unstable")
  expect_equal(sm$stability[sm$disease == "d2" & sm$method == "published"], "no_mean")
  dg <- cm_diagnose(mc)
  expect_s3_class(dg, "cm_diagnosis")
  expect_true(any(dg$stability == "no_mean" & dg$method == "published"))
  expect_false(any(dg$stability == "no_mean" & dg$method == "simultaneous"))
  expect_output(print(dg), "Suggestion")
  expect_silent(summary(mc, diagnose = FALSE))
})

test_that("importance sampling with a suggested proposal is unbiased", {
  s <- cm_sampler(example_supplement(), impacts = list("yield:d2" = dist_normal(0.2, 1.5)))
  mc <- cm_monte_carlo(s, 200, method = "simultaneous", seed = 5)
  expect_message(prop <- cm_suggest_proposal(mc, "yield", "d2"), "Proposal for impact:yield:d2")
  expect_equal(names(prop), "impact:yield:d2")
  expect_s3_class(prop[[1]], "cm_dist")
  mc_is <- cm_monte_carlo(s, 400, method = "simultaneous", proposal = prop, seed = 6)
  expect_lt(mc_is$ess, 400)
  expect_gt(mc_is$ess, 100)
  expect_output(print(mc_is), "importance sampling")
  # Target: the exact method is linear, so the true mean is the adjustment
  # of the mean raw impacts (0.025, 0.002, 0.075).
  A <- deconflate(example_supplement())$conflation$yield$A
  target <- solve(A, c(0.025, 0.002, 0.075))[2]
  est <- summary(mc_is, diagnose = FALSE)
  expect_lt(abs(est$mean[est$disease == "d2"] - target), 0.005)
  expect_error(cm_monte_carlo(s, 5, proposal = list("impact:yield:zz" = dist_fixed(1))),
               "No sampled input")
})

test_that("Latin hypercube sampling stratifies each input", {
  s <- cm_sampler(example_supplement(), impacts = list("yield:d1" = dist_uniform(2, 3)))
  mc <- cm_monte_carlo(s, 50, sampling = "lhs", seed = 7)
  x <- mc$params[["impact:yield:d1"]]
  expect_equal(sort(floor((x - 2) * 50 + 1e-9)), 0:49)
  expect_equal(mc$sampling, "lhs")
  expect_error(cm_monte_carlo(function(i) example_supplement(), 5, sampling = "lhs"),
               "cm_sampler")
})
