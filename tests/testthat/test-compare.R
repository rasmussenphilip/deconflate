# Method comparison tables (R/compare.R) and the Monte Carlo stability tools
# they summarise. Reference values from inst/validation/reference_v02.py and
# reference_v020_tests.py (supplement impacts in percent: 2.5, 5, 7.5).

supp_valuation <- list(observed = 10000, direction = "decrease", effect = "percent", unit_value = 0.3)

test_that("compare_methods tabulates the methods for a model, with valuation", {
  cmp <- compare_methods(example_supplement(), valuation = supp_valuation)
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_equal(cmp$source, "model")
  expect_equal(cmp$units, "%")
  expect_equal(cmp$label, "yield")
  expect_length(cmp$failed, 0)
  expect_equal(names(cmp$impacts), c("disease", "raw", "published", "simultaneous", "global"))
  expect_equal(cmp$impacts$raw, c(2.5, 5, 7.5))
  expect_equal(cmp$impacts$published, c(2.06500099913358, 3.69929353514358, 6.74847307107448),
               tolerance = 1e-10)
  expect_equal(cmp$impacts$simultaneous, c(2.14325047838512, 3.38707958918446, 6.93420945118756),
               tolerance = 1e-10)
  expect_equal(cmp$impacts$global, cmp$impacts$simultaneous, tolerance = 1e-8)
  expect_equal(cmp$change$published, cmp$impacts$published / cmp$impacts$raw - 1)
  expect_equal(nrow(cmp$long), 9)
  expect_true(all(c("disease", "method", "raw", "adjusted", "change", "contribution", "sign_change") %in%
                    names(cmp$long)))
  expect_false(any(cmp$long$sign_change))
  expect_equal(nrow(cmp$diagnostics), 3)
  expect_equal(cmp$diagnostics$method, cmp$methods)
  expect_equal(names(cmp$results), cmp$methods)

  tt <- cmp$totals
  expect_equal(names(tt), c("method", "raw_sum", "adjusted_total", "gap", "value"))
  expect_equal(tt$method, cmp$methods)
  expect_equal(tt$raw_sum, rep(2.5, 3))
  expect_equal(tt$adjusted_total[1:2], c(2.1110887443997908, 2.109228876453694), tolerance = 1e-10)
  expect_equal(tt$gap[1:2], c(215.6616839763883, 215.46759232202567), tolerance = 1e-9)
  expect_equal(tt$value, 0.3 * tt$gap)
  expect_equal(tt$gap[3], tt$gap[2], tolerance = 1e-8)
  # Without valuation there is no gap or value.
  expect_equal(names(compare_methods(example_supplement())$totals), c("method", "raw_sum", "adjusted_total"))
  expect_equal(compare_methods(example_supplement(), methods = "published")$methods, "published")
  expect_output(print(cmp), "Adjusted values")
  expect_output(print(cmp), "Totals")
})

test_that("methods that cannot be run are reported in `failed`", {
  m_adj <- supp_model(c(2.2, 5, 7.5), estimand = c("adjusted_linear", "crude", "crude"),
                      adjusted_for = c("d2", NA, NA))
  cmp <- compare_methods(m_adj)
  expect_equal(cmp$methods, c("simultaneous", "global"))
  expect_equal(names(cmp$failed), "published")
  expect_match(cmp$failed[["published"]], "crude estimates only")
  expect_equal(names(cmp$impacts), c("disease", "raw", "simultaneous", "global"))
  expect_output(print(cmp), "Not run")

  # Interactions need the global method.
  m_int <- supp_model(c(2.5, 5, 7.5), interactions = cm_interactions("d1", "d2", 1))
  cmp2 <- compare_methods(m_int)
  expect_equal(cmp2$methods, "global")
  expect_equal(names(cmp2$failed), c("published", "simultaneous"))
  expect_match(cmp2$failed[["simultaneous"]], "require method = 'global'")
  expect_error(compare_methods(m_int, methods = c("published", "simultaneous")), "All methods failed")
})

test_that("sign changes are flagged in the long table and printed", {
  cmp <- compare_methods(supp_model(c(0.1, 5, 7.5)))
  sc <- cmp$long[cmp$long$sign_change, ]
  expect_equal(unique(sc$disease), "d1")
  expect_false("published" %in% sc$method)
  expect_true("simultaneous" %in% sc$method)
  expect_lt(cmp$impacts$simultaneous[1], 0)
  expect_gt(cmp$impacts$published[1], 0)
  expect_output(print(cmp), "Sign changes")
})

test_that("compare_methods stacks the analyses of a cm_analyses object", {
  a <- cm_analyses(supp_population(),
                   yield = cm_impacts(ids3, c(2.5, 5, 7.5), units = "%"),
                   fertility = cm_impacts(ids3, c(1, 2, 0), units = "%"))
  val <- list(yield = supp_valuation,
              fertility = list(observed = 400, direction = "increase", effect = "percent", unit_value = 4))
  cmp <- compare_methods(a, valuation = val)
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$source, "analyses")
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_equal(names(cmp$impacts), c("analysis", "disease", "raw", "published", "simultaneous", "global"))
  expect_equal(cmp$impacts$analysis, rep(c("yield", "fertility"), each = 3))
  expect_equal(nrow(cmp$long), 2 * 3 * 3)
  expect_equal(names(cmp$comparisons), c("yield", "fertility"))
  expect_length(cmp$failed, 0)
  tt <- cmp$totals
  expect_equal(names(tt), c("analysis", "method", "raw_sum", "adjusted_total", "gap", "value"))
  expect_equal(nrow(tt), 6)
  y <- compare_methods(a$models$yield, valuation = val$yield)
  expect_equal(tt$gap[tt$analysis == "yield"], y$totals$gap)
  expect_equal(tt$gap[tt$analysis == "yield" & tt$method == "simultaneous"], 215.46759232202567,
               tolerance = 1e-9)
  f <- tt[tt$analysis == "fertility", ]
  expect_equal(f$value, 4 * f$gap)
  expect_true(all(f$gap > 0))
  # A zero raw impact has no relative change.
  expect_true(is.na(cmp$change$published[cmp$change$analysis == "fertility" & cmp$change$disease == "d3"]))
  expect_output(print(cmp), "Adjusted values")
})

test_that("valuations are checked, and one that fails for a method is reported", {
  expect_error(compare_methods(example_supplement(), valuation = list(observed = 10000)), "direction")
  # Impacts in percent read as proportions: an aggregate loss above 100%.
  bad <- list(observed = 10000, direction = "decrease", effect = "proportion")
  cmp <- compare_methods(example_supplement(), valuation = bad)
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_equal(names(cmp$failed), cmp$methods)
  expect_match(cmp$failed[["simultaneous"]], "^valuation: ")
  expect_equal(names(cmp$totals), c("method", "raw_sum", "adjusted_total", "gap"))
  expect_true(all(is.na(cmp$totals$gap)))
  expect_equal(cmp$impacts, compare_methods(example_supplement())$impacts)
})

test_that("stacked analyses fill columns that only some analyses have", {
  a <- cm_analyses(supp_population(),
                   yield = cm_impacts(ids3, c(2.5, 5, 7.5)),
                   fertility = cm_impacts(ids3, c(1, 2, 0.5), estimand = c("adjusted_linear", "crude", "crude"),
                                          adjusted_for = c("d2", NA, NA)))
  # Valuation for one analysis only; the published method fails for the other.
  cmp <- compare_methods(a, valuation = list(yield = supp_valuation))
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_equal(names(cmp$failed), "fertility: published")
  expect_equal(names(cmp$impacts), c("analysis", "disease", "raw", "published", "simultaneous", "global"))
  expect_true(all(is.na(cmp$impacts$published[cmp$impacts$analysis == "fertility"])))
  expect_true(all(is.finite(cmp$impacts$published[cmp$impacts$analysis == "yield"])))
  expect_equal(nrow(cmp$long), 3 * 3 + 2 * 3)
  tt <- cmp$totals
  expect_equal(names(tt), c("analysis", "method", "raw_sum", "adjusted_total", "gap", "value"))
  expect_equal(nrow(tt), 5)
  expect_true(all(is.na(tt$gap[tt$analysis == "fertility"])))
  expect_equal(tt$gap[tt$analysis == "yield" & tt$method == "simultaneous"], 215.46759232202567,
               tolerance = 1e-9)
  expect_output(print(cmp), "fertility: published")
  # A single valuation list (not one per analysis) is an error, not ignored.
  expect_error(compare_methods(a, valuation = supp_valuation), "named after the analyses")
})

test_that("compare_methods reproduces the UK dairy totals (Rasmussen et al. 2022)", {
  skip_on_cran()
  eco <- uk_dairy_2022_economics()
  cmp <- compare_methods(example_uk_dairy_2022(), valuation = eco$valuation)
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_equal(nrow(cmp$impacts), 3 * 13)
  expect_equal(nrow(cmp$long), 3 * 3 * 13)
  expect_equal(cmp$impacts$global, cmp$impacts$simultaneous, tolerance = 1e-6)
  tt <- cmp$totals
  total <- vapply(c("published", "simultaneous"), function(m) sum(tt$value[tt$method == m]), numeric(1))
  expect_equal(unname(total) + sum(eco$additional), c(402.2475530293503, 393.1318185491835),
               tolerance = 1e-9)
})

test_that("compare_methods compares the hazard-ratio methods", {
  cmp <- compare_methods(supp_hr_model(), overall_risk = 0.25)
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$source, "hr")
  expect_equal(cmp$units, "hazard ratio")
  expect_equal(cmp$methods, c("published", "first_order", "snapshot"))
  expect_equal(names(cmp$impacts), c("disease", "raw", "published", "first_order", "snapshot"))
  expect_equal(cmp$impacts$published, c(1.4130001998, 1.9090403704, 1.1927014244), tolerance = 1e-8)
  expect_equal(cmp$impacts$first_order, c(1.4029251095, 1.8874397472, 1.1691237145), tolerance = 1e-8)
  expect_equal(cmp$impacts$snapshot, c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 1e-7)
  expect_equal(cmp$change$snapshot, cmp$impacts$snapshot / cmp$impacts$raw - 1)
  tt <- cmp$totals
  expect_equal(names(tt), c("method", "overall_risk", "disease_free_risk", "attributable"))
  expect_equal(tt$attributable[tt$method == "snapshot"], 0.03664898509534775, tolerance = 1e-7)
  expect_equal(tt$attributable[tt$method == "snapshot"] + tt$disease_free_risk[tt$method == "snapshot"], 0.25)
  expect_equal(tt$attributable,
               vapply(cmp$results, function(r) attributable_risk(r, 0.25, allocate = FALSE)$summary$attributable,
                      numeric(1), USE.NAMES = FALSE),
               tolerance = 1e-10)
  expect_output(print(cmp), "Totals")
  expect_null(compare_methods(supp_hr_model())$totals)
  expect_output(print(compare_methods(supp_hr_model())), "hazard ratio")

  # The published approach cannot handle adjusted hazard ratios.
  hr_adj <- cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = c("adjusted", "crude", "crude"),
                             adjusted_for = c("d2", NA, NA))
  cmp2 <- compare_methods(cm_hr_model(supp_population(), hr_adj), overall_risk = 0.25)
  expect_equal(cmp2$methods, c("first_order", "snapshot"))
  expect_equal(names(cmp2$failed), "published")
  expect_equal(nrow(cmp2$totals), 2)
  expect_output(print(cmp2), "Not run")

  # A method whose adjusted hazard ratios cannot be valued keeps its column
  # and gets NA attributable risk, with the reason in `failed`.
  cmp3 <- compare_methods(cm_hr_model(supp_population(), cm_hazard_ratios(ids3, c(0.5, 4, 1))),
                          overall_risk = 0.25)
  expect_equal(cmp3$methods, c("published", "first_order", "snapshot"))
  expect_equal(names(cmp3$failed), "published")
  expect_match(cmp3$failed[["published"]], "^attributable risk: ")
  expect_true(is.na(cmp3$totals$attributable[cmp3$totals$method == "published"]))
  expect_true(all(is.finite(cmp3$totals$attributable[cmp3$totals$method != "published"])))
  expect_error(compare_methods(supp_hr_model(), overall_risk = 1.5), "between 0 and 1")
})

test_that("Monte Carlo with several methods compares them on the same draws", {
  s <- cm_sampler(example_supplement(),
                  impacts = list(d1 = dist_normal(2.5, 0.5), d3 = dist_normal(7.5, 1)))
  mc <- cm_monte_carlo(s, 200, method = c("published", "simultaneous"), seed = 3)
  expect_equal(mc$method, c("published", "simultaneous"))
  expect_equal(mc$n_rejected, 0)
  expect_equal(nrow(mc$draws), 200 * 2 * 3)
  sm <- summary(mc, diagnose = FALSE)
  expect_equal(nrow(sm), 6)
  expect_true(all(c("method", "trimmed_mean", "rel_mcse", "tail_share", "stability") %in% names(sm)))
  cmp <- compare_methods(mc)
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$source, "mc")
  expect_equal(names(cmp$impacts), c("disease", "raw_mean", "published", "simultaneous"))
  expect_equal(names(cmp$stability), c("disease", "published", "simultaneous"))
  expect_true(all(c("quantity", "method", "mean", "mcse", "stability") %in% names(cmp$totals)))
  # The exact method is linear in the impacts (associations fixed), so its
  # mean equals the adjustment of the mean raw impacts.
  A <- deconflate(example_supplement())$conflation$A
  expect_equal(cmp$impacts$simultaneous, as.vector(solve(A, cmp$impacts$raw_mean)), tolerance = 1e-10)
  expect_equal(cmp$impacts$raw_mean[2], 5)
  expect_equal(compare_methods(mc, stat = "median")$impacts$simultaneous,
               sm$q0.5[sm$method == "simultaneous"])
  expect_output(print(cmp), "Monte Carlo")
  expect_error(compare_methods(cm_monte_carlo(s, 10, seed = 1)), "one method")
})

test_that("a published estimate with a pole inside the inputs is flagged", {
  # The published denominator of d2 is m2 + 1.758 (reference_v020_tests.py), so
  # m2 ~ N(0.2, 1.5) reaches both sides of the pole.
  s <- cm_sampler(example_supplement(), impacts = list(d2 = dist_normal(0.2, 1.5)))
  mc <- cm_monte_carlo(s, 400, method = c("published", "simultaneous"), seed = 4)
  expect_message(sm <- summary(mc), "unstable")
  expect_equal(sm$stability[sm$disease == "d2" & sm$method == "published"], "possible_pole")
  dg <- cm_diagnose(mc)
  expect_s3_class(dg, "cm_diagnosis")
  expect_true(any(dg$stability == "possible_pole" & dg$method == "published" & dg$item == "d2"))
  expect_false(any(dg$stability == "possible_pole" & dg$method == "simultaneous"))
  expect_output(print(dg), "Suggestion")
  expect_silent(summary(mc, diagnose = FALSE))
})

test_that("importance sampling with a suggested proposal is unbiased", {
  s <- cm_sampler(example_supplement(), impacts = list(d2 = dist_normal(0.2, 1.5)))
  mc <- cm_monte_carlo(s, 200, method = "simultaneous", seed = 5)
  expect_message(prop <- cm_suggest_proposal(mc, "d2"), "Proposal for impact:d2")
  expect_equal(names(prop), "impact:d2")
  expect_s3_class(prop[[1]], "cm_dist")
  mc_is <- cm_monte_carlo(s, 400, method = "simultaneous", proposal = prop, seed = 6)
  expect_lt(mc_is$ess, 400)
  expect_gt(mc_is$ess, 100)
  expect_output(print(mc_is), "importance sampling")
  # Target: the exact method is linear, so the true mean is the adjustment
  # of the mean raw impacts (2.5, 0.2, 7.5).
  A <- deconflate(example_supplement())$conflation$A
  target <- solve(A, c(2.5, 0.2, 7.5))[2]
  est <- summary(mc_is, diagnose = FALSE)
  expect_lt(abs(est$mean[est$disease == "d2"] - target), 0.5)
  expect_error(cm_monte_carlo(s, 5, proposal = list("impact:zz" = dist_fixed(1))), "No sampled input")
})

test_that("Latin hypercube sampling stratifies each input within replicate blocks", {
  s <- cm_sampler(example_supplement(), impacts = list(d1 = dist_uniform(2, 3)))
  mc <- cm_monte_carlo(s, 50, sampling = "lhs", seed = 7)
  expect_equal(mc$sampling, "lhs")
  expect_equal(mc$n_rejected, 0)
  expect_equal(sort(unique(mc$block)), 1:10)
  x <- mc$params[["impact:d1"]]
  for (b in 1:10) {
    expect_equal(sort(floor((x[mc$block == b] - 2) * 5)), 0:4)
  }
  expect_output(print(mc), "Latin hypercube")
  expect_error(cm_monte_carlo(function(i) example_supplement(), 5, sampling = "lhs"), "cm_sampler")
})
