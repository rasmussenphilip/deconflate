# Method comparison tables (R/compare.R). Reference values from
# inst/validation/reference_v02.py, reference_v020_tests.py and
# reference_v040.py (supplement impacts in percent: 2.5, 5, 7.5).

test_that("compare_methods tabulates the methods for additive impacts", {
  cmp <- compare_methods(example_supplement())
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$methods, c("published", "simultaneous", "global"))
  expect_false(cmp$event_model)
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
  expect_null(cmp$draws)

  # Totals: naive and adjusted aggregates only (no gap or value).
  tt <- cmp$totals
  expect_equal(names(tt), c("method", "raw_sum", "adjusted_total"))
  expect_equal(tt$method, cmp$methods)
  expect_equal(tt$raw_sum, rep(2.5, 3))
  expect_equal(tt$adjusted_total[1:2], c(2.1110887443997908, 2.109228876453694), tolerance = 1e-10)
  expect_equal(tt$adjusted_total[3], tt$adjusted_total[2], tolerance = 1e-8)
  expect_equal(compare_methods(example_supplement(), methods = "published")$methods, "published")
  expect_output(print(cmp), "Adjusted values")
  expect_output(print(cmp), "Totals")
})

test_that("the published approximation is not a method of deconflate()", {
  m <- example_supplement()
  expect_error(deconflate(m, method = "published"), class = "deconflate_unsupported")
  expect_error(deconflate(m, method = "published"), "compare_methods")
  expect_equal(deconflate(m)$method, "simultaneous")
  # compare_methods() still includes it by default, through the internal
  # adjust_impacts().
  cmp <- compare_methods(m)
  expect_true("published" %in% cmp$methods)
  expect_equal(cmp$results$published$adjusted$adjusted,
               adjust_impacts(m, method = "published")$adjusted$adjusted)
  expect_error(compare_methods(m, methods = "auto"))
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

  # Interactions need the global method; compare_methods() does not switch.
  m_int <- supp_model(c(2.5, 5, 7.5), interactions = cm_interactions("d1", "d2", 1))
  cmp2 <- compare_methods(m_int)
  expect_equal(cmp2$methods, "global")
  expect_equal(names(cmp2$failed), c("published", "simultaneous"))
  expect_match(cmp2$failed[["simultaneous"]], "need the global method")
  expect_error(compare_methods(m_int, methods = c("published", "simultaneous")), "All methods failed")

  # Unknown pairs: only the global method runs.
  m_unk <- cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5)))
  cmp3 <- compare_methods(m_unk)
  expect_equal(cmp3$methods, "global")
  expect_match(cmp3$failed[["simultaneous"]], "No association for d1-d3")
  expect_equal(cmp3$impacts$global, c(1.9889421686, 3.4043343379, 6.9066256849), tolerance = 1e-7)
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

test_that("compare_methods on the UK dairy example (Rasmussen et al. 2022)", {
  skip_on_cran()
  uk <- example_uk_dairy_2022()
  # Pairs without an association are unknown: only the global method runs.
  cmp <- compare_methods(uk)
  expect_equal(cmp$methods, "global")
  expect_equal(names(cmp$failed), c("published", "simultaneous"))
  expect_equal(nrow(cmp$impacts), 13)
  expect_equal(cmp$impacts$global,
               c(-0.724739, 1.95468, 2.803416, 7.33, 3.28, 4.922689, 3.898285, 1.876286,
                 -0.689519, 4.2, 4.07841, 6.335223, 1.641737), tolerance = 1e-5)
  expect_equal(cmp$totals$adjusted_total, 5.8434734878, tolerance = 1e-7)
  # With the unlisted pairs independent (as in the paper), all three run, and
  # the global method equals the simultaneous one.
  ind <- cm_model(with_independent_pairs(as_population(uk)), uk$impacts)
  cmp2 <- compare_methods(ind)
  expect_equal(cmp2$methods, c("published", "simultaneous", "global"))
  expect_equal(cmp2$impacts$global, cmp2$impacts$simultaneous, tolerance = 1e-6)
  expect_equal(cmp2$totals$adjusted_total[cmp2$totals$method == "published"],
               adjust_impacts(ind, method = "published", warn = FALSE)$totals$adjusted_total, tolerance = 1e-12)
  expect_equal(cmp2$totals$adjusted_total[cmp2$totals$method == "simultaneous"],
               suppressWarnings(deconflate(ind, n_draws = 0))$totals$adjusted_total, tolerance = 1e-12)
})

test_that("compare_methods compares the event-model methods", {
  cmp <- compare_methods(supp_hr_model(), event_model = TRUE, overall_risk = 0.25)
  expect_s3_class(cmp, "cm_comparison")
  expect_true(cmp$event_model)
  expect_equal(cmp$methods, c("published", "first_order", "snapshot"))
  expect_equal(names(cmp$impacts), c("disease", "measure", "raw", "published", "first_order", "snapshot"))
  expect_equal(cmp$impacts$published, c(1.4130001998, 1.9090403704, 1.1927014244), tolerance = 1e-8)
  expect_equal(cmp$impacts$first_order, c(1.4029251095, 1.8874397472, 1.1691237145), tolerance = 1e-8)
  expect_equal(cmp$impacts$snapshot, c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 1e-7)
  expect_equal(cmp$change$snapshot, cmp$impacts$snapshot / cmp$impacts$raw - 1)
  tt <- cmp$totals
  expect_equal(names(tt), c("method", "overall_risk", "disease_free_risk", "attributable"))
  expect_equal(tt$attributable[tt$method == "snapshot"], 0.03664898509534775, tolerance = 1e-7)
  expect_equal(tt$attributable[tt$method == "snapshot"] + tt$disease_free_risk[tt$method == "snapshot"], 0.25)
  expect_equal(tt$attributable[tt$method == "snapshot"],
               deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.25)$attributable$summary$attributable,
               tolerance = 1e-10)
  expect_output(print(cmp), "Totals")
  expect_output(print(cmp), "Adjusted hazard ratios")
  expect_error(compare_methods(supp_hr_model(), event_model = TRUE), "overall_risk")
  expect_error(compare_methods(supp_hr_model()), class = "deconflate_unsupported")

  # The published approach cannot handle adjusted hazard ratios.
  adj <- supp_hr_model(estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                       adjusted_for = c("d2", NA, NA))
  cmp2 <- compare_methods(adj, event_model = TRUE, overall_risk = 0.25)
  expect_equal(cmp2$methods, c("first_order", "snapshot"))
  expect_equal(names(cmp2$failed), "published")
  expect_equal(nrow(cmp2$totals), 2)
  expect_output(print(cmp2), "Not run")

  # Risk-based measures need the snapshot model.
  rr <- supp_hr_model(values = c(1.592, 1.771, 0.0572), measure = c("OR", "RR", "RD"))
  cmp_rr <- compare_methods(rr, event_model = TRUE, overall_risk = 0.25)
  expect_equal(cmp_rr$methods, "snapshot")
  expect_match(cmp_rr$failed[["first_order"]], "hazard or rate ratios")
  expect_equal(cmp_rr$impacts$snapshot, c(1.383112894, 1.8909782708, 1.1439101733), tolerance = 1e-7)

  # A method that gives a non-positive adjusted hazard ratio is undefined: it
  # is reported in `failed` and left out of the comparison.
  cmp3 <- compare_methods(supp_hr_model(c(0.5, 4, 1)), event_model = TRUE, overall_risk = 0.25)
  expect_equal(cmp3$methods, c("first_order", "snapshot"))
  expect_equal(names(cmp3$failed), "published")
  expect_match(cmp3$failed[["published"]], "^undefined")
  expect_true(all(is.finite(cmp3$totals$attributable)))
  expect_error(compare_methods(supp_hr_model(), event_model = TRUE, overall_risk = 1.5), "between 0 and 1")
})

test_that("with draws, the methods are compared on the same draws", {
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impact:d1" = dist_normal(2.5, 0.5), "impact:d3" = dist_normal(7.5, 1)))
  cmp <- compare_methods(m, methods = c("published", "simultaneous"), n_draws = 200, seed = 3)
  dr <- cmp$draws
  expect_equal(dr$n_draws, 200)
  expect_equal(dr$n_rejected, 0)
  expect_true(all(c("method", "quantity", "mean", "lower", "upper", "mcse", "stability") %in% names(dr$summary)))
  expect_setequal(unique(dr$summary$method), c("published", "simultaneous"))
  # The exact method is linear in the impacts (associations fixed), so its
  # mean equals the adjustment of the mean raw impacts of the draws.
  A <- deconflate(example_supplement())$conflation$A
  s <- dr$summary[dr$summary$method == "simultaneous", ]
  raw_mean <- colMeans(dr$params[, c("impact:d1", "impact:d3")])
  expect_equal(s$mean[match(paste0("adjusted:", ids3), s$quantity)],
               as.vector(solve(A, c(raw_mean[1], 5, raw_mean[2]))), tolerance = 1e-8)
  expect_output(print(cmp), "Uncertainty \\(200 draws, 0 rejected\\)")
  # Without distributions, a note.
  expect_match(compare_methods(example_supplement(), n_draws = 10)$notes, "No input has a distribution")
})
