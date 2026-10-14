# Contribution tables, summaries and plots (R/report.R, R/plot.R).
# Supplement impacts are in percent (2.5, 5, 7.5); reference values from
# inst/validation/reference_v020_tests.py.

test_that("contribution tables have the documented columns and add up", {
  res <- deconflate(example_supplement())
  ct <- contribution_table(res)
  expect_equal(names(ct), c("disease", "raw", "adjusted", "change", "contribution", "main",
                            "interaction", "share"))
  expect_equal(ct$disease, ids3)
  expect_equal(ct$adjusted, c(2.14325047838512, 3.38707958918446, 6.93420945118756), tolerance = 1e-10)
  expect_equal(ct$contribution, ct$adjusted * c(0.10, 0.15, 0.20))
  expect_equal(ct$interaction, c(0, 0, 0))
  expect_equal(sum(ct$contribution), res$totals$adjusted_total)
  expect_equal(sum(ct$share), 1)
  # Valuation is no longer part of the package.
  expect_equal(names(formals(contribution_table)), "result")
  expect_error(contribution_table(res, list(observed = 10000)), "unused argument")
  expect_error(contribution_table(example_supplement()), "deconflate")
})

test_that("contribution tables of event results and of results with draws", {
  ev <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.25)
  ct <- contribution_table(ev)
  expect_equal(names(ct), c("disease", "measure", "raw", "adjusted_hr", "change", "risk_alone", "rr_alone",
                            "rd_alone", "attributable", "main", "interaction", "share"))
  expect_equal(ct$adjusted_hr, c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 1e-7)
  expect_equal(sum(ct$attributable), 0.03664898509534775, tolerance = 1e-7)
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impact:d1" = dist_normal(2.5, 0.5)))
  ctd <- contribution_table(deconflate(m, n_draws = 50, seed = 1))
  expect_true(all(c("adjusted_lower", "adjusted_upper", "contribution_lower", "contribution_upper") %in% names(ctd)))
  expect_true(all(ctd$adjusted_lower <= ctd$adjusted & ctd$adjusted <= ctd$adjusted_upper))
})

test_that("zero aggregates give NA shares", {
  # All impacts zero.
  res0 <- deconflate(supp_model(c(0, 0, 0)))
  expect_equal(res0$totals$adjusted_total, 0)
  expect_equal(res0$adjusted$adjusted, c(0, 0, 0))
  expect_equal(res0$contributions$total, c(0, 0, 0))
  expect_true(all(is.na(res0$contributions$share)))
  expect_true(all(is.na(contribution_table(res0)$share)))
  # A zero naive sum has no defined reduction.
  expect_true(is.na(summary(res0)$totals$reduction))

  # Impacts that cancel: two independent diseases (an odds ratio of 1) with
  # equal probability and impacts +1 and -1.
  m <- cm_model(cm_population(cm_diseases(c("a", "b"), c(0.2, 0.2)), cm_associations("a", "b", 1)),
                cm_impacts(c("a", "b"), c(1, -1)))
  res <- deconflate(m)
  expect_equal(res$adjusted$adjusted, c(1, -1))
  expect_equal(res$totals$adjusted_total, 0)
  expect_equal(res$contributions$total, c(0.2, -0.2))
  expect_true(all(is.na(res$contributions$share)))
  s <- summary(res)
  expect_equal(s$totals$raw_total, 0)
  expect_equal(s$totals$adjusted_total, 0)
  expect_true(is.na(s$totals$reduction))
})

test_that("summary() reports totals and contributions", {
  res <- deconflate(example_supplement())
  s <- summary(res)
  expect_s3_class(s, "summary.cm_result")
  expect_equal(s$method, "pairwise shortcut")
  expect_equal(s$label, "yield")
  expect_equal(s$units, "%")
  expect_equal(names(s$totals), c("raw_total", "adjusted_total", "interaction_total", "reduction"))
  expect_equal(s$totals$raw_total, 2.5)
  total <- sum(c(0.10, 0.15, 0.20) * c(2.14325047838512, 3.38707958918446, 6.93420945118756))
  expect_equal(s$totals$adjusted_total, total, tolerance = 1e-10)
  expect_equal(s$totals$reduction, 1 - total / 2.5, tolerance = 1e-10)
  expect_equal(s$contributions, contribution_table(res))
  expect_false("valuation" %in% names(formals(summary.cm_result)))
  expect_output(print(s), "Comorbidity adjustment \\(pairwise shortcut\\) - yield \\[%\\]")
  expect_output(print(s), "Feasibility")
  # Sign changes are printed.
  ss <- summary(deconflate(supp_model(c(0.1, 5, 7.5)), warn = FALSE))
  expect_output(print(ss), "Sign changes: d1")
  # Event results.
  se <- summary(deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.25))
  expect_s3_class(se, "summary.cm_result")
  expect_equal(se$method, "snapshot")
  expect_equal(se$totals$attributable, 0.03664898509534775, tolerance = 1e-7)
  expect_output(print(se), "Comorbidity adjustment \\(snapshot\\)")
  # Notes are carried over.
  su <- summary(deconflate(cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5),
                                                                          adjusted_for = c("parity", NA, NA)))))
  expect_true(any(grepl("factors outside the model", su$notes)))
  expect_output(print(su), "Note: Impact row 1 is adjusted only")
})

test_that("plots draw without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(deconflate(example_supplement())))
  expect_invisible(plot(deconflate(read_five("impacts_welfare.csv"), n_draws = 0)))
  expect_invisible(plot(screen_associations(example_supplement())))
  expect_invisible(plot(sensitivity_oat(example_supplement())))
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impact:d1" = dist_normal(2.5, 0.5)))
  expect_invisible(plot(deconflate(m, n_draws = 20, seed = 1)))
  ev <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = dist_beta(25, 75), n_draws = 20, seed = 1)
  expect_invisible(plot(ev))
})
