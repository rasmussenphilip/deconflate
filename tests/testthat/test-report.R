# Contribution tables, summaries and plots (R/report.R, R/plot.R).
# Supplement impacts are in percent (2.5, 5, 7.5); reference values from
# inst/validation/reference_v020_tests.py.

test_that("contribution tables have the documented columns and add up", {
  res <- deconflate(example_supplement())
  ct <- contribution_table(res)
  expect_equal(names(ct), c("disease", "raw", "adjusted", "adjusted_lower", "adjusted_upper", "change",
                            "raw_contribution", "contribution", "contribution_lower", "contribution_upper",
                            "main", "interaction", "share"))
  # Raw contributions: p * raw here (single-disease rows, diseases unrelated).
  expect_equal(ct$raw_contribution, c(2.5, 5, 7.5) * c(0.10, 0.15, 0.20))
  # Without draws, the interval columns are there but empty.
  expect_true(all(is.na(ct$adjusted_lower) & is.na(ct$contribution_upper)))
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
  expect_equal(names(ct), c("disease", "measure", "raw", "adjusted", "adjusted_lower", "adjusted_upper",
                            "change", "risk_alone", "rr_alone", "rd_alone", "raw_contribution", "contribution",
                            "contribution_lower", "contribution_upper", "main", "interaction", "share"))
  expect_equal(ct$adjusted, c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 1e-7)
  expect_equal(sum(ct$contribution), 0.03664898509534775, tolerance = 1e-7)
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
  expect_equal(res0$contributions$contribution, c(0, 0, 0))
  expect_true(all(is.na(res0$contributions$share)))
  expect_true(all(is.na(contribution_table(res0)$share)))
  # A zero naive sum has no defined reduction.
  expect_true(is.na(total_value(summary(res0), "reduction")))

  # Impacts that cancel: two independent diseases (an odds ratio of 1) with
  # equal probability and impacts +1 and -1.
  m <- cm_model(cm_population(cm_diseases(c("a", "b"), c(0.2, 0.2)), cm_associations("a", "b", 1)),
                cm_impacts(c("a", "b"), c(1, -1)))
  res <- deconflate(m)
  expect_equal(res$adjusted$adjusted, c(1, -1))
  expect_equal(res$totals$adjusted_total, 0)
  expect_equal(res$contributions$contribution, c(0.2, -0.2))
  expect_true(all(is.na(res$contributions$share)))
  s <- summary(res)
  expect_equal(total_value(s, "raw total"), 0)
  expect_equal(total_value(s, "adjusted total"), 0)
  expect_true(is.na(total_value(s, "reduction")))
})

test_that("summary() reports totals and contributions", {
  res <- deconflate(example_supplement())
  s <- summary(res)
  expect_s3_class(s, "summary.cm_result")
  expect_equal(s$method, "pairwise shortcut")
  expect_equal(s$outcome, "yield")
  expect_equal(s$units, "%")
  expect_equal(names(s$totals), c("quantity", "value", "lower", "upper"))
  expect_equal(s$totals$quantity, c("raw total", "adjusted total", "interactions", "reduction"))
  expect_equal(total_value(s, "raw total"), 2.5)
  total <- sum(c(0.10, 0.15, 0.20) * c(2.14325047838512, 3.38707958918446, 6.93420945118756))
  expect_equal(total_value(s, "adjusted total"), total, tolerance = 1e-10)
  expect_equal(total_value(s, "reduction"), 1 - total / 2.5, tolerance = 1e-10)
  expect_equal(s$contributions$contribution, contribution_table(res)$contribution)
  expect_equal(names(s$contributions), c("disease", "contribution", "main", "interaction", "share"))
  expect_equal(names(s$impacts), c("disease", "raw", "adjusted", "change"))
  expect_null(s$alone)
  expect_false("valuation" %in% names(formals(summary.cm_result)))
  expect_output(print(s), "Comorbidity adjustment \\(pairwise shortcut\\) - yield \\[%\\]")
  expect_output(print(s), "Feasibility")
  out <- capture.output(print(s))
  expect_true(all(c("Impacts", "Contributions to the adjusted total (Shapley values)", "Totals", "Diagnostics")
                  %in% out))
  expect_true(all(nchar(out) <= 80))
  # Sign changes are printed, and print() points to them.
  rs <- deconflate(supp_model(c(0.1, 5, 7.5)), warn = FALSE)
  expect_match(printed(summary(rs)), "sign changes 1 (d1)", fixed = TRUE)
  expect_match(printed(rs), "Check: sign changes for d1; see summary().", fixed = TRUE)
  expect_false(grepl("Check:", printed(res)))
  # Event results.
  se <- summary(deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.25))
  expect_s3_class(se, "summary.cm_result")
  expect_equal(se$method, "snapshot")
  expect_equal(total_value(se, "attributable risk"), 0.03664898509534775, tolerance = 1e-7)
  expect_equal(names(se$alone), c("disease", "risk_alone", "rr_alone", "rd_alone"))
  expect_equal(names(se$impacts), c("disease", "measure", "raw", "adjusted", "change"))
  expect_output(print(se), "Comorbidity adjustment \\(snapshot\\)")
  expect_output(print(se), "Contributions to the attributable risk")
  # Notes are carried over.
  su <- summary(deconflate(cm_model(supp_unknown_population(), cm_impacts(ids3, c(2.5, 5, 7.5),
                                                                          adjusted_for = c("parity", NA, NA)))))
  expect_true(any(grepl("factors outside the model", su$notes)))
  expect_match(printed(su), "* Impact row 1 is adjusted only", fixed = TRUE)
})

test_that("plots draw without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(deconflate(example_supplement())))
  expect_invisible(plot(deconflate(read_five("impacts_welfare.csv"), n_draws = 0)))
  expect_invisible(plot(cm_sensitivity(deconflate(example_supplement(), n_draws = 0))))
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impact:d1" = dist_normal(2.5, 0.5)))
  expect_invisible(plot(deconflate(m, n_draws = 20, seed = 1)))
  ev <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = dist_beta(25, 75), n_draws = 20, seed = 1)
  expect_invisible(plot(ev))
})

test_that("printed numbers have 3 significant figures and tables fit in 80 columns", {
  expect_equal(fmt_num(c(4.8, 0.1, -0.1451, 0.00158, 1234.5, 0, NA, 3.87e-05, 2e5, Inf)),
               c("4.80", "0.100", "-0.145", "0.00158", "1230", "0", "", "3.87e-05", "2.00e+05", "Inf"))
  withr_width <- options(width = 80)
  on.exit(options(withr_width))
  m <- read_five("impacts_yield.csv")
  r <- deconflate(m, n_draws = 30, seed = 1)
  out <- capture.output(print(r))
  expect_true(all(nchar(out) <= 80))
  expect_true(any(grepl("^ disease +raw +adjusted +lower +upper +change +contribution +share$", out)))
  expect_true(all(nchar(capture.output(print(summary(r)))) <= 80))
  ev <- cull_five(overall_risk = 0.25)
  out <- capture.output(print(ev))
  expect_true(all(nchar(out) <= 80))
  expect_true(any(grepl("^ disease +measure +raw +adjusted +change +contribution +share$", out)))
  expect_true(all(nchar(capture.output(print(summary(ev)))) <= 80))
})

test_that("event results have a contributions table like additive results", {
  ev <- cull_five(overall_risk = 0.25)
  expect_equal(names(ev$contributions), c("disease", "contribution", "main", "interaction", "share", "raw_contribution"))
  expect_null(ev$attributable$by_disease)
  expect_equal(sum(ev$contributions$contribution), ev$attributable$summary$attributable, tolerance = 1e-10)
  ad <- deconflate(read_five("impacts_yield.csv"), n_draws = 0)
  expect_equal(names(ad$contributions), names(ev$contributions))
})
