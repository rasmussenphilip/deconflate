# Contribution tables, summaries and plots (R/report.R, R/plot.R), and the
# internal gap and valuation helpers used only to reproduce Rasmussen et al.
# (2022) (legacy_gap() and legacy_value(), R/legacy.R). Supplement impacts
# are in percent (2.5, 5, 7.5); reference values from
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

test_that("legacy gaps follow the direction and effect of the impacts", {
  pub <- adjust_impacts(example_supplement(), method = "published")
  g <- legacy_gap(pub, observed = 10000, direction = "decrease", effect = "percent")
  expect_equal(names(g$summary), c("observed", "disease_free", "gap", "aggregate", "direction", "effect"))
  expect_equal(names(g$attribution), c("disease", "gap", "gap_main", "gap_interaction"))
  # Supplementary File (unrounded): 10,215.7 and 21.1, 56.7, 137.9.
  expect_equal(g$summary$disease_free, 10215.661683976388, tolerance = 1e-10)
  expect_equal(g$summary$gap, 215.6616839763883, tolerance = 1e-9)
  expect_equal(g$summary$aggregate, 2.1110887443997908 / 100, tolerance = 1e-10)
  expect_equal(g$attribution$gap, c(21.09535158422185, 56.68609683712177, 137.88023555504404),
               tolerance = 1e-9)
  expect_equal(g$attribution$gap_interaction, c(0, 0, 0))
  # Percent and proportion impacts give the same gap.
  pub_prop <- adjust_impacts(supp_model(c(0.025, 0.05, 0.075)), method = "published")
  expect_equal(legacy_gap(pub_prop, 10000, "decrease", "proportion")$summary$gap, g$summary$gap,
               tolerance = 1e-10)
  # An outcome that disease increases.
  gi <- legacy_gap(pub, 10000, direction = "increase", effect = "percent")
  expect_equal(gi$summary$disease_free, 9793.255681595545, tolerance = 1e-10)
  expect_equal(gi$summary$gap, 10000 - 9793.255681595545, tolerance = 1e-9)
  expect_equal(sum(gi$attribution$gap), gi$summary$gap, tolerance = 1e-10)
  # Impacts in the outcome's own units.
  ga <- legacy_gap(pub, 10000, direction = "decrease", effect = "absolute")
  expect_equal(ga$summary$disease_free, 10000 + 2.1110887443997908, tolerance = 1e-12)
  expect_equal(ga$attribution$gap, pub$contributions$total)
  gai <- legacy_gap(pub, 10000, direction = "increase", effect = "absolute")
  expect_equal(gai$summary$gap, 2.1110887443997908, tolerance = 1e-10)

  expect_error(legacy_gap(pub, 10000, effect = "proportion"), "below 1")
  expect_error(legacy_gap(deconflate(supp_model(c(-250, -500, -750))), 10000, "increase", "percent"),
               "above -1")
  # The exported helpers were removed with valuation.
  ns <- asNamespace("deconflate")
  expect_false(exists("productivity_gap", envir = ns, inherits = FALSE))
  expect_false(exists("value_losses", envir = ns, inherits = FALSE))
})

test_that("legacy_value values a gap", {
  pub <- adjust_impacts(example_supplement(), method = "published")
  g <- legacy_gap(pub, 10000, "decrease", "percent")
  vl <- legacy_value(g, unit_value = 0.3)
  expect_equal(names(vl), c("by_disease", "value"))
  expect_equal(names(vl$by_disease), c("disease", "gap", "value"))
  expect_equal(vl$value, 0.3 * 215.6616839763883, tolerance = 1e-9)
  expect_equal(vl$by_disease$value, 0.3 * g$attribution$gap)
  expect_equal(sum(vl$by_disease$value), vl$value, tolerance = 1e-10)
})

test_that("zero aggregates give NA shares but valid legacy gaps and values", {
  # All impacts zero.
  res0 <- deconflate(supp_model(c(0, 0, 0)))
  expect_equal(res0$totals$adjusted_total, 0)
  expect_true(all(is.na(res0$contributions$share)))
  g0 <- legacy_gap(res0, 10000, "decrease", "percent")
  expect_equal(g0$summary$gap, 0)
  expect_equal(g0$summary$disease_free, 10000)
  expect_equal(g0$attribution$gap, c(0, 0, 0))
  expect_false(anyNA(unlist(g0$summary[c("observed", "disease_free", "gap", "aggregate")])))
  v0 <- legacy_value(g0, 0.3)
  expect_equal(v0$value, 0)
  expect_equal(v0$by_disease$value, c(0, 0, 0))
  expect_equal(adjust_impacts(supp_model(c(0, 0, 0)), method = "published")$adjusted$adjusted, c(0, 0, 0))
  expect_true(all(is.na(contribution_table(res0)$share)))
  # A zero naive sum has no defined reduction.
  expect_true(is.na(summary(res0)$totals$reduction))

  # Impacts that cancel: two independent diseases with equal probability and
  # impacts +1 and -1.
  m <- cm_model(cm_population(cm_diseases(c("a", "b"), c(0.2, 0.2))),
                cm_impacts(c("a", "b"), c(1, -1)))
  res <- deconflate(m)
  expect_equal(res$adjusted$adjusted, c(1, -1))
  expect_equal(res$totals$adjusted_total, 0)
  expect_true(all(is.na(res$contributions$share)))
  ga <- legacy_gap(res, 100, "decrease", "absolute")
  expect_equal(ga$summary$gap, 0)
  expect_equal(ga$attribution$gap, c(0.2, -0.2))
  expect_true(all(is.finite(ga$attribution$gap)))
  gp <- legacy_gap(res, 100, "decrease", "proportion")
  expect_equal(gp$attribution$gap, c(20, -20))
  expect_equal(sum(gp$attribution$gap), gp$summary$gap)
  expect_equal(legacy_value(ga, 2)$by_disease$value, c(0.4, -0.4))
  expect_equal(legacy_value(ga, 2)$value, 0)
  s <- summary(res)
  expect_equal(s$totals$raw_sum, 0)
  expect_equal(s$totals$adjusted_total, 0)
  expect_true(is.na(s$totals$reduction))
})

test_that("summary() reports totals and contributions", {
  res <- adjust_impacts(example_supplement(), method = "published")
  s <- summary(res)
  expect_s3_class(s, "summary.cm_result")
  expect_equal(s$method, "published")
  expect_equal(s$label, "yield")
  expect_equal(s$units, "%")
  expect_equal(names(s$totals), c("raw_sum", "adjusted_total", "interaction_total", "reduction"))
  expect_equal(s$totals$raw_sum, 2.5)
  expect_equal(s$totals$adjusted_total, 2.1110887443997908, tolerance = 1e-10)
  expect_equal(s$totals$reduction, 1 - 2.1110887443997908 / 2.5, tolerance = 1e-10)
  expect_equal(s$contributions, contribution_table(res))
  expect_false("valuation" %in% names(formals(summary.cm_result)))
  expect_output(print(s), "Comorbidity adjustment \\(method: published\\) - yield \\[%\\]")
  expect_output(print(s), "Feasibility")
  ssim <- summary(deconflate(example_supplement()))
  expect_equal(ssim$method, "simultaneous")
  expect_output(print(ssim), "Comorbidity adjustment \\(method: simultaneous\\) - yield \\[%\\]")
  # Sign changes are printed.
  ss <- summary(deconflate(supp_model(c(0.1, 5, 7.5)), warn = FALSE))
  expect_output(print(ss), "Sign changes: d1")
})

test_that("the legacy helpers reproduce the published total value (2022)", {
  res <- adjust_analyses(uk_dairy_2022_analyses(3.05, culling = TRUE), method = "published")
  eco <- uk_dairy_2022_economics()
  eco$valuation$culling <- uk_dairy_2022_culling_valuation()
  vls <- lapply(names(res), function(nm) {
    v <- eco$valuation[[nm]]
    legacy_value(legacy_gap(res[[nm]], v$observed, v$direction, v$effect), v$unit_value)
  })
  names(vls) <- names(res)
  vals <- vapply(vls, function(v) v$value, numeric(1))
  expect_equal(unname(sum(vals)) + sum(eco$additional), 402.2475530293503, tolerance = 1e-9)
  expect_equal(eco$additional, c(veterinary = 71.09))
  expect_equal(nrow(vls$yield$by_disease), 13)
  expect_equal(sum(vls$yield$by_disease$value), unname(vals["yield"]), tolerance = 1e-10)
  g <- legacy_gap(res$yield, 8737, "decrease", "percent")
  expect_equal(legacy_value(g, 0.3022)$value, unname(vals["yield"]), tolerance = 1e-10)
})

test_that("plots draw without error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  res <- adjust_analyses(example_uk_dairy_2022(), method = "published")
  expect_invisible(plot(res$yield))
  expect_invisible(plot(res))
  expect_invisible(plot_burden(res))
  expect_invisible(plot_burden(res$yield))
  expect_invisible(plot_burden(deconflate(example_supplement())))
  expect_false("valuation" %in% names(formals(plot_burden)))
  expect_invisible(plot(screen_associations(example_supplement())))
  expect_invisible(plot(sensitivity_oat(example_supplement())))
  s <- cm_sampler(example_supplement(), impacts = list(d1 = dist_normal(2.5, 0.5)))
  expect_invisible(plot(cm_monte_carlo(s, 20, seed = 1)))
  expect_error(plot_burden(example_supplement()), "deconflate")
})
