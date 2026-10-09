# Rasmussen et al. (2024) global inputs at their central values. Impacts are
# in percent. Reference values computed independently in Python
# (inst/validation/python_reference.py, reference_2024_analysis.py and
# reference_v040.py).

ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")

# The published analyses gave every pair not in Table 3 an odds ratio of 1.
independent_pairs <- function(m) cm_model(with_independent_pairs(as_population(m)), m$impacts)

test_that("incidence is converted with 1 - exp(-I)", {
  m <- example_global_dairy(inputs = "tables")
  d <- m$diseases
  expect_equal(d$prob[d$id == "SCK"], 1 - exp(-0.4789))
  expect_equal(d$prob[d$id == "PTB"], 0.1001)
  expect_equal(d$id, ids)
  expect_equal(nrow(m$associations), 38)
  expect_equal(attr(m$impacts, "units"), "% decrease")
  expect_equal(attr(example_global_dairy("fertility", inputs = "tables")$impacts, "label"),
               "calving interval increase")
})

test_that("published approximation at the means (Tables 2-4) matches the reference", {
  # The published approximation is internal (adjust_impacts()); it is not a
  # method of deconflate().
  y <- adjust_impacts(independent_pairs(example_global_dairy(inputs = "tables")), method = "published")
  f <- adjust_impacts(independent_pairs(example_global_dairy("fertility", inputs = "tables")),
                      method = "published")
  expect_equal(y$adjusted$disease, ids)
  expect_equal(y$adjusted$adjusted,
               c(0.024548, 1.345727, 0.793385, 3.57018, 2.539522, 2.838661,
                 0.070094, 2.639434, 3.227319, 2.259406, 7.10634, 5.559423),
               tolerance = 1e-5)
  expect_equal(f$adjusted$adjusted,
               c(0.328155, 6.206903, 0.152224, 1.064007, 1.118749, 10.755051,
                 1.056262, 7.85979, 3.989669, 3.705326, 0.348395, 0.03033),
               tolerance = 1e-5)
  expect_equal(y$units, "% decrease")
})

test_that("with unlisted pairs independent, the simultaneous solution flags sign changes the published method hides", {
  m <- independent_pairs(example_global_dairy(inputs = "tables"))
  w <- capture_warnings(res <- deconflate(m, n_draws = 0))
  expect_length(w, 1)
  expect_match(w[1], "CK, CM, DA, MF")
  expect_equal(res$method, "simultaneous")
  expect_length(res$notes, 0)
  expect_equal(res$adjusted$adjusted,
               c(-5.363694, -0.249627, -2.672514, 4.313586, 1.987313, 2.774729,
                 -1.520856, 3.218063, 3.926516, 2.459847, 8.274455, 6.572452),
               tolerance = 1e-5)
  expect_equal(res$diagnostics$n_sign_changes, 4)
  expect_equal(res$diagnostics$sign_changes, "CK, CM, DA, MF")
  # The aggregate burden is similar across methods.
  pub <- adjust_impacts(m, method = "published", warn = FALSE)
  expect_equal(pub$totals$adjusted_total, 6.855, tolerance = 1e-3)
  expect_equal(res$totals$adjusted_total, 7.015, tolerance = 1e-3)
})

test_that("the 12-disease joint distribution is feasible and converges", {
  j <- fit_joint(example_global_dairy(inputs = "tables"))
  expect_true(j$converged)
  expect_equal(sum(j$prob), 1, tolerance = 1e-10)
  expect_equal(length(j$prob), 2^12)
  cp <- combination_probs(j, min_prob = 1e-3)
  expect_true(all(diff(cp$prob) <= 0))
})

test_that("analysis inputs use the fixed de-conflation probabilities", {
  gd <- global_dairy_analyses("analysis", culling = TRUE)
  expect_equal(names(gd), c("yield", "fertility", "culling_hr_minus_1"))
  d <- gd$yield$diseases
  expect_equal(d$prob,
               c(0.0301644249, 0.2629368525, 0.0213945, 0.05928092105, 0.2245319578,
                 0.09135551457, 0.02389287642, 0.1065666588, 0.1039640738, 0.1161504785,
                 0.3806298397, 0.4094148319),
               tolerance = 1e-8)
  # SCM was entered without converting incidence to a probability.
  expect_equal(d$type[d$id == "SCM"], "probability")
  # The analyses gave every unlisted pair an odds ratio of 1.
  expect_equal(sum(pair_tables(gd$yield)$status == "unknown"), 0)
  cu <- gd$culling_hr_minus_1$impacts
  expect_equal(cu$value[cu$disease == "CK"], 0.5001)
  hr <- example_global_dairy("culling")
  expect_equal(hr$impacts$value[hr$impacts$disease == "CK"], 1.5001)
  expect_equal(unique(hr$impacts$measure), "HR")
})

test_that("published approximation at the analysis central values matches the reference", {
  gd <- global_dairy_analyses("analysis", culling = TRUE)
  res <- lapply(gd, adjust_impacts, method = "published")
  expect_equal(res$yield$adjusted$adjusted,
               c(0.024732965, 1.3319499, 0.79378291, 3.568819, 2.5324074,
                 2.8409515, 0.068975316, 2.6385934, 3.2291707, 2.2577111,
                 7.1030719, 5.589619),
               tolerance = 1e-6)
  expect_equal(res$fertility$adjusted$adjusted,
               c(0.32648798, 6.1977756, 0.15284518, 1.0618825, 1.1202463,
                 10.760338, 1.0576634, 7.853776, 3.9895865, 3.7076766,
                 0.34933077, 0.032205164),
               tolerance = 1e-6)
  # The 2024 analysis adjusted HR - 1 additively (legacy; reproduction only).
  expect_equal(res$culling_hr_minus_1$adjusted$adjusted,
               c(0.17758038, 0.90394063, 1.1979299, 0.098376671, 0.3806827,
                 0.012410582, 1.6476374, 0.45864412, 1.0472349, 0.28449635,
                 0.67525327, 0.25492756),
               tolerance = 1e-6)
  # The same numbers as the published method of the event model
  # (internal adjust_event()).
  hr <- adjust_event(independent_pairs(example_global_dairy("culling")), method = "published",
                     overall_risk = 0.2366)
  expect_equal(hr$adjusted$adjusted, 1 + res$culling_hr_minus_1$adjusted$adjusted,
               tolerance = 1e-10)
})

test_that("simultaneous solution at the analysis central values", {
  gd <- global_dairy_analyses("analysis", culling = TRUE)
  w <- capture_warnings(res <- lapply(gd, deconflate, n_draws = 0))
  expect_length(w, 3)
  expect_equal(vapply(res, function(r) r$method, character(1)),
               c(yield = "simultaneous", fertility = "simultaneous", culling_hr_minus_1 = "simultaneous"))
  expect_equal(vapply(res, function(r) r$diagnostics$sign_changes, character(1)),
               c(yield = "CK, CM, DA, MF", fertility = "CK, DA, LAM, SCK, SCM",
                 culling_hr_minus_1 = "CK, MET"))
  expect_equal(res$culling_hr_minus_1$adjusted$adjusted,
               c(-0.11369765, 1.0214547, 1.2210513, 0.054545722, 0.27720849,
                 -0.45506511, 1.8664962, 0.58193553, 1.1146889, 0.1870114,
                 0.77426951, 0.19881744),
               tolerance = 1e-6)
})

test_that("by default the 28 unlisted pairs are unknown and the global method fills them in", {
  skip_on_cran()
  m <- example_global_dairy()
  expect_equal(sum(pair_tables(m)$status == "unknown"), 28)
  w <- capture_warnings(y <- deconflate(m, n_draws = 0))
  expect_length(w, 1)
  expect_equal(y$method, "global")
  expect_match(y$notes, "28 pairs without an association")
  expect_equal(y$adjusted$adjusted,
               c(-5.030744, -0.139559, -2.61526, 3.98935, 2.05439, 2.382588, -2.210161,
                 2.680359, 3.121419, 2.144482, 7.886341, 5.949534), tolerance = 1e-5)
  expect_equal(y$totals$raw_sum, 9.9311739814, tolerance = 1e-8)
  expect_equal(y$totals$adjusted_total, 6.9150664186, tolerance = 1e-7)
  expect_equal(nrow(y$unknown_pairs), 28)
  expect_equal(range(y$unknown_pairs$fitted_or), c(1.084282, 2.370312), tolerance = 1e-5)
  expect_equal(y$unknown_pairs$fitted_or[y$unknown_pairs$disease1 == "CK" & y$unknown_pairs$disease2 == "DA"],
               2.370312023777597, tolerance = 1e-6)
  f <- suppressWarnings(deconflate(example_global_dairy("fertility"), n_draws = 0, joint = y$joint))
  expect_equal(f$adjusted$adjusted,
               c(-1.322907, 7.433706, -3.241946, 0.551242, -1.642836, 12.882467, 1.329921,
                 8.504154, 3.913382, 2.721181, -0.2837, -1.647566), tolerance = 1e-5)
  expect_equal(f$totals$adjusted_total, 3.5644537148, tolerance = 1e-7)
})
