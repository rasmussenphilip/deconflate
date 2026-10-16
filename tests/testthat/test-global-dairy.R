# Rasmussen et al. (2024) global inputs at their central values. Impacts are
# in percent. Reference values computed independently in Python
# (inst/validation/python_reference.py, reference_2024_analysis.py and
# reference_v040.py).

ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")

# Every pair not in Table 3 given an odds ratio of 1, as in the 2024 analyses.
independent_pairs <- function(m) {
  pop <- as_population(m)
  pt <- pair_tables(pop)
  u <- pt[pt$status == "unknown", , drop = FALSE]
  add <- cm_associations(u$disease1, u$disease2, 1, measure = "OR")
  add$row <- max(pop$associations$row) + seq_len(nrow(add))
  pop$associations <- rbind(pop$associations, add)
  cm_model(pop, m$impacts)
}

test_that("incidence is converted with 1 - exp(-I)", {
  m <- example_global_dairy(inputs = "tables")
  d <- m$diseases
  expect_equal(d$prob[d$disease == "SCK"], 1 - exp(-0.4789))
  expect_equal(d$prob[d$disease == "PTB"], 0.1001)
  expect_equal(d$disease, ids)
  expect_equal(nrow(m$associations), 38)
  expect_equal(attr(m$impacts, "units"), "% decrease")
  expect_equal(attr(example_global_dairy("fertility", inputs = "tables")$impacts, "label"),
               "calving interval increase")
})

test_that("with unlisted pairs independent, the model flags sign changes", {
  m <- independent_pairs(example_global_dairy(inputs = "tables"))
  w <- capture_warnings(res <- deconflate(m, n_draws = 0))
  expect_length(w, 1)
  expect_match(w[1], "CK, CM, DA, MF")
  expect_equal(res$method, "pairwise shortcut")
  expect_equal(grep("n_draws = 0", res$notes, invert = TRUE, value = TRUE), character(0))
  expect_equal(res$adjusted$adjusted,
               c(-5.363694, -0.249627, -2.672514, 4.313586, 1.987313, 2.774729,
                 -1.520856, 3.218063, 3.926516, 2.459847, 8.274455, 6.572452),
               tolerance = 1e-5)
  expect_equal(res$diagnostics$n_sign_changes, 4)
  expect_equal(res$diagnostics$sign_changes, "CK, CM, DA, MF")
  expect_equal(res$totals$adjusted_total, 7.015, tolerance = 1e-3)
})

test_that("the 12-disease joint distribution is feasible and converges", {
  j <- fit_joint(example_global_dairy(inputs = "tables"))
  expect_true(j$converged)
  expect_equal(sum(j$prob), 1, tolerance = 1e-10)
  expect_equal(length(j$prob), 2^12)
})

test_that("the culling example has hazard-ratio impacts", {
  hr <- example_global_dairy("culling")
  expect_equal(hr$impacts$value[hr$impacts$disease == "CK"], 1.5001)
  expect_equal(unique(hr$impacts$measure), "HR")
})

test_that("by default the 28 unlisted pairs are unknown and the model fills them in", {
  skip_on_cran()
  m <- example_global_dairy()
  expect_equal(sum(pair_tables(m)$status == "unknown"), 28)
  w <- capture_warnings(y <- deconflate(m, n_draws = 0))
  expect_length(w, 1)
  expect_equal(y$method, "full model")
  expect_equal(grep("n_draws = 0", y$notes, invert = TRUE, value = TRUE), character(0))
  expect_equal(y$adjusted$adjusted,
               c(-5.030744, -0.139559, -2.61526, 3.98935, 2.05439, 2.382588, -2.210161,
                 2.680359, 3.121419, 2.144482, 7.886341, 5.949534), tolerance = 1e-5)
  expect_equal(y$totals$raw_total, 9.9311739814, tolerance = 1e-8)
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
