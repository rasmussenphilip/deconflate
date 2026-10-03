# Rasmussen et al. (2024) global inputs at their central values.
# Reference values computed independently in Python
# (inst/validation/python_reference.py and reference_2024_analysis.py).
# The first tests use Tables 2-4 as printed (inputs = "tables").

ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")

test_that("incidence is converted with 1 - exp(-I)", {
  m <- example_global_dairy(inputs = "tables", culling = FALSE)
  expect_equal(m$diseases$prob[m$diseases$id == "SCK"], 1 - exp(-0.4789))
  expect_equal(m$diseases$prob[m$diseases$id == "PTB"], 0.1001)
})

test_that("published method at the means matches the reference (close to Table 5)", {
  res <- deconflate(example_global_dairy(inputs = "tables", culling = FALSE), method = "published")
  y <- res$adjusted[res$adjusted$outcome == "yield", ]
  expect_equal(y$disease, ids)
  expect_equal(y$adjusted * 100,
               c(0.024548, 1.345727, 0.793385, 3.57018, 2.539522, 2.838661,
                 0.070094, 2.639434, 3.227319, 2.259406, 7.10634, 5.559423),
               tolerance = 1e-5)
  f <- res$adjusted[res$adjusted$outcome == "fertility", ]
  expect_equal(f$adjusted * 100,
               c(0.328155, 6.206903, 0.152224, 1.064007, 1.118749, 10.755051,
                 1.056262, 7.85979, 3.989669, 3.705326, 0.348395, 0.03033),
               tolerance = 1e-5)
})

test_that("simultaneous solution flags sign changes the published method hides", {
  m <- example_global_dairy(inputs = "tables", culling = FALSE)
  # One sign-change warning per outcome (yield and fertility).
  w <- capture_warnings(res <- deconflate(m))
  expect_length(w, 2)
  expect_match(w[1], "CK, CM, DA, MF")
  y <- res$adjusted[res$adjusted$outcome == "yield", ]
  expect_equal(y$adjusted * 100,
               c(-5.363694, -0.249627, -2.672514, 4.313586, 1.987313, 2.774729,
                 -1.520856, 3.218063, 3.926516, 2.459847, 8.274455, 6.572452),
               tolerance = 1e-5)
  d <- res$diagnostics[res$diagnostics$outcome == "yield", ]
  expect_equal(d$n_sign_changes, 4)
  expect_equal(d$sign_changes, "CK, CM, DA, MF")
  # Aggregate burden is similar across methods.
  P <- m$diseases$prob
  pub <- deconflate(m, method = "published", warn = FALSE)$adjusted
  expect_equal(sum(pub$adjusted[pub$outcome == "yield"] * P), 0.06855, tolerance = 1e-3)
  expect_equal(sum(y$adjusted * P), 0.07015, tolerance = 1e-3)
})

test_that("the 12-disease joint distribution is feasible and converges", {
  j <- fit_joint(example_global_dairy(inputs = "tables", culling = FALSE))
  expect_true(j$converged)
  expect_equal(sum(j$prob), 1, tolerance = 1e-10)
  expect_equal(length(j$prob), 2^12)
  cp <- combination_probs(j, min_prob = 1e-3)
  expect_true(all(diff(cp$prob) <= 0))
})

test_that("analysis inputs use the fixed de-conflation probabilities", {
  m <- example_global_dairy()
  expect_equal(unique(m$impacts$outcome), c("yield", "fertility", "culling"))
  expect_equal(m$diseases$prob,
               c(0.0301644249, 0.2629368525, 0.0213945, 0.05928092105, 0.2245319578,
                 0.09135551457, 0.02389287642, 0.1065666588, 0.1039640738, 0.1161504785,
                 0.3806298397, 0.4094148319),
               tolerance = 1e-8)
  # SCM was entered without converting incidence to a probability.
  expect_equal(m$diseases$type[m$diseases$id == "SCM"], "probability")
  cull <- m$impacts[m$impacts$outcome == "culling", ]
  expect_equal(cull$value[cull$disease == "CK"], 0.5001)
  expect_equal(cull$scale[1], "absolute")
})

test_that("published method at the analysis central values matches the reference", {
  res <- deconflate(example_global_dairy(), method = "published")
  a <- res$adjusted
  expect_equal(a$adjusted[a$outcome == "yield"],
               c(0.00024732965, 0.013319499, 0.0079378291, 0.03568819, 0.025324074,
                 0.028409515, 0.00068975316, 0.026385934, 0.032291707, 0.022577111,
                 0.071030719, 0.05589619),
               tolerance = 1e-6)
  expect_equal(a$adjusted[a$outcome == "fertility"],
               c(0.0032648798, 0.061977756, 0.0015284518, 0.010618825, 0.011202463,
                 0.10760338, 0.010576634, 0.07853776, 0.039895865, 0.037076766,
                 0.0034933077, 0.00032205164),
               tolerance = 1e-6)
  expect_equal(a$adjusted[a$outcome == "culling"],
               c(0.17758038, 0.90394063, 1.1979299, 0.098376671, 0.3806827,
                 0.012410582, 1.6476374, 0.45864412, 1.0472349, 0.28449635,
                 0.67525327, 0.25492756),
               tolerance = 1e-6)
  hr <- adjusted_hr(res, method = "excess_hr")
  expect_equal(hr$hr_adjusted, 1 + a$adjusted[a$outcome == "culling"])
  expect_equal(hr$hr[hr$disease == "CK"], 1.5001)
})

test_that("simultaneous solution at the analysis central values", {
  w <- capture_warnings(res <- deconflate(example_global_dairy()))
  expect_length(w, 3)
  d <- res$diagnostics
  expect_equal(d$sign_changes, c("CK, CM, DA, MF", "CK, DA, LAM, SCK, SCM", "CK, MET"))
  cu <- res$adjusted[res$adjusted$outcome == "culling", ]
  expect_equal(cu$adjusted,
               c(-0.11369765, 1.0214547, 1.2210513, 0.054545722, 0.27720849,
                 -0.45506511, 1.8664962, 0.58193553, 1.1146889, 0.1870114,
                 0.77426951, 0.19881744),
               tolerance = 1e-6)
})
