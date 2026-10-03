# Rasmussen et al. (2024) global inputs at their means.
# Reference values computed independently in Python.

ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")

test_that("incidence is converted with 1 - exp(-I)", {
  m <- example_global_dairy()
  expect_equal(m$diseases$prob[m$diseases$id == "SCK"], 1 - exp(-0.4789))
  expect_equal(m$diseases$prob[m$diseases$id == "PTB"], 0.1001)
})

test_that("published method at the means matches the reference (close to Table 5)", {
  res <- deconflate(example_global_dairy(), method = "published")
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
  m <- example_global_dairy()
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
  j <- fit_joint(example_global_dairy())
  expect_true(j$converged)
  expect_equal(sum(j$prob), 1, tolerance = 1e-10)
  expect_equal(length(j$prob), 2^12)
  cp <- combination_probs(j, min_prob = 1e-3)
  expect_true(all(diff(cp$prob) <= 0))
})
