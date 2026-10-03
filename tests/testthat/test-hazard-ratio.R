# Hazard-ratio outcomes: published (HR - 1), first-order log-linear and exact
# multiplicative adjustment, and attributable risk. Reference values from
# inst/validation/reference_v02.py.

test_that("hazard-ratio impacts are validated", {
  expect_error(cm_impacts("d1", 0, outcome = "culling", scale = "hazard_ratio"), "positive")
  i <- cm_impacts("d1", 1.5, outcome = "culling", scale = "hazard_ratio")
  expect_equal(i$direction, "increase")
  expect_equal(i$units, "hazard ratio")
})

test_that("the three methods for hazard ratios match the reference", {
  m <- supp_hr()
  cu <- function(method) {
    r <- deconflate(m, method = method)
    r$adjusted$adjusted[r$adjusted$outcome == "culling"]
  }
  expect_equal(cu("published"), c(1.4130001998, 1.9090403704, 1.1927014244), tolerance = 1e-8)
  expect_equal(cu("simultaneous"), c(1.4029251095, 1.8874397472, 1.1691237145), tolerance = 1e-8)
  expect_equal(cu("global"), c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 1e-7)
  # Adding a hazard-ratio outcome leaves the yield adjustment unchanged.
  y <- deconflate(m)$adjusted
  expect_equal(y$adjusted[y$outcome == "yield"],
               deconflate(example_supplement())$adjusted$adjusted)
  expect_equal(adjusted_hr(deconflate(m))$hr, c(1.5, 2.0, 1.3))
})

test_that("the exact multiplicative solution reproduces the raw hazard ratios", {
  res <- deconflate(supp_hr(), method = "global")
  j <- res$joint
  beta <- log(res$adjusted$adjusted[res$adjusted$outcome == "culling"])
  w <- j$prob * exp(as.vector(j$cells %*% beta))
  crude <- vapply(1:3, function(i) {
    (sum(w[j$cells[, i] == 1]) / sum(j$prob[j$cells[, i] == 1])) /
      (sum(w[j$cells[, i] == 0]) / sum(j$prob[j$cells[, i] == 0]))
  }, numeric(1))
  expect_equal(crude, c(1.5, 2.0, 1.3), tolerance = 1e-9)
  d <- res$diagnostics[res$diagnostics$outcome == "culling", ]
  expect_lt(d$max_reconstruction_residual, 1e-9)
})

test_that("attributable risk and its Shapley allocation match the reference", {
  res <- deconflate(supp_hr(), method = "global")
  ar <- attributable_risk(res, overall_risk = 0.25, unit_value = 1000)
  expect_s3_class(ar, "cm_attributable")
  expect_equal(ar$summary$disease_free_risk, 0.21335101490465225, tolerance = 1e-8)
  expect_equal(ar$summary$attributable, 0.03664898509534775, tolerance = 1e-7)
  expect_equal(ar$by_disease$attributable, c(0.0073707436, 0.0234877275, 0.0057905139),
               tolerance = 1e-7)
  expect_equal(sum(ar$by_disease$attributable), ar$summary$attributable, tolerance = 1e-10)
  expect_equal(ar$summary$value, 1000 * ar$summary$attributable)
  expect_output(print(ar), "attributable")
  expect_error(attributable_risk(res, 0.25, outcome = "yield"), "hazard-ratio")
  expect_error(attributable_risk(res, 1.2), "between 0 and 1")
  expect_error(productivity_gap(res, c(culling = 0.25)), "attributable_risk")
  # Hazard-ratio outcomes are left out of the additive burden.
  expect_false("culling" %in% attribute_burden(res)$outcome)
})

test_that("summary values hazard-ratio outcomes through attributable risk", {
  res <- deconflate(supp_hr(), method = "global")
  eco <- list(observed = c(yield = 10000, culling = 0.25),
              unit_value = c(yield = 0.3, culling = 1000))
  s <- summary(res, economics = eco)
  cu <- s$totals[s$totals$outcome == "culling", ]
  expect_equal(cu$gap, 0.03664898509534775, tolerance = 1e-7)
  expect_equal(cu$value, 1000 * cu$gap)
  ct <- contribution_table(res, eco)
  expect_equal(sum(ct$value[ct$outcome == "culling"]), cu$value, tolerance = 1e-8)
})

test_that("global dairy culling: the three hazard-ratio methods", {
  skip_on_cran()
  m <- example_global_dairy(culling_scale = "hazard_ratio")
  cu <- function(r) r$adjusted$adjusted[r$adjusted$outcome == "culling"]
  pub <- deconflate(m, method = "published", warn = FALSE)
  # Identical to adjusting HR - 1 as in Rasmussen et al. (2024).
  expect_equal(cu(pub), 1 + cu(deconflate(example_global_dairy(), method = "published")))
  expect_equal(cu(pub), c(1.1775803801, 1.9039406263, 2.1979299277, 1.0983766711,
                          1.3806827011, 1.0124105823, 2.6476374499, 1.4586441181,
                          2.0472349212, 1.2844963465, 1.6752532746, 1.2549275636),
               tolerance = 1e-8)
  sim <- deconflate(m, method = "simultaneous", warn = FALSE)
  expect_equal(cu(sim), c(0.9619070191, 1.8350863583, 1.825054203, 1.1165198124,
                          1.2864922835, 0.7369878787, 2.6233669144, 1.5487477113,
                          2.0203081931, 1.2305973733, 1.7412069142, 1.2520841682),
               tolerance = 1e-8)
  glob <- deconflate(m, method = "global", warn = FALSE)
  expect_equal(cu(glob), c(0.9821306495, 1.8193535917, 1.8476664195, 1.050543232,
                           1.2230209916, 0.7198426768, 2.6601485634, 1.5579660463,
                           2.0165455108, 1.1816994757, 1.7107120034, 1.210900523),
               tolerance = 1e-6)
  ar <- attributable_risk(glob, overall_risk = 0.2366, allocate = FALSE)
  expect_equal(ar$summary$attributable, 0.11657219075293432, tolerance = 1e-6)
  expect_null(ar$by_disease)
})

test_that("hazard-ratio draws must be positive", {
  s <- cm_sampler(supp_hr(), impacts = list("culling:d1" = dist_normal(0.2, 1)))
  mc <- cm_monte_carlo(s, 50, seed = 1)
  expect_gt(mc$n_rejected, 0)
  expect_match(mc$rejections$reason[1], "not a valid hazard ratio")
})
