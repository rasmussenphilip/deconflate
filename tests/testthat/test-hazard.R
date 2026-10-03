# Hazard-ratio adapter (R/hazard.R): cm_hazard_ratios(), cm_hr_model(),
# deconflate_hr() (snapshot, first_order, published) and attributable_risk(),
# plus the internal legacy conversions (R/legacy.R). Reference values from
# inst/validation/reference_v02.py (its "simultaneous" is first_order and its
# "global" is snapshot) and reference_v020_tests.py.

supp_hr_ref <- list(
  published = c(1.4130001998, 1.9090403704, 1.1927014244),
  first_order = c(1.4029251095, 1.8874397472, 1.1691237145),
  snapshot = c(1.3831529713, 1.8909315489, 1.1439870467)
)

# Crude snapshot hazard ratio of each disease over a joint distribution:
# ratio of the mean hazard multiplier among animals with and without it.
snapshot_crude_hr <- function(joint, beta) {
  w <- joint$prob * exp(as.vector(joint$cells %*% beta))
  vapply(seq_len(ncol(joint$cells)), function(i) {
    x <- joint$cells[, i]
    (sum(w[x == 1]) / sum(joint$prob[x == 1])) / (sum(w[x == 0]) / sum(joint$prob[x == 0]))
  }, numeric(1))
}

# Stratified (Mantel-Haenszel-type weighted) log ratio of disease i within
# strata of the diseases in s, as in the snapshot solver.
snapshot_stratified_log_hr <- function(joint, beta, i, s) {
  w <- joint$prob * exp(as.vector(joint$cells %*% beta))
  x <- joint$cells[, i]
  stratum <- as.integer(joint$cells[, s, drop = FALSE] %*% 2^(seq_along(s) - 1))
  num <- 0
  den <- 0
  for (k in unique(stratum)) {
    in_k <- stratum == k
    d1 <- sum(joint$prob[in_k & x == 1])
    d0 <- sum(joint$prob[in_k & x == 0])
    r <- log(sum(w[in_k & x == 1]) / d1) - log(sum(w[in_k & x == 0]) / d0)
    wt <- d1 * d0 / (d1 + d0)
    num <- num + wt * r
    den <- den + wt
  }
  num / den
}

test_that("hazard ratios and hazard-ratio models are validated", {
  expect_error(cm_hazard_ratios("d1", 0, estimand = "snapshot_crude"), "positive")
  expect_error(cm_hazard_ratios(c("d1", "d1"), c(1.5, 2), estimand = "snapshot_crude"), "only one")
  expect_error(cm_hazard_ratios(ids3, c(1.5, 2), estimand = "snapshot_crude"), "one entry per disease")
  expect_error(cm_hazard_ratios("d1", 1.5, adjusted_for = "d2", estimand = "snapshot_crude"), class = "deconflate_unsupported")
  expect_error(cm_hazard_ratios("d1", 1.5, estimand = "snapshot_stratified"), "adjusted_for")
  expect_error(cm_hazard_ratios("d1", 1.5, estimand = "cox"), "estimand")
  h <- cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = "snapshot_crude")
  expect_s3_class(h, "cm_hazard_ratios")
  expect_equal(h$estimand, rep("snapshot_crude", 3))

  expect_error(cm_hr_model(supp_population(), cm_hazard_ratios(c("d1", "d2"), c(1.5, 2), estimand = "snapshot_crude")),
               "No hazard ratio for: d3")
  expect_error(cm_hr_model(supp_population(), cm_hazard_ratios(c(ids3, "zz"), c(1.5, 2, 1.3, 1), estimand = "snapshot_crude")),
               "unknown diseases")
  expect_error(cm_hr_model(supp_population(),
                           cm_hazard_ratios(ids3, c(1.5, 2, 1.3), estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                                            adjusted_for = c("zz", NA, NA))),
               "unknown diseases")
  expect_error(cm_hr_model(supp_population(), data.frame(disease = ids3, value = 1)), "cm_hazard_ratios")

  m <- supp_hr_model()
  expect_s3_class(m, "cm_hr_model")
  expect_s3_class(m$population, "cm_population")
  # Hazard ratios are put in the population's disease order.
  m2 <- cm_hr_model(supp_population(), cm_hazard_ratios(rev(ids3), c(1.3, 2.0, 1.5), estimand = "snapshot_crude"))
  expect_equal(m2$hazard_ratios$disease, ids3)
  expect_equal(m2$hazard_ratios$value, c(1.5, 2.0, 1.3))
  # A cm_model can be given; its impacts are dropped.
  m3 <- cm_hr_model(example_supplement(), cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = "snapshot_crude"))
  expect_null(m3$population$impacts)
  expect_equal(class(m3$population), "cm_population")
  expect_output(print(m), "cm_hr_model")
  expect_error(deconflate_hr(example_supplement()), "cm_hr_model")
})

test_that("the three hazard-ratio methods match the reference (supplement)", {
  m <- supp_hr_model()
  for (meth in names(supp_hr_ref)) {
    r <- deconflate_hr(m, method = meth)
    expect_s3_class(r, "cm_hr_result")
    expect_equal(r$method, meth)
    expect_equal(r$adjusted$disease, ids3)
    expect_equal(r$adjusted$raw, c(1.5, 2.0, 1.3))
    expect_equal(r$adjusted$adjusted, supp_hr_ref[[meth]], tolerance = 1e-7)
    expect_equal(r$adjusted$change, r$adjusted$adjusted / r$adjusted$raw - 1)
    expect_equal(r$diagnostics$n_sign_changes, 0)
  }
  expect_equal(deconflate_hr(m)$method, "snapshot")
  expect_null(deconflate_hr(m, method = "first_order")$joint)
  expect_s3_class(deconflate_hr(m)$joint, "cm_joint")
  expect_output(print(deconflate_hr(m)), "cm_hr_result")
  # The first-order solution solves log(HR_raw) = A beta exactly.
  fo <- deconflate_hr(m, method = "first_order")
  expect_lt(fo$diagnostics$max_reconstruction_residual, 1e-12)
})

test_that("the snapshot method reproduces the raw hazard ratios over the fitted joint", {
  res <- deconflate_hr(supp_hr_model())
  beta <- log(res$adjusted$adjusted)
  expect_equal(snapshot_crude_hr(res$joint, beta), c(1.5, 2.0, 1.3), tolerance = 1e-9)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-9)
  # A supplied joint gives the same answer; a joint of another population
  # (here with a three-way term) is rejected.
  j <- fit_joint(supp_population())
  expect_equal(deconflate_hr(supp_hr_model(), joint = j)$adjusted$adjusted,
               res$adjusted$adjusted, tolerance = 1e-12)
  j3 <- fit_joint(supp_population(cm_three_way("d1", "d2", "d3", 2)))
  expect_error(deconflate_hr(supp_hr_model(), joint = j3), "three-way terms")
})

test_that("three-way terms change the snapshot method only", {
  m3 <- cm_hr_model(supp_population(cm_three_way("d1", "d2", "d3", 2)),
                    cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = "snapshot_crude"))
  expect_equal(deconflate_hr(m3, method = "first_order")$adjusted$adjusted,
               supp_hr_ref$first_order, tolerance = 1e-8)
  expect_equal(deconflate_hr(m3, method = "published")$adjusted$adjusted,
               supp_hr_ref$published, tolerance = 1e-8)
  snap <- deconflate_hr(m3)
  expect_gt(max(abs(snap$adjusted$adjusted - supp_hr_ref$snapshot)), 1e-3)
  # Independent Python solution with the three-way start (reference_v020_tests.py).
  expect_equal(snap$adjusted$adjusted, c(1.3798925799054613, 1.8913281873921548, 1.1399165814382808),
               tolerance = 1e-7)
  # It still reproduces the raw hazard ratios over its own joint.
  expect_equal(snapshot_crude_hr(snap$joint, log(snap$adjusted$adjusted)), c(1.5, 2.0, 1.3),
               tolerance = 1e-9)
})

test_that("adjusted hazard ratios are matched within strata of their adjustment set", {
  hr <- cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                         adjusted_for = c("d2", NA, NA))
  m <- cm_hr_model(supp_population(), hr)
  res <- deconflate_hr(m)
  expect_equal(res$adjusted$estimand, c("snapshot_stratified", "snapshot_crude", "snapshot_crude"))
  expect_equal(res$adjusted$adjusted_for, c("d2", NA, NA))
  beta <- log(res$adjusted$adjusted)
  expect_equal(snapshot_stratified_log_hr(res$joint, beta, 1, 2), log(1.5), tolerance = 1e-9)
  expect_equal(snapshot_crude_hr(res$joint, beta)[2:3], c(2.0, 1.3), tolerance = 1e-9)
  # Independent Python solution of the same equations (reference_v020_tests.py).
  expect_equal(res$adjusted$adjusted, c(1.504714547766569, 1.8751068533634567, 1.1457314393824285),
               tolerance = 1e-7)
  # The published approach is defined for crude hazard ratios only.
  expect_error(deconflate_hr(m, method = "published"), class = "deconflate_unsupported")
  expect_error(deconflate_hr(m, method = "published"), "crude hazard ratios only")
  # The first-order approximation accepts adjusted estimands.
  fo <- deconflate_hr(m, method = "first_order")
  expect_true(all(is.finite(fo$adjusted$adjusted)))
})

test_that("hazard ratios adjusted for all other diseases are used as they are", {
  hr_all <- cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = "snapshot_stratified", adjusted_for = "all")
  m <- cm_hr_model(supp_population(), hr_all)
  for (meth in c("snapshot", "first_order")) {
    r <- deconflate_hr(m, method = meth)
    expect_equal(r$adjusted$adjusted, c(1.5, 2.0, 1.3), tolerance = 1e-10)
    expect_equal(r$adjusted$change, c(0, 0, 0), tolerance = 1e-10)
  }
  expect_error(deconflate_hr(m, method = "published"), class = "deconflate_unsupported")
  # Mixed: d1 adjusted for all, the others crude. d1 keeps its value.
  hr_mix <- cm_hazard_ratios(ids3, c(1.5, 2.0, 1.3), estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                             adjusted_for = c("all", NA, NA))
  mix <- deconflate_hr(cm_hr_model(supp_population(), hr_mix))
  expect_equal(mix$adjusted$adjusted[1], 1.5, tolerance = 1e-10)
  expect_equal(snapshot_crude_hr(mix$joint, log(mix$adjusted$adjusted))[2:3], c(2.0, 1.3),
               tolerance = 1e-9)
})

test_that("hazard ratios of 1 stay 1 and give no attributable risk", {
  m <- cm_hr_model(supp_population(), cm_hazard_ratios(ids3, c(1, 1, 1), estimand = "snapshot_crude"))
  for (meth in c("snapshot", "first_order", "published")) {
    r <- deconflate_hr(m, method = meth)
    expect_equal(r$adjusted$adjusted, c(1, 1, 1), tolerance = 1e-12)
    expect_equal(r$diagnostics$n_sign_changes, 0)
  }
  ar <- attributable_risk(deconflate_hr(m), overall_risk = 0.25)
  expect_equal(ar$summary$disease_free_risk, 0.25, tolerance = 1e-10)
  expect_equal(ar$summary$attributable, 0, tolerance = 1e-10)
  expect_equal(ar$by_disease$attributable, c(0, 0, 0), tolerance = 1e-12)
  expect_true(all(is.na(ar$by_disease$share)))
})

test_that("attributable risk and its Shapley allocation match the reference", {
  res <- deconflate_hr(supp_hr_model())
  ar <- attributable_risk(res, overall_risk = 0.25, unit_value = 1000)
  expect_s3_class(ar, "cm_attributable")
  expect_equal(ar$summary$overall_risk, 0.25)
  expect_equal(ar$summary$disease_free_risk, 0.21335101490465225, tolerance = 1e-8)
  expect_equal(ar$summary$attributable, 0.03664898509534775, tolerance = 1e-7)
  expect_equal(ar$summary$attributable_fraction, 0.03664898509534775 / 0.25, tolerance = 1e-7)
  expect_equal(ar$by_disease$disease, ids3)
  expect_equal(ar$by_disease$hr_adjusted, res$adjusted$adjusted)
  expect_equal(ar$by_disease$attributable, c(0.0073707436, 0.0234877275, 0.0057905139),
               tolerance = 1e-7)
  expect_equal(sum(ar$by_disease$attributable), ar$summary$attributable, tolerance = 1e-10)
  expect_equal(sum(ar$by_disease$share), 1, tolerance = 1e-12)
  expect_equal(ar$summary$unallocated, 0, tolerance = 1e-10)
  expect_equal(ar$summary$value, 1000 * ar$summary$attributable)
  expect_equal(ar$by_disease$value, 1000 * ar$by_disease$attributable)
  expect_equal(ar$baseline_hazard, -log(1 - ar$summary$disease_free_risk), tolerance = 1e-12)
  expect_output(print(ar), "attributable")

  # Without allocation, and with a supplied joint.
  ar0 <- attributable_risk(res, 0.25, allocate = FALSE, joint = fit_joint(supp_population()))
  expect_null(ar0$by_disease)
  expect_true(is.na(ar0$summary$unallocated))
  expect_equal(ar0$summary$attributable, ar$summary$attributable, tolerance = 1e-10)
  # First-order results have no joint; one is fitted.
  ar_fo <- attributable_risk(deconflate_hr(supp_hr_model(), method = "first_order"), 0.25)
  expect_equal(sum(ar_fo$by_disease$attributable), ar_fo$summary$attributable, tolerance = 1e-10)

  expect_error(attributable_risk(res, 1.2), "between 0 and 1")
  expect_error(attributable_risk(res, c(0.2, 0.3)), "between 0 and 1")
  expect_error(attributable_risk(deconflate(example_supplement()), 0.25), "deconflate_hr")
  expect_error(attributable_risk(res, 0.25, joint = fit_joint(supp_population(cm_three_way("d1", "d2", "d3", 2)))),
               "does not match")
})

test_that("attributable risk reports what it cannot allocate", {
  res <- deconflate_hr(supp_hr_model())
  expect_warning(ar <- attributable_risk(res, 0.25, max_present = 2),
                 class = "deconflate_incomplete_allocation")
  s <- ar$summary
  expect_gt(s$unallocated, 0)
  expect_equal(sum(ar$by_disease$attributable) + s$unallocated, s$attributable, tolerance = 1e-12)
  # The unallocated part is the expected excess risk of animals with all three.
  j <- res$joint
  all3 <- rowSums(j$cells) == 3
  excess <- 1 - exp(-ar$baseline_hazard * prod(res$adjusted$adjusted)) - s$disease_free_risk
  expect_equal(s$unallocated, j$prob[all3] * excess, tolerance = 1e-10)
  expect_output(print(ar), "Unallocated")
  # Allocating everything is silent.
  expect_silent(attributable_risk(res, 0.25, max_present = 3))
})

test_that("a joint supplied to the pairwise methods is checked and kept", {
  j <- fit_joint(supp_population())
  fo <- deconflate_hr(supp_hr_model(), method = "first_order", joint = j)
  expect_identical(fo$joint, j)
  expect_equal(fo$adjusted$adjusted, supp_hr_ref$first_order, tolerance = 1e-8)
  expect_equal(attributable_risk(fo, 0.25)$summary$attributable,
               attributable_risk(deconflate_hr(supp_hr_model(), method = "first_order"), 0.25)$summary$attributable,
               tolerance = 1e-12)
  j3 <- fit_joint(supp_population(cm_three_way("d1", "d2", "d3", 2)))
  expect_error(deconflate_hr(supp_hr_model(), method = "published", joint = j3), "three-way terms")
})

test_that("non-positive published hazard ratios are flagged and cannot be valued", {
  # HR 0.5 for d1 with a strong d2: the published denominator of d1 is
  # negative (reference_v02.py hr_methods).
  m <- cm_hr_model(supp_population(), cm_hazard_ratios(ids3, c(0.5, 4, 1), estimand = "snapshot_crude"))
  expect_warning(pub <- deconflate_hr(m, method = "published"), class = "deconflate_nonfinite_warning")
  expect_equal(pub$adjusted$adjusted, c(-0.3585461449, 4.0376405055, 1), tolerance = 1e-8)
  expect_equal(pub$diagnostics$n_sign_changes, 0)
  expect_silent(deconflate_hr(m, method = "published", warn = FALSE))
  expect_error(attributable_risk(pub, 0.25), class = "deconflate_nonfinite")
  # The other methods keep the hazard ratios positive.
  expect_equal(deconflate_hr(m, method = "first_order")$adjusted$adjusted,
               c(0.4268305992, 4.4915295004, 0.7780752535), tolerance = 1e-8)
  snap <- deconflate_hr(m)
  expect_equal(snap$adjusted$adjusted, c(0.4021045208, 4.487174924, 0.7095065051), tolerance = 1e-7)
  expect_equal(attributable_risk(snap, 0.25, allocate = FALSE)$summary$attributable,
               0.03568836025623201, tolerance = 1e-6)
})

test_that("global dairy culling: the three hazard-ratio methods", {
  skip_on_cran()
  m <- example_global_dairy_hr()
  expect_s3_class(m, "cm_hr_model")
  pub <- deconflate_hr(m, method = "published", warn = FALSE)
  expect_equal(pub$adjusted$adjusted,
               c(1.1775803801, 1.9039406263, 2.1979299277, 1.0983766711, 1.3806827011, 1.0124105823,
                 2.6476374499, 1.4586441181, 2.0472349212, 1.2844963465, 1.6752532746, 1.2549275636),
               tolerance = 1e-8)
  # Identical to adjusting HR - 1 additively with the published method, as
  # in Rasmussen et al. (2024).
  legacy <- deconflate(global_dairy_analyses("analysis", culling = TRUE), method = "published",
                       warn = FALSE)$culling_hr_minus_1
  expect_equal(pub$adjusted$adjusted, 1 + legacy$adjusted$adjusted, tolerance = 1e-12)

  expect_warning(deconflate_hr(m, method = "first_order"), class = "deconflate_sign_change")
  fo <- deconflate_hr(m, method = "first_order", warn = FALSE)
  expect_equal(fo$adjusted$adjusted,
               c(0.9619070191, 1.8350863583, 1.825054203, 1.1165198124, 1.2864922835, 0.7369878787,
                 2.6233669144, 1.5487477113, 2.0203081931, 1.2305973733, 1.7412069142, 1.2520841682),
               tolerance = 1e-8)
  expect_equal(fo$diagnostics$sign_changes, "CK, MET")

  snap <- deconflate_hr(m, method = "snapshot", warn = FALSE)
  expect_equal(snap$adjusted$adjusted,
               c(0.9821306495, 1.8193535917, 1.8476664195, 1.050543232, 1.2230209916, 0.7198426768,
                 2.6601485634, 1.5579660463, 2.0165455108, 1.1816994757, 1.7107120034, 1.210900523),
               tolerance = 1e-6)
  expect_equal(snapshot_crude_hr(snap$joint, log(snap$adjusted$adjusted)), m$hazard_ratios$value,
               tolerance = 1e-8)
  expect_equal(snap$diagnostics$n_sign_changes, 2)

  ar <- attributable_risk(snap, overall_risk = 0.2366, allocate = FALSE)
  expect_equal(ar$summary$disease_free_risk, 0.12002780924706569, tolerance = 1e-6)
  expect_equal(ar$summary$attributable, 0.11657219075293432, tolerance = 1e-6)
  expect_null(ar$by_disease)
})

test_that("legacy conversions reproduce the published culling inputs", {
  # Rasmussen et al. (2022), Table 6: hazard ratio treated as an odds ratio
  # with an overall culling rate of 0.27 (DA: 3.83, P = 0.03; LAM: 3.40,
  # P = 0.30).
  ex <- legacy_hr_as_or_excess(c(3.83, 3.40, 1), c(0.03, 0.30, 0.10), 0.27)
  expect_equal(ex, c(0.3138414572431465, 0.2556484948872295, 0), tolerance = 1e-10)
  uk <- uk_dairy_2022_analyses(3.05, culling = TRUE)
  cu <- uk$models$culling$impacts
  expect_equal(cu$value[cu$disease %in% c("DA", "LAM")], ex[1:2], tolerance = 1e-12)
  expect_equal(cu$value[cu$disease == "CO"], 0)
  h <- attr(uk, "hazard_ratios")
  expect_s3_class(h, "cm_hazard_ratios")
  expect_equal(h$value[h$disease == "DA"], 3.83)
  # Eq. 23 rescaling, and the 2024 overall-odds excess.
  expect_equal(legacy_eq23(c(3.83, 1), c(0.31, 0), c(0.2, 0)), c(0.2 * 3.83 / 0.31, 1))
  expect_equal(legacy_overall_odds_excess(c(2.75, 1), 0.27), c(0.23424448217317484, 0),
               tolerance = 1e-12)
})
