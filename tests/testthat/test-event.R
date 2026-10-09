# Event impacts (R/hazard.R): cm_impacts(measure = ...) and
# deconflate(event_model = TRUE, overall_risk = ...) with the snapshot hazard
# model; its first-order and published (Rasmussen et al. 2024)
# approximations (internal adjust_event(), compare_methods()); the
# attributable risk and its Shapley allocation; and the internal legacy
# conversions (R/legacy.R). Reference values: inst/validation/reference_v02.py
# (its "simultaneous" is first_order and its "global" is snapshot),
# reference_v020_tests.py and reference_v040.py (sections 2, 4, 5, 8, 9).

t2_ids <- c("d1", "d2", "d3")

# Supplementary File population (Rasmussen et al. 2022): P = 0.10, 0.15,
# 0.20; OR d1:d2 = 2, d1:d3 = 1, d2:d3 = 3.
t2_pop <- function(three_way = NULL) {
  cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)),
                three_way = three_way)
}

# Event model on a population (default: culling hazard ratios 1.5, 2.0, 1.3).
t2_event <- function(values = c(1.5, 2.0, 1.3), measure = "HR", estimand = "snapshot_crude",
                     adjusted_for = NA_character_, pop = t2_pop()) {
  cm_model(pop, cm_impacts(t2_ids, values, measure = measure, estimand = estimand,
                           adjusted_for = adjusted_for, label = "culling"))
}

t2_hr_ref <- list(
  published = c(1.4130001998, 1.9090403704, 1.1927014244),
  first_order = c(1.4029251095, 1.8874397472, 1.1691237145),
  snapshot = c(1.3831529713, 1.8909315489, 1.1439870467)
)
t2_ar_ref <- 0.03664898509534775   # attributable risk of the snapshot solution at 0.25

# The snapshot measure of disease i over a joint distribution, on the
# measure's own scale (log for ratios, the difference for RD): crude, or
# within strata of the diseases in columns s combined with weights
# d1 d0 / (d1 + d0). h0: baseline hazard (risk-based measures).
t2_snapshot_measure <- function(joint, beta, i, s = integer(0), kind = "HR", h0 = NULL) {
  mult <- exp(as.vector(joint$cells %*% beta))
  val <- if (kind %in% c("HR", "rate_ratio")) mult else 1 - exp(-h0 * mult)
  x <- joint$cells[, i]
  stratum <- if (length(s)) as.integer(joint$cells[, s, drop = FALSE] %*% 2^(seq_along(s) - 1)) else
    rep(0L, length(x))
  num <- 0
  den <- 0
  for (k in unique(stratum)) {
    in_k <- stratum == k
    d1 <- sum(joint$prob[in_k & x == 1])
    d0 <- sum(joint$prob[in_k & x == 0])
    if (d1 <= 0 || d0 <= 0) next
    m1 <- sum((joint$prob * val)[in_k & x == 1]) / d1
    m0 <- sum((joint$prob * val)[in_k & x == 0]) / d0
    e <- switch(kind, HR = , rate_ratio = , RR = log(m1 / m0),
                OR = log(m1 / (1 - m1)) - log(m0 / (1 - m0)), RD = m1 - m0)
    wt <- d1 * d0 / (d1 + d0)
    num <- num + wt * e
    den <- den + wt
  }
  num / den
}

# Five-disease example (inst/extdata/five_diseases, point values): culling
# with mixed measures, crude and stratified estimands, and the LAM:MAS:SCK
# three-way ratio 1.5.
t2_five_culling <- function() {
  ids <- c("LAM", "MAS", "MET", "SCK", "RP")
  pop <- cm_population(
    cm_diseases(ids, c(0.25, 0.30, 0.10, 0.35, 0.06),
                type = c("prevalence", "incidence_rate", "prevalence", "probability", "prevalence")),
    cm_associations(c("LAM", "LAM", "MET", "MET", "MAS", "MAS", "LAM", "RP", "LAM", "MAS"),
                    c("MAS", "SCK", "RP", "SCK", "SCK", "MET", "MET", "SCK", "RP", "RP"),
                    c(1.8, 2.0, 6.0, 1.5, 0.35, 0.08, 0.05, 1.4296875, 1, 1.5),
                    measure = c("OR", "OR", "OR", "RR", "cond_prob", "RD", "phi", "OR", "OR", "OR")),
    three_way = cm_three_way("LAM", "MAS", "SCK", 1.5))
  cm_model(pop, cm_impacts(ids, c(1.74, 1.6, 1.23, 1.45, 0.10),
                           measure = c("HR", "HR", "RR", "OR", "RD"),
                           estimand = c("snapshot_crude", "snapshot_stratified", "snapshot_crude",
                                        "snapshot_crude", "snapshot_stratified"),
                           adjusted_for = c(NA, "LAM; SCK", NA, NA, "all"), label = "culling"))
}

test_that("event impacts are validated by cm_impacts()", {
  ev <- function(...) cm_impacts(..., measure = "HR", estimand = "snapshot_crude")
  expect_error(ev("d1", 0), "positive")
  expect_error(ev(c("d1", "d1"), c(1.5, 2)), "only one")
  expect_error(ev(t2_ids, c(1.5, 2)), "one entry per disease")
  expect_error(cm_impacts("d1", 1.5, measure = "HR", estimand = "snapshot_crude", adjusted_for = "d2"),
               class = "deconflate_unsupported")
  expect_error(cm_impacts("d1", 1.5, measure = "HR", estimand = "snapshot_stratified"), "adjusted_for")
  expect_error(cm_impacts("d1", 1.5, measure = "HR", estimand = "cox"), "estimand")
  expect_error(cm_impacts("d1", 1.5, measure = "HR", estimand = "crude"), "estimand")
  # The estimand of event impacts has no default.
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR"), "State the estimand")
  expect_error(cm_impacts("d1", 1.5, measure = "hazard", estimand = "snapshot_crude"), "measure")
  # Ratios are positive; risk differences lie in (-1, 1).
  for (m in c("HR", "rate_ratio", "RR", "OR")) {
    expect_error(cm_impacts("d1", -0.5, measure = m, estimand = "snapshot_crude"), "positive")
  }
  expect_error(cm_impacts("d1", 1, measure = "RD", estimand = "snapshot_crude"), "between -1 and 1")
  expect_error(cm_impacts("d1", -1, measure = "RD", estimand = "snapshot_crude"), "between -1 and 1")
  expect_silent(cm_impacts(c("d1", "d2"), c(-0.2, 0), measure = "RD", estimand = "snapshot_crude"))
  # Snapshot estimands belong to event impacts.
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), estimand = "snapshot_crude"),
               class = "deconflate_unsupported")
  expect_error(cm_impacts(t2_ids, c(1.5, 2, 1.3), estimand = "snapshot_stratified", adjusted_for = "all"),
               class = "deconflate_unsupported")

  h <- cm_impacts(t2_ids, c(1.5, 0.04, 1.3), measure = c("HR", "RD", "OR"), estimand = "snapshot_crude")
  expect_s3_class(h, "cm_impacts")
  expect_equal(attr(h, "kind"), "event")
  expect_equal(h$measure, c("HR", "RD", "OR"))
  expect_equal(h$estimand, rep("snapshot_crude", 3))
  expect_equal(impact_kind(h), "event")
  expect_equal(attr(cm_impacts(t2_ids, c(1, 2, 3)), "kind"), "additive")
  expect_true(all(is.na(cm_impacts(t2_ids, c(1, 2, 3))$measure)))
})

test_that("event models are validated by cm_model()", {
  ev <- function(ids, v, ...) cm_impacts(ids, v, measure = "HR", estimand = "snapshot_crude", ...)
  expect_error(cm_model(t2_pop(), ev(c("d1", "d2"), c(1.5, 2))), "No impact for: d3")
  expect_error(cm_model(t2_pop(), ev(c("d1", "d2"), c(1.5, 2))), "1 for ratios, 0 for risk differences")
  expect_error(cm_model(t2_pop(), ev(c(t2_ids, "zz"), c(1.5, 2, 1.3, 1))), "unknown diseases")
  expect_error(cm_model(t2_pop(), cm_impacts(t2_ids, c(1.5, 2, 1.3), measure = "HR",
                                             estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                                             adjusted_for = c("zz", NA, NA))),
               "unknown diseases")
  # Impacts are put in the population's disease order.
  m <- cm_model(t2_pop(), cm_impacts(rev(t2_ids), c(0.05, 2.0, 1.5), measure = c("RD", "HR", "OR"),
                                     estimand = "snapshot_crude"))
  expect_equal(m$impacts$disease, t2_ids)
  expect_equal(m$impacts$value, c(1.5, 2.0, 0.05))
  expect_equal(m$impacts$measure, c("OR", "HR", "RD"))
  expect_equal(impact_kind(m$impacts), "event")
  # Interactions apply to additive impacts only.
  expect_error(cm_model(t2_pop(), ev(t2_ids, c(1.5, 2, 1.3)), cm_interactions("d1", "d2", 0.5)),
               class = "deconflate_unsupported")
  expect_error(cm_model(t2_pop(), ev(t2_ids, c(1.5, 2, 1.3)), cm_interactions("d1", "d2", 0.5)),
               "additive impacts only")
  expect_error(set_interaction(t2_event(), "d1", "d2", 0.5), class = "deconflate_unsupported")
  expect_output(print(t2_event()), "event impacts: use event_model = TRUE")
  expect_output(print(t2_event()), "Measures: HR: 3")
})

test_that("additive and event impacts must match event_model, which needs overall_risk", {
  m <- t2_event()
  expect_error(deconflate(m), class = "deconflate_unsupported")
  expect_error(deconflate(m), "event_model = TRUE")
  expect_error(deconflate(example_supplement(), event_model = TRUE, overall_risk = 0.25),
               class = "deconflate_unsupported")
  expect_error(deconflate(example_supplement(), overall_risk = 0.25), class = "deconflate_unsupported")
  expect_error(deconflate(example_supplement(), overall_risk = 0.25), "only with event_model")
  expect_error(deconflate(m, event_model = NA, overall_risk = 0.25), "TRUE or FALSE")
  # overall_risk is required, also when every row is a hazard ratio.
  expect_error(deconflate(m, event_model = TRUE), "overall_risk")
  for (bad in list(1.2, 0, 1, c(0.2, 0.3), NA_real_, "0.25")) {
    expect_error(deconflate(m, event_model = TRUE, overall_risk = bad), "between 0 and 1")
  }
  expect_error(deconflate(m, event_model = TRUE, overall_risk = dist_normal(0.25, 0.1)),
               "must lie between 0 and 1")
  # The internal engines check the kind of impacts too.
  expect_error(adjust_impacts(m), class = "deconflate_unsupported")
  expect_error(adjust_event(example_supplement(), overall_risk = 0.25), class = "deconflate_unsupported")
  expect_error(compare_methods(m), class = "deconflate_unsupported")
  expect_error(screen_associations(m), class = "deconflate_unsupported")
  expect_error(screen_associations(example_supplement(), event_model = TRUE, overall_risk = 0.25),
               class = "deconflate_unsupported")
})

test_that("the snapshot model matches the reference (supplement hazard ratios)", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25)
  expect_s3_class(res, "cm_event_result")
  expect_equal(res$method, "snapshot")
  expect_equal(res$overall_risk, 0.25)
  expect_equal(res$adjusted$disease, t2_ids)
  expect_equal(res$adjusted$measure, rep("HR", 3))
  expect_equal(res$adjusted$raw, c(1.5, 2.0, 1.3))
  expect_equal(res$adjusted$adjusted, t2_hr_ref$snapshot, tolerance = 1e-7)
  expect_equal(res$adjusted$change, res$adjusted$adjusted / res$adjusted$raw - 1)
  expect_equal(res$adjusted$estimand, rep("snapshot_crude", 3))
  expect_equal(res$diagnostics$n_sign_changes, 0)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-9)
  expect_s3_class(res$joint, "cm_joint")
  expect_equal(nrow(res$unknown_pairs), 0)
  expect_null(res$draws)
  expect_equal(res$label, "culling")
  # No input has a distribution: no draws, and a note says so.
  expect_length(res$notes, 1)
  expect_match(res$notes, "No input has a distribution")
  expect_length(deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)$notes, 0)
  expect_output(print(res), "cm_event_result")
  expect_output(print(res), "Attributable risk by disease")
  # The method argument is not used for event impacts.
  expect_equal(deconflate(t2_event(), method = "global", event_model = TRUE, overall_risk = 0.25)$method,
               "snapshot")
  rs <- deconflate(t2_event(), method = "simultaneous", event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(rs$adjusted$adjusted, res$adjusted$adjusted)
  expect_match(rs$notes, "was not used", all = FALSE)
  # The published approximation is not a method of deconflate().
  expect_error(deconflate(t2_event(), method = "published", event_model = TRUE, overall_risk = 0.25),
               class = "deconflate_unsupported")
  expect_error(deconflate(t2_event(), method = "published", event_model = TRUE, overall_risk = 0.25),
               "compare_methods")
})

test_that("first-order and published approximations match the reference", {
  m <- t2_event()
  for (meth in c("published", "first_order", "snapshot")) {
    r <- adjust_event(m, method = meth, overall_risk = 0.25)
    expect_s3_class(r, "cm_event_result")
    expect_equal(r$method, meth)
    expect_equal(r$adjusted$adjusted, t2_hr_ref[[meth]], tolerance = 1e-7)
    expect_equal(r$diagnostics$n_sign_changes, 0)
    # Every method has a joint distribution (for the attributable risk).
    expect_s3_class(r$joint, "cm_joint")
    expect_equal(sum(r$attributable$by_disease$attributable), r$attributable$summary$attributable,
                 tolerance = 1e-10)
  }
  # The first-order solution solves log(HR_raw) = A beta exactly.
  fo <- adjust_event(m, method = "first_order", overall_risk = 0.25)
  expect_lt(fo$diagnostics$max_reconstruction_residual, 1e-12)

  cmp <- compare_methods(m, event_model = TRUE, overall_risk = 0.25)
  expect_s3_class(cmp, "cm_comparison")
  expect_equal(cmp$methods, c("published", "first_order", "snapshot"))
  expect_equal(names(cmp$impacts), c("disease", "measure", "raw", "published", "first_order", "snapshot"))
  for (meth in names(t2_hr_ref)) {
    expect_equal(cmp$impacts[[meth]], t2_hr_ref[[meth]], tolerance = 1e-7)
  }
  expect_equal(names(cmp$totals), c("method", "overall_risk", "disease_free_risk", "attributable"))
  expect_equal(cmp$totals$attributable[cmp$totals$method == "snapshot"], t2_ar_ref, tolerance = 1e-7)
  expect_equal(cmp$totals$overall_risk, rep(0.25, 3))
  expect_false(any(cmp$long$sign_change))
  expect_length(cmp$failed, 0)
  expect_output(print(cmp), "Adjusted hazard ratios")
  expect_error(compare_methods(m, event_model = TRUE), "overall_risk")
})

test_that("the snapshot model reproduces the raw hazard ratios over the fitted joint", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  beta <- log(res$adjusted$adjusted)
  expect_equal(exp(vapply(1:3, function(i) t2_snapshot_measure(res$joint, beta, i), numeric(1))),
               c(1.5, 2.0, 1.3), tolerance = 1e-9)
  # A supplied joint gives the same answer; a joint of another population
  # (here with a three-way term) is rejected.
  j <- fit_joint(t2_pop())
  expect_equal(deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, joint = j)$adjusted$adjusted,
               res$adjusted$adjusted, tolerance = 1e-12)
  j3 <- fit_joint(t2_pop(cm_three_way("d1", "d2", "d3", 2)))
  expect_error(deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, joint = j3), "three-way terms")
  # A joint supplied to the pairwise approximations is checked and kept.
  fo <- adjust_event(t2_event(), method = "first_order", overall_risk = 0.25, joint = j)
  expect_identical(fo$joint, j)
  expect_equal(fo$adjusted$adjusted, t2_hr_ref$first_order, tolerance = 1e-8)
  expect_error(adjust_event(t2_event(), method = "published", overall_risk = 0.25, joint = j3), "three-way terms")
})

test_that("hazard-ratio results do not depend on the overall risk; the attributable risk does", {
  r1 <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.1, n_draws = 0)
  r4 <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.4, n_draws = 0)
  expect_equal(r1$adjusted$adjusted, t2_hr_ref$snapshot, tolerance = 1e-7)
  expect_equal(r4$adjusted$adjusted, r1$adjusted$adjusted, tolerance = 1e-9)
  expect_gt(r4$attributable$summary$attributable, r1$attributable$summary$attributable + 0.01)
  # The same for rate ratios, which are treated like hazard ratios.
  rr <- deconflate(t2_event(measure = c("rate_ratio", "HR", "rate_ratio")), event_model = TRUE,
                   overall_risk = 0.25, n_draws = 0)
  expect_equal(rr$adjusted$measure, c("rate_ratio", "HR", "rate_ratio"))
  expect_equal(rr$adjusted$adjusted, t2_hr_ref$snapshot, tolerance = 1e-7)
  expect_equal(rr$adjusted$change, rr$adjusted$adjusted / rr$adjusted$raw - 1)
  expect_equal(rr$attributable$summary$attributable, t2_ar_ref, tolerance = 1e-7)
})

test_that("risk ratios, odds ratios and risk differences reproduce the hazard-ratio solution", {
  # The crude OR of d1, RR of d2 and RD of d3 implied by the hazard-ratio
  # solution at an overall risk of 0.25 (reference_v040.py, section 5).
  exact <- t2_event(c(1.592048950260, 1.770988234178, 0.057214512738), measure = c("OR", "RR", "RD"))
  re <- deconflate(exact, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(re$adjusted$adjusted, t2_hr_ref$snapshot, tolerance = 1e-8)
  expect_equal(re$attributable$summary$attributable, t2_ar_ref, tolerance = 1e-8)
  # The change column is defined for hazard and rate ratios only.
  expect_true(all(is.na(re$adjusted$change)))
  expect_equal(re$adjusted$measure, c("OR", "RR", "RD"))
  # Rounded inputs (reference_v040.py, section 5).
  mixed <- t2_event(c(1.592, 1.771, 0.0572), measure = c("OR", "RR", "RD"))
  rm <- deconflate(mixed, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(rm$adjusted$adjusted, c(1.383112894, 1.8909782708, 1.1439101733), tolerance = 1e-8)
  # The model reproduces each raw estimate on its own scale.
  h0 <- rm$attributable$baseline_hazard
  beta <- log(rm$adjusted$adjusted)
  expect_equal(c(exp(t2_snapshot_measure(rm$joint, beta, 1, kind = "OR", h0 = h0)),
                 exp(t2_snapshot_measure(rm$joint, beta, 2, kind = "RR", h0 = h0)),
                 t2_snapshot_measure(rm$joint, beta, 3, kind = "RD", h0 = h0)),
               c(1.592, 1.771, 0.0572), tolerance = 1e-9)
  # Risk-based rows depend on the overall risk.
  r15 <- deconflate(mixed, event_model = TRUE, overall_risk = 0.15, n_draws = 0)
  expect_gt(max(abs(r15$adjusted$adjusted - rm$adjusted$adjusted)), 1e-3)
  expect_error(deconflate(mixed, event_model = TRUE), "overall_risk")
  # Only the snapshot model handles risk-based rows.
  for (meth in c("first_order", "published")) {
    expect_error(adjust_event(mixed, method = meth, overall_risk = 0.25), class = "deconflate_unsupported")
  }
  cmp <- compare_methods(mixed, event_model = TRUE, overall_risk = 0.25)
  expect_equal(cmp$methods, "snapshot")
  expect_setequal(names(cmp$failed), c("published", "first_order"))
  expect_match(cmp$failed[["first_order"]], "hazard or rate ratios")
})

test_that("mixed measures, crude and stratified, match the five-disease reference", {
  m <- t2_five_culling()
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(res$adjusted$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(res$adjusted$measure, c("HR", "HR", "RR", "OR", "RD"))
  expect_equal(res$adjusted$adjusted, c(1.58344491, 1.57873827, 1.08072624, 1.19114596, 1.50176712),
               tolerance = 1e-7)
  expect_equal(is.na(res$adjusted$change), c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(res$adjusted$adjusted_for, c(NA, "LAM; SCK", NA, NA, "all"))
  a <- res$attributable
  expect_equal(a$summary$attributable, 0.0697502441, tolerance = 1e-8)
  expect_equal(a$by_disease$attributable,
               c(0.0248718025, 0.0256191332, 0.0015698203, 0.0123310667, 0.0053584214), tolerance = 1e-8)
  expect_equal(sum(a$by_disease$attributable), a$summary$attributable, tolerance = 1e-10)
  expect_equal(a$summary$unallocated, 0, tolerance = 1e-10)
  # Each row is reproduced on its own scale: crude, stratified by LAM and
  # SCK, and stratified by all other diseases.
  beta <- log(res$adjusted$adjusted)
  h0 <- a$baseline_hazard
  j <- res$joint
  got <- c(t2_snapshot_measure(j, beta, 1, kind = "HR"),
           t2_snapshot_measure(j, beta, 2, s = c(1, 4), kind = "HR"),
           t2_snapshot_measure(j, beta, 3, kind = "RR", h0 = h0),
           t2_snapshot_measure(j, beta, 4, kind = "OR", h0 = h0),
           t2_snapshot_measure(j, beta, 5, s = 1:4, kind = "RD", h0 = h0))
  expect_equal(got, c(log(c(1.74, 1.6, 1.23, 1.45)), 0.10), tolerance = 1e-9)
})

test_that("three-way terms change the snapshot model only", {
  m3 <- t2_event(pop = t2_pop(cm_three_way("d1", "d2", "d3", 2)))
  expect_equal(adjust_event(m3, method = "first_order", overall_risk = 0.25)$adjusted$adjusted,
               t2_hr_ref$first_order, tolerance = 1e-8)
  expect_equal(adjust_event(m3, method = "published", overall_risk = 0.25)$adjusted$adjusted,
               t2_hr_ref$published, tolerance = 1e-8)
  snap <- deconflate(m3, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_gt(max(abs(snap$adjusted$adjusted - t2_hr_ref$snapshot)), 1e-3)
  # Independent Python solutions (reference_v020_tests.py; reference_v040.py, section 8).
  expect_equal(snap$adjusted$adjusted, c(1.3798925799054613, 1.8913281873921548, 1.1399165814382808),
               tolerance = 1e-7)
  expect_equal(snap$attributable$summary$attributable, 0.036447342123939075, tolerance = 1e-8)
  # It still reproduces the raw hazard ratios over its own joint.
  beta <- log(snap$adjusted$adjusted)
  expect_equal(exp(vapply(1:3, function(i) t2_snapshot_measure(snap$joint, beta, i), numeric(1))),
               c(1.5, 2.0, 1.3), tolerance = 1e-9)
})

test_that("stratified hazard ratios are matched within strata of their adjustment set", {
  m <- t2_event(estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                adjusted_for = c("d2", NA, NA))
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(res$adjusted$estimand, c("snapshot_stratified", "snapshot_crude", "snapshot_crude"))
  expect_equal(res$adjusted$adjusted_for, c("d2", NA, NA))
  beta <- log(res$adjusted$adjusted)
  expect_equal(t2_snapshot_measure(res$joint, beta, 1, s = 2), log(1.5), tolerance = 1e-9)
  expect_equal(exp(vapply(2:3, function(i) t2_snapshot_measure(res$joint, beta, i), numeric(1))),
               c(2.0, 1.3), tolerance = 1e-9)
  # Independent Python solution of the same equations (reference_v020_tests.py).
  expect_equal(res$adjusted$adjusted, c(1.504714547766569, 1.8751068533634567, 1.1457314393824285),
               tolerance = 1e-7)
  # The published approach is defined for crude hazard ratios only.
  expect_error(adjust_event(m, method = "published", overall_risk = 0.25), class = "deconflate_unsupported")
  expect_error(adjust_event(m, method = "published", overall_risk = 0.25), "crude hazard ratios only")
  # The first-order approximation accepts stratified estimands.
  fo <- adjust_event(m, method = "first_order", overall_risk = 0.25)
  expect_true(all(is.finite(fo$adjusted$adjusted)))
  cmp <- compare_methods(m, event_model = TRUE, overall_risk = 0.25)
  expect_equal(cmp$methods, c("first_order", "snapshot"))
  expect_match(cmp$failed[["published"]], "crude hazard ratios only")
})

test_that("hazard ratios stratified by all other diseases are used as they are", {
  m <- t2_event(estimand = "snapshot_stratified", adjusted_for = "all")
  for (meth in c("snapshot", "first_order")) {
    r <- adjust_event(m, method = meth, overall_risk = 0.25)
    expect_equal(r$adjusted$adjusted, c(1.5, 2.0, 1.3), tolerance = 1e-10)
    expect_equal(r$adjusted$change, c(0, 0, 0), tolerance = 1e-10)
  }
  expect_error(adjust_event(m, method = "published", overall_risk = 0.25), class = "deconflate_unsupported")
  # Mixed: d1 stratified by all, the others crude. d1 keeps its value.
  mix <- deconflate(t2_event(estimand = c("snapshot_stratified", "snapshot_crude", "snapshot_crude"),
                             adjusted_for = c("all", NA, NA)),
                    event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(mix$adjusted$adjusted[1], 1.5, tolerance = 1e-10)
  expect_equal(exp(vapply(2:3, function(i) t2_snapshot_measure(mix$joint, log(mix$adjusted$adjusted), i),
                          numeric(1))),
               c(2.0, 1.3), tolerance = 1e-9)
})

test_that("hazard ratios of 1 stay 1 and give no attributable risk", {
  m <- t2_event(c(1, 1, 1))
  for (meth in c("snapshot", "first_order", "published")) {
    r <- adjust_event(m, method = meth, overall_risk = 0.25)
    expect_equal(r$adjusted$adjusted, c(1, 1, 1), tolerance = 1e-12)
    expect_equal(r$diagnostics$n_sign_changes, 0)
  }
  ar <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)$attributable
  expect_equal(ar$summary$disease_free_risk, 0.25, tolerance = 1e-10)
  expect_equal(ar$summary$attributable, 0, tolerance = 1e-10)
  expect_equal(ar$by_disease$attributable, c(0, 0, 0), tolerance = 1e-12)
  expect_true(all(is.na(ar$by_disease$share)))
})

test_that("the attributable risk and its Shapley allocation match the reference", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  ar <- res$attributable
  expect_equal(names(ar$summary), c("overall_risk", "disease_free_risk", "attributable",
                                    "attributable_fraction", "unallocated"))
  expect_equal(ar$summary$overall_risk, 0.25)
  expect_equal(ar$summary$disease_free_risk, 0.21335101490465225, tolerance = 1e-8)
  expect_equal(ar$summary$attributable, t2_ar_ref, tolerance = 1e-7)
  expect_equal(ar$summary$attributable_fraction, t2_ar_ref / 0.25, tolerance = 1e-7)
  expect_equal(names(ar$by_disease), c("disease", "attributable", "share"))
  expect_equal(ar$by_disease$disease, t2_ids)
  expect_equal(ar$by_disease$attributable, c(0.0073707436, 0.0234877275, 0.0057905139), tolerance = 1e-7)
  expect_equal(sum(ar$by_disease$attributable), ar$summary$attributable, tolerance = 1e-10)
  expect_equal(sum(ar$by_disease$share), 1, tolerance = 1e-12)
  expect_equal(ar$summary$unallocated, 0, tolerance = 1e-10)
  expect_equal(ar$baseline_hazard, -log(1 - ar$summary$disease_free_risk), tolerance = 1e-12)
  # attribute_burden() and contribution_table() show the allocation.
  expect_equal(attribute_burden(res), ar$by_disease)
  ct <- contribution_table(res)
  expect_equal(ct$adjusted_hr, res$adjusted$adjusted)
  expect_equal(ct$attributable, ar$by_disease$attributable)
  expect_equal(ct$measure, rep("HR", 3))

  # Without allocation (internal), and the first-order approximation.
  ea <- event_attributable(log(res$adjusted$adjusted), res$joint, 0.25, allocate = FALSE)
  expect_null(ea$by_disease)
  expect_true(is.na(ea$summary$unallocated))
  expect_equal(ea$summary$attributable, ar$summary$attributable, tolerance = 1e-12)
  fo <- adjust_event(t2_event(), method = "first_order", overall_risk = 0.25)
  expect_equal(sum(fo$attributable$by_disease$attributable), fo$attributable$summary$attributable,
               tolerance = 1e-10)
  # An overall risk given as a distribution is used at its mean.
  rd <- deconflate(t2_event(), event_model = TRUE, overall_risk = dist_beta(250, 750), n_draws = 0)
  expect_equal(rd$overall_risk, 0.25)
  expect_equal(rd$attributable$summary$attributable, ar$summary$attributable, tolerance = 1e-12)
})

test_that("the allocation reports what it cannot allocate", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  beta <- log(res$adjusted$adjusted)
  expect_silent(ea <- event_attributable(beta, res$joint, 0.25, max_present = 2))
  s <- ea$summary
  expect_gt(s$unallocated, 0)
  expect_equal(sum(ea$by_disease$attributable) + s$unallocated, s$attributable, tolerance = 1e-12)
  # The unallocated part is the expected excess risk of animals with all three.
  j <- res$joint
  all3 <- rowSums(j$cells) == 3
  excess <- 1 - exp(-ea$baseline_hazard * prod(res$adjusted$adjusted)) - s$disease_free_risk
  expect_equal(s$unallocated, j$prob[all3] * excess, tolerance = 1e-10)
  expect_equal(event_attributable(beta, res$joint, 0.25, max_present = 3)$summary$unallocated, 0,
               tolerance = 1e-12)
})

test_that("adjusted hazard ratios on the other side of 1 are warned about", {
  # d1's raw hazard ratio 1.05 is weaker than its association with d2
  # (hazard ratio 3) alone would produce (reference_v040.py, section 9).
  m <- t2_event(c(1.05, 3, 1.3))
  expect_warning(res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0),
                 class = "deconflate_sign_change")
  expect_equal(res$adjusted$adjusted, c(0.9009116545, 3.0081830439, 1.0227021972), tolerance = 1e-8)
  expect_equal(res$diagnostics$n_sign_changes, 1)
  expect_equal(res$diagnostics$sign_changes, "d1")
  expect_equal(res$attributable$summary$attributable, 0.04305334120549864, tolerance = 1e-8)
  expect_silent(deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0, warn = FALSE))
  # A positive risk difference with an adjusted hazard ratio below 1.
  rd <- t2_event(c(0.005, 3, 1.3), measure = c("RD", "HR", "HR"))
  expect_warning(r2 <- deconflate(rd, event_model = TRUE, overall_risk = 0.25, n_draws = 0),
                 class = "deconflate_sign_change")
  expect_equal(r2$adjusted$adjusted, c(0.8857581689, 3.0118114622, 1.0224662318), tolerance = 1e-8)
  expect_equal(r2$diagnostics$sign_changes, "d1")
  # compare_methods() flags sign changes without warning.
  # (The published approximation keeps the raw side of 1.)
  expect_silent(cmp <- compare_methods(m, event_model = TRUE, overall_risk = 0.25))
  sc <- stats::setNames(cmp$long$sign_change[cmp$long$disease == "d1"],
                        cmp$long$method[cmp$long$disease == "d1"])
  expect_true(sc[["snapshot"]])
  expect_false(sc[["published"]])
})

test_that("non-positive published hazard ratios are flagged and have no attributable risk", {
  # HR 0.5 for d1 with a strong d2: the published denominator of d1 is
  # negative (reference_v02.py hr_methods).
  m <- t2_event(c(0.5, 4, 1))
  expect_warning(pub <- adjust_event(m, method = "published", overall_risk = 0.25),
                 class = "deconflate_nonfinite_warning")
  expect_equal(pub$adjusted$adjusted, c(-0.3585461449, 4.0376405055, 1), tolerance = 1e-8)
  expect_equal(pub$diagnostics$n_sign_changes, 0)
  expect_null(pub$attributable)
  expect_false(result_is_finite(pub))
  expect_silent(adjust_event(m, method = "published", overall_risk = 0.25, warn = FALSE))
  cmp <- compare_methods(m, event_model = TRUE, overall_risk = 0.25, methods = c("published", "snapshot"))
  expect_equal(cmp$methods, "snapshot")
  expect_match(cmp$failed[["published"]], "^undefined")
  expect_s3_class(cmp$undefined$published, "cm_event_result")
  # The other methods keep the hazard ratios positive.
  expect_equal(adjust_event(m, method = "first_order", overall_risk = 0.25)$adjusted$adjusted,
               c(0.4268305992, 4.4915295004, 0.7780752535), tolerance = 1e-8)
  snap <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(snap$adjusted$adjusted, c(0.4021045208, 4.487174924, 0.7095065051), tolerance = 1e-7)
  expect_equal(snap$attributable$summary$attributable, 0.03568836025623201, tolerance = 1e-6)
})

test_that("unknown pairs are filled in by the snapshot model; the approximations need every pair", {
  pop <- cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                       cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)))
  m <- t2_event(pop = pop)
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(res$method, "snapshot")
  expect_equal(nrow(res$unknown_pairs), 1)
  expect_equal(c(res$unknown_pairs$disease1, res$unknown_pairs$disease2), c("d1", "d3"))
  expect_true(is.finite(res$unknown_pairs$fitted_or))
  # Unknown is not independent: the fit gives d1:d3 an association through d2.
  expect_gt(abs(res$unknown_pairs$fitted_or - 1), 1e-3)
  expect_gt(max(abs(res$adjusted$adjusted - t2_hr_ref$snapshot)), 1e-4)
  for (meth in c("first_order", "published")) {
    expect_error(adjust_event(m, method = meth, overall_risk = 0.25), class = "deconflate_unknown_pairs")
  }
  cmp <- compare_methods(m, event_model = TRUE, overall_risk = 0.25)
  expect_equal(cmp$methods, "snapshot")
  expect_match(cmp$failed[["first_order"]], "every pair")
})

test_that("draws of the overall risk give intervals for the attributable risk only", {
  skip_on_cran()
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = dist_beta(250, 750),
                    n_draws = 20, seed = 1)
  expect_equal(res$draws$seed, 1)
  expect_equal(res$draws$n_draws, 20L)
  expect_equal(res$draws$n_rejected, 0L)
  expect_true(all(c("adjusted:d1", "contribution:d1", "total", "fraction") %in% res$draws$summary$quantity))
  # The central estimate is the point estimate at the mean overall risk.
  expect_equal(res$attributable$summary$attributable, t2_ar_ref, tolerance = 1e-7)
  # Hazard ratios do not depend on the overall risk: no spread.
  expect_equal(res$adjusted$lower, res$adjusted$adjusted, tolerance = 1e-8)
  expect_equal(res$adjusted$upper, res$adjusted$adjusted, tolerance = 1e-8)
  s <- res$attributable$summary
  expect_lt(s$attributable_lower, s$attributable)
  expect_gt(s$attributable_upper, s$attributable)
  expect_true(all(c("lower", "upper") %in% names(res$attributable$by_disease)))
})

test_that("global dairy culling: unknown pairs, overall risk 0.2366", {
  skip_on_cran()
  m <- example_global_dairy("culling")
  expect_equal(unique(m$impacts$measure), "HR")
  expect_equal(unique(m$impacts$estimand), "snapshot_crude")
  expect_warning(res <- deconflate(m, event_model = TRUE, overall_risk = 0.2366, n_draws = 0),
                 class = "deconflate_sign_change")
  expect_equal(res$adjusted$disease, c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM"))
  expect_equal(res$adjusted$adjusted,
               c(0.9771, 1.818179, 1.786972, 0.954031, 1.240105, 0.69279, 2.42472, 1.40223, 1.920838,
                 1.183173, 1.645226, 1.154687), tolerance = 1e-5)
  expect_equal(res$diagnostics$sign_changes, "CK, DYS, MET")
  expect_equal(res$attributable$summary$attributable, 0.1086401602, tolerance = 1e-7)
  expect_equal(res$attributable$summary$attributable_fraction, 0.4591722747, tolerance = 1e-7)
  expect_equal(sum(res$attributable$by_disease$attributable), res$attributable$summary$attributable,
               tolerance = 1e-10)
  expect_equal(nrow(res$unknown_pairs), 28)
  expect_true(all(is.finite(res$unknown_pairs$fitted_or)))
  expect_length(res$notes, 0)
})

test_that("global dairy culling with unlisted pairs independent: the three methods (as in 0.3)", {
  skip_on_cran()
  m0 <- example_global_dairy("culling")
  m <- cm_model(with_independent_pairs(as_population(m0)), m0$impacts)
  expect_equal(nrow(pair_tables(m)), 66)
  expect_true(all(pair_tables(m)$status == "specified"))
  pub <- adjust_event(m, method = "published", overall_risk = 0.2366, warn = FALSE, allocate = FALSE)
  expect_equal(pub$adjusted$adjusted,
               c(1.1775803801, 1.9039406263, 2.1979299277, 1.0983766711, 1.3806827011, 1.0124105823,
                 2.6476374499, 1.4586441181, 2.0472349212, 1.2844963465, 1.6752532746, 1.2549275636),
               tolerance = 1e-8)
  # Identical to adjusting HR - 1 additively with the published method, as
  # in Rasmussen et al. (2024).
  legacy <- adjust_impacts(global_dairy_analyses("analysis", culling = TRUE)$culling_hr_minus_1,
                           method = "published", warn = FALSE)
  expect_equal(pub$adjusted$adjusted, 1 + legacy$adjusted$adjusted, tolerance = 1e-12)

  expect_warning(adjust_event(m, method = "first_order", overall_risk = 0.2366, allocate = FALSE),
                 class = "deconflate_sign_change")
  fo <- adjust_event(m, method = "first_order", overall_risk = 0.2366, warn = FALSE, allocate = FALSE)
  expect_equal(fo$adjusted$adjusted,
               c(0.9619070191, 1.8350863583, 1.825054203, 1.1165198124, 1.2864922835, 0.7369878787,
                 2.6233669144, 1.5487477113, 2.0203081931, 1.2305973733, 1.7412069142, 1.2520841682),
               tolerance = 1e-8)
  expect_equal(fo$diagnostics$sign_changes, "CK, MET")

  snap <- adjust_event(m, method = "snapshot", overall_risk = 0.2366, warn = FALSE, allocate = FALSE)
  expect_equal(snap$adjusted$adjusted,
               c(0.9821306495, 1.8193535917, 1.8476664195, 1.050543232, 1.2230209916, 0.7198426768,
                 2.6601485634, 1.5579660463, 2.0165455108, 1.1816994757, 1.7107120034, 1.210900523),
               tolerance = 1e-6)
  beta <- log(snap$adjusted$adjusted)
  expect_equal(exp(vapply(1:12, function(i) t2_snapshot_measure(snap$joint, beta, i), numeric(1))),
               m$impacts$value, tolerance = 1e-8)
  expect_equal(snap$diagnostics$n_sign_changes, 2)
  expect_equal(snap$attributable$summary$disease_free_risk, 0.12002780924706569, tolerance = 1e-6)
  expect_equal(snap$attributable$summary$attributable, 0.11657219075293432, tolerance = 1e-6)
  expect_null(snap$attributable$by_disease)
})

test_that("legacy conversions reproduce the published culling inputs", {
  # Rasmussen et al. (2022), Table 6: hazard ratio treated as an odds ratio
  # with an overall culling rate of 0.27 (DA: 3.83, P = 0.03; LAM: 3.40,
  # P = 0.30).
  ex <- legacy_hr_as_or_excess(c(3.83, 3.40, 1), c(0.03, 0.30, 0.10), 0.27)
  expect_equal(ex, c(0.3138414572431465, 0.2556484948872295, 0), tolerance = 1e-10)
  uk <- uk_dairy_2022_analyses(3.05, culling = TRUE)
  cu <- uk$culling$impacts
  expect_equal(cu$value[cu$disease %in% c("DA", "LAM")], ex[1:2], tolerance = 1e-12)
  expect_equal(cu$value[cu$disease == "CO"], 0)
  expect_equal(impact_kind(cu), "additive")
  h <- attr(uk, "hazard_ratios")
  expect_equal(h[["DA"]], 3.83)
  expect_equal(h[["CO"]], 1)
  # The same hazard ratios are the event impacts of the culling example.
  ukc <- example_uk_dairy_2022("culling")
  expect_equal(ukc$impacts$value, unname(h[ukc$impacts$disease]))
  expect_equal(unique(ukc$impacts$measure), "HR")
  # Eq. 23 rescaling, and the 2024 overall-odds excess.
  expect_equal(legacy_eq23(c(3.83, 1), c(0.31, 0), c(0.2, 0)), c(0.2 * 3.83 / 0.31, 1))
  expect_equal(legacy_overall_odds_excess(c(2.75, 1), 0.27), c(0.23424448217317484, 0),
               tolerance = 1e-12)
})
