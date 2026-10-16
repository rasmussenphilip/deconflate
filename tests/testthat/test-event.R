# Event impacts (R/hazard.R): cm_impacts(measure = ...) and
# deconflate(event_model = TRUE, overall_risk = ...) with the snapshot hazard
# model (internal adjust_event()); and the attributable risk and its Shapley
# allocation. Reference values: inst/validation/reference_v02.py (its
# "global" is snapshot), reference_v020_tests.py, reference_v040.py
# (sections 2, 4, 5, 8, 9) and reference_v041.py (the five-disease culling
# table).

t2_ids <- c("d1", "d2", "d3")

# Supplementary File population (Rasmussen et al. 2022): P = 0.10, 0.15,
# 0.20; OR d1:d2 = 2, d1:d3 = 1, d2:d3 = 3.
t2_pop <- function(extra = NULL) supp_population(extra)

# Event model on a population (default: culling hazard ratios 1.5, 2.0, 1.3).
t2_event <- function(values = c(1.5, 2.0, 1.3), measure = "HR", adjusted_for = NA_character_, pop = t2_pop()) {
  cm_model(pop, cm_impacts(t2_ids, values, measure = measure, adjusted_for = adjusted_for, label = "culling"))
}

# The snapshot model on the model's own joint distribution (internal).
t2_snap <- function(m, risk = 0.25, ...) adjust_event(m, overall_risk = risk, joint = fit_joint(m), ...)
one_term <- function(ids) stats::setNames(as.list(ids), ids)

t2_hr_ref <- list(snapshot = c(1.3831529713, 1.8909315489, 1.1439870467))
t2_ar_ref <- 0.03664898509534775   # attributable risk of the snapshot solution at 0.25

# The snapshot measure of disease i over a joint distribution, on the
# measure's own scale (log for ratios, the difference for RD): crude, or
# within strata of the diseases in columns s combined with weights
# d1 d0 / (d1 + d0). h0: baseline hazard (risk-based measures).
t2_snapshot_measure <- function(joint, beta, i, s = integer(0), kind = "HR", h0 = NULL, lp = NULL) {
  mult <- exp(if (is.null(lp)) as.vector(joint$cells %*% beta) else lp)
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

test_that("event impacts are validated by cm_impacts()", {
  ev <- function(...) cm_impacts(..., measure = "HR")
  expect_error(ev("d1", 0), "positive")
  expect_error(ev(t2_ids, c(1.5, 2)), "one entry per row")
  expect_error(cm_impacts("d1", 1.5, measure = "HR", estimand = "snapshot_crude"), "no longer used",
               class = "deconflate_unsupported")
  expect_error(cm_impacts("d1", 1.5, measure = "hazard"), "measure")
  # Ratios are positive; risk differences lie in (-1, 1).
  for (m in c("HR", "rate_ratio", "RR", "OR")) {
    expect_error(cm_impacts("d1", -0.5, measure = m), "positive")
  }
  expect_error(cm_impacts("d1", 1, measure = "RD"), "between -1 and 1")
  expect_error(cm_impacts("d1", -1, measure = "RD"), "between -1 and 1")
  expect_silent(cm_impacts(c("d1", "d2"), c(-0.2, 0), measure = "RD"))

  h <- cm_impacts(t2_ids, c(1.5, 0.04, 1.3), measure = c("HR", "RD", "OR"))
  expect_s3_class(h, "cm_impacts")
  expect_equal(attr(h, "kind"), "event")
  expect_equal(h$measure, c("HR", "RD", "OR"))
  expect_equal(impact_kind(h), "event")
  expect_equal(attr(cm_impacts(t2_ids, c(1, 2, 3)), "kind"), "additive")
  expect_true(all(is.na(cm_impacts(t2_ids, c(1, 2, 3))$measure)))
  # Several rows per disease are allowed (fitted together).
  expect_s3_class(ev(c("d1", "d1"), c(1.5, 2)), "cm_impacts")
})

test_that("event models are validated by cm_model()", {
  ev <- function(ids, v, ...) cm_impacts(ids, v, measure = "HR", ...)
  expect_error(cm_model(t2_pop(), ev(c("d1", "d2"), c(1.5, 2))), "No impact for: d3")
  expect_error(cm_model(t2_pop(), ev(c("d1", "d2"), c(1.5, 2))), "1 for ratios, 0 for risk differences")
  expect_error(cm_model(t2_pop(), ev(c(t2_ids, "zz"), c(1.5, 2, 1.3, 1))), "unknown disease")
  expect_error(cm_model(t2_pop(), ev(c(t2_ids, "d1;zz"), c(1.5, 2, 1.3, 1))), "unknown disease")
  # Rows keep their order; results are in disease order.
  m <- cm_model(t2_pop(), cm_impacts(rev(t2_ids), c(0.05, 2.0, 1.5), measure = c("RD", "HR", "OR")))
  expect_equal(m$impacts$disease, rev(t2_ids))
  expect_equal(impact_kind(m$impacts), "event")
  r <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(r$adjusted$disease, t2_ids)
  expect_equal(r$adjusted$measure, c("OR", "HR", "RD"))
  # Interactions: combination rows; set_interaction() adds an HR row adjusted for all.
  expect_equal(utils::tail(set_interaction(t2_event(), "d1", "d2", 0.5)$impacts$measure, 1), "HR")
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
  expect_error(adjust_event(example_supplement(), overall_risk = 0.25, joint = fit_joint(example_supplement())),
               class = "deconflate_unsupported")
  expect_error(cm_sensitivity(m), class = "deconflate_unsupported")
  expect_error(cm_sensitivity(example_supplement(), event_model = TRUE, overall_risk = 0.25),
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
  expect_match(printed(res), "contribution share", fixed = TRUE)
  expect_match(printed(res), "attributable to disease", fixed = TRUE)
  # There is no method argument.
  expect_error(deconflate(t2_event(), method = "global", event_model = TRUE, overall_risk = 0.25),
               "has been removed")
})

test_that("the internal snapshot engine matches the reference", {
  r <- t2_snap(t2_event())
  expect_s3_class(r, "cm_event_result")
  expect_equal(r$method, "snapshot")
  expect_equal(r$adjusted$adjusted, t2_hr_ref$snapshot, tolerance = 1e-7)
  expect_equal(r$diagnostics$n_sign_changes, 0)
  expect_s3_class(r$joint, "cm_joint")
  expect_equal(r$attributable$summary$attributable, t2_ar_ref, tolerance = 1e-7)
  expect_equal(sum(r$contributions$contribution), r$attributable$summary$attributable,
               tolerance = 1e-10)
})

test_that("the snapshot model reproduces the raw hazard ratios over the fitted joint", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  beta <- log(res$adjusted$adjusted)
  expect_equal(exp(vapply(1:3, function(i) t2_snapshot_measure(res$joint, beta, i), numeric(1))),
               c(1.5, 2.0, 1.3), tolerance = 1e-9)
  # A supplied joint gives the same answer; a joint of another population
  # (here with an OR_ratio row) is rejected.
  j <- fit_joint(t2_pop())
  expect_equal(deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, joint = j)$adjusted$adjusted,
               res$adjusted$adjusted, tolerance = 1e-12)
  j3 <- fit_joint(t2_pop(cm_associations("d1", "d2", 2, measure = "OR_ratio", among = "d3")))
  expect_error(deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, joint = j3), "association rows")
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
  # The change is on each row's measure: OR, RR and RD alone vs disease-free.
  expect_equal(re$adjusted$change, c(re$adjusted$or_alone[1], re$adjusted$rr_alone[2], re$adjusted$rd_alone[3]) /
                 re$adjusted$raw - 1)
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
})

test_that("mixed measures, crude, stratified and a combination match the five-disease reference", {
  m <- read_five("culling.csv")
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(res$adjusted$disease, c("LAM", "MAS", "MET", "SCK", "RP", "LAM;MAS"))
  expect_equal(res$adjusted$kind, c(rep("disease", 5), "interaction"))
  expect_equal(res$adjusted$measure, c("HR", "HR", "RR", "OR", "RD", "HR"))
  # Python reference (reference_v041.py).
  expect_equal(res$adjusted$adjusted, c(1.060298968, 1.2342866774, 1.080381173, 1.1753668593, 1.509273874,
                                        2.2779452712), tolerance = 1e-7)
  expect_equal(is.na(res$adjusted$change), c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(res$estimates$adjusted_for[1:6], c(NA, "LAM;SCK", NA, NA, "all", NA))
  a <- res$attributable
  expect_equal(a$summary$attributable, 0.0541469811, tolerance = 1e-8)
  expect_equal(res$contributions$contribution,
               c(0.0137316223, 0.0218514738, 0.0015789582, 0.0114739881, 0.0055109386), tolerance = 1e-8)
  expect_equal(sum(res$contributions$contribution), a$summary$attributable, tolerance = 1e-10)
  expect_equal(a$summary$unallocated, 0, tolerance = 1e-10)
  # Each single-disease row is reproduced on its own scale: crude, stratified
  # by LAM and SCK, and stratified by all other diseases; with the LAM:MAS
  # interaction in the linear predictor.
  j <- res$joint
  lp <- as.vector(term_matrix(j$cells, res$terms) %*% log(res$adjusted$adjusted))
  h0 <- a$baseline_hazard
  got <- c(t2_snapshot_measure(j, NULL, 1, kind = "HR", lp = lp),
           t2_snapshot_measure(j, NULL, 2, s = c(1, 4), kind = "HR", lp = lp),
           t2_snapshot_measure(j, NULL, 3, kind = "RR", h0 = h0, lp = lp),
           t2_snapshot_measure(j, NULL, 4, kind = "OR", h0 = h0, lp = lp),
           t2_snapshot_measure(j, NULL, 5, s = 1:4, kind = "RD", h0 = h0, lp = lp))
  expect_equal(got, c(log(c(1.74, 1.6, 1.23, 1.45)), 0.10), tolerance = 1e-9)
  # The combination row: cows with both LAM and MAS against cows with neither.
  both <- j$cells[, "LAM"] == 1 & j$cells[, "MAS"] == 1
  neither <- j$cells[, "LAM"] == 0 & j$cells[, "MAS"] == 0
  mult <- exp(lp)
  hr <- (sum(j$prob[both] * mult[both]) / sum(j$prob[both])) / (sum(j$prob[neither] * mult[neither]) / sum(j$prob[neither]))
  expect_equal(hr, 3.2, tolerance = 1e-9)
})

test_that("OR_ratio rows change the snapshot model", {
  # With three diseases, an OR_ratio row is the log-linear three-way term of
  # earlier versions.
  m3 <- t2_event(pop = t2_pop(cm_associations("d1", "d2", 2, measure = "OR_ratio", among = "d3")))
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
  m <- t2_event(adjusted_for = c("d2", NA, NA))
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(res$estimates$adjusted_for[1:3], c("d2", NA, NA))
  beta <- log(res$adjusted$adjusted)
  expect_equal(t2_snapshot_measure(res$joint, beta, 1, s = 2), log(1.5), tolerance = 1e-9)
  expect_equal(exp(vapply(2:3, function(i) t2_snapshot_measure(res$joint, beta, i), numeric(1))),
               c(2.0, 1.3), tolerance = 1e-9)
  # Independent Python solution of the same equations (reference_v020_tests.py).
  expect_equal(res$adjusted$adjusted, c(1.504714547766569, 1.8751068533634567, 1.1457314393824285),
               tolerance = 1e-7)
})

test_that("hazard ratios adjusted for all other diseases are used as they are", {
  m <- t2_event(adjusted_for = "all")
  r <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(r$adjusted$adjusted, c(1.5, 2.0, 1.3), tolerance = 1e-10)
  expect_equal(r$adjusted$change, c(0, 0, 0), tolerance = 1e-10)
  # Mixed: d1 adjusted for all, the others crude. d1 keeps its value.
  mix <- deconflate(t2_event(adjusted_for = c("all", NA, NA)),
                    event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(mix$adjusted$adjusted[1], 1.5, tolerance = 1e-10)
  expect_equal(exp(vapply(2:3, function(i) t2_snapshot_measure(mix$joint, log(mix$adjusted$adjusted), i),
                          numeric(1))),
               c(2.0, 1.3), tolerance = 1e-9)
})

test_that("hazard ratios of 1 stay 1 and give no attributable risk", {
  m <- t2_event(c(1, 1, 1))
  r <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(r$adjusted$adjusted, c(1, 1, 1), tolerance = 1e-12)
  expect_equal(r$diagnostics$n_sign_changes, 0)
  ar <- r$attributable
  expect_equal(ar$summary$disease_free_risk, 0.25, tolerance = 1e-10)
  expect_equal(ar$summary$attributable, 0, tolerance = 1e-10)
  expect_equal(r$contributions$contribution, c(0, 0, 0), tolerance = 1e-12)
  expect_true(all(is.na(r$contributions$share)))
})

test_that("the attributable risk and its Shapley allocation match the reference", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  ar <- res$attributable
  expect_equal(names(ar$summary), c("overall_risk", "disease_free_risk", "attributable",
                                    "attributable_fraction", "unallocated", "raw_attributable"))
  expect_equal(ar$summary$overall_risk, 0.25)
  expect_equal(ar$summary$disease_free_risk, 0.21335101490465225, tolerance = 1e-8)
  expect_equal(ar$summary$attributable, t2_ar_ref, tolerance = 1e-7)
  expect_equal(ar$summary$attributable_fraction, t2_ar_ref / 0.25, tolerance = 1e-7)
  expect_equal(names(ar), c("summary", "baseline_hazard"))
  ctb <- res$contributions
  expect_equal(names(ctb), c("disease", "contribution", "main", "interaction", "share"))
  expect_equal(ctb$disease, t2_ids)
  expect_equal(ctb$contribution, c(0.0073707436, 0.0234877275, 0.0057905139), tolerance = 1e-7)
  expect_equal(sum(ctb$contribution), ar$summary$attributable, tolerance = 1e-10)
  expect_equal(sum(ctb$share), 1, tolerance = 1e-12)
  expect_equal(ar$summary$unallocated, 0, tolerance = 1e-10)
  expect_equal(ar$baseline_hazard, -log(1 - ar$summary$disease_free_risk), tolerance = 1e-12)
  # contribution_table() shows the allocation.
  ct <- contribution_table(res)
  expect_equal(ct$adjusted, res$adjusted$adjusted)
  expect_equal(ct$contribution, ctb$contribution)
  expect_equal(ct$measure, rep("HR", 3))

  # Without allocation (internal).
  ea <- event_attributable(log(res$adjusted$adjusted), res$joint, 0.25, one_term(t2_ids), allocate = FALSE)
  expect_null(ea$contributions)
  expect_true(is.na(ea$summary$unallocated))
  expect_equal(ea$summary$attributable, ar$summary$attributable, tolerance = 1e-12)
  # An overall risk given as a distribution is used at its mean.
  rd <- deconflate(t2_event(), event_model = TRUE, overall_risk = dist_beta(250, 750), n_draws = 0)
  expect_equal(rd$overall_risk, 0.25)
  expect_equal(rd$attributable$summary$attributable, ar$summary$attributable, tolerance = 1e-12)
})

test_that("the allocation reports what it cannot allocate", {
  res <- deconflate(t2_event(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  beta <- log(res$adjusted$adjusted)
  expect_silent(ea <- event_attributable(beta, res$joint, 0.25, one_term(t2_ids), max_present = 2))
  s <- ea$summary
  expect_gt(s$unallocated, 0)
  expect_equal(sum(ea$contributions$contribution) + s$unallocated, s$attributable, tolerance = 1e-12)
  # The unallocated part is the expected excess risk of animals with all three.
  j <- res$joint
  all3 <- rowSums(j$cells) == 3
  excess <- 1 - exp(-ea$baseline_hazard * prod(res$adjusted$adjusted)) - s$disease_free_risk
  expect_equal(s$unallocated, j$prob[all3] * excess, tolerance = 1e-10)
  expect_equal(event_attributable(beta, res$joint, 0.25, one_term(t2_ids), max_present = 3)$summary$unallocated, 0,
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
})

test_that("a hazard ratio below 1 next to a strong associated disease stays positive", {
  # HR 0.5 for d1 with a strong d2 (reference_v02.py hr_methods).
  m <- t2_event(c(0.5, 4, 1))
  snap <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(snap$adjusted$adjusted, c(0.4021045208, 4.487174924, 0.7095065051), tolerance = 1e-7)
  expect_equal(snap$attributable$summary$attributable, 0.03568836025623201, tolerance = 1e-6)
})

test_that("unknown pairs are filled in by the snapshot model", {
  m <- t2_event(pop = supp_unknown_population())
  res <- deconflate(m, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(res$method, "snapshot")
  expect_equal(nrow(res$unknown_pairs), 1)
  expect_equal(c(res$unknown_pairs$disease1, res$unknown_pairs$disease2), c("d1", "d3"))
  expect_true(is.finite(res$unknown_pairs$fitted_or))
  # Unknown is not independent: the fit gives d1:d3 an association through d2.
  expect_gt(abs(res$unknown_pairs$fitted_or - 1), 1e-3)
  expect_gt(max(abs(res$adjusted$adjusted - t2_hr_ref$snapshot)), 1e-4)
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
  expect_true(all(c("lower", "upper") %in% names(res$contributions)))
})

test_that("global dairy culling: unknown pairs, overall risk 0.2366", {
  skip_on_cran()
  m <- example_global_dairy("culling")
  expect_equal(unique(m$impacts$measure), "HR")
  expect_warning(res <- deconflate(m, event_model = TRUE, overall_risk = 0.2366, n_draws = 0),
                 class = "deconflate_sign_change")
  expect_equal(res$adjusted$disease, c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM"))
  expect_equal(res$adjusted$adjusted,
               c(0.9771, 1.818179, 1.786972, 0.954031, 1.240105, 0.69279, 2.42472, 1.40223, 1.920838,
                 1.183173, 1.645226, 1.154687), tolerance = 1e-5)
  expect_equal(res$diagnostics$sign_changes, "CK, DYS, MET")
  expect_equal(res$attributable$summary$attributable, 0.1086401602, tolerance = 1e-7)
  expect_equal(res$attributable$summary$attributable_fraction, 0.4591722747, tolerance = 1e-7)
  expect_equal(sum(res$contributions$contribution), res$attributable$summary$attributable,
               tolerance = 1e-10)
  expect_equal(nrow(res$unknown_pairs), 28)
  expect_true(all(is.finite(res$unknown_pairs$fitted_or)))
  # The only note: the inputs have distributions but no draws were run.
  expect_length(res$notes, 1)
  expect_match(res$notes, "but n_draws = 0: the results are point estimates without intervals", fixed = TRUE)
})

test_that("global dairy culling with unlisted pairs independent (as in 0.3)", {
  skip_on_cran()
  m0 <- example_global_dairy("culling")
  # Every unlisted pair gets an OR = 1 row.
  pop <- as_population(m0)
  pt <- pair_tables(pop)
  u <- pt[pt$status == "unknown", , drop = FALSE]
  add <- cm_associations(u$disease1, u$disease2, 1, measure = "OR")
  add$row <- max(pop$associations$row) + seq_len(nrow(add))
  pop$associations <- rbind(pop$associations, add)
  m <- cm_model(pop, m0$impacts)
  expect_equal(nrow(pair_tables(m)), 66)
  expect_true(all(pair_tables(m)$status == "specified"))
  snap <- t2_snap(m, 0.2366, warn = FALSE, allocate = FALSE)
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
  expect_null(snap$contributions)
})
