# Estimands, identifiability and edge cases of the additive engine.
# Reference values: inst/validation/reference_v020.py.

t2_ids <- c("d1", "d2", "d3")
truth <- c(0.02, 0.04, 0.06)

t2_pop <- function(extra = NULL) supp_population(extra)
t2_model <- function(values, adjusted_for = NA_character_, interactions = NULL, units = NULL) {
  supp_model(values, adjusted_for = adjusted_for, interactions = interactions, units = units)
}
# Bob's infeasible triple: p = 0.5, ORs 20 (a:b), 20 (a:c), 0.05 (b:c).
t2_bad_pop <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

test_that("crude estimates: the simultaneous solution recovers the known truth", {
  pop <- t2_pop()
  raw <- simulate_raw_impacts(pop, stats::setNames(truth, t2_ids))
  expect_true(all(is.na(raw$adjusted_for)))
  res <- deconflate(cm_model(pop, raw))
  expect_equal(res$method, "pairwise shortcut")
  expect_equal(res$adjusted$adjusted, truth, tolerance = 1e-9)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-12)
})

test_that("Bob's counterexample: an adjusted coefficient needs the projection mapping", {
  # d1's estimate is the coefficient of d1 in an additive regression adjusted
  # for d2 (in the same population); d2 and d3 are crude.
  raw <- c(0.018664808346, 0.054064376914, 0.066681750858)
  m <- t2_model(raw, adjusted_for = c("d2", NA, NA))
  res <- deconflate(m)
  expect_equal(res$adjusted$adjusted, truth, tolerance = 1e-9)
  expect_equal(res$estimates$adjusted_for[1:3], c("d2", NA, NA))
  A <- res$conflation$A
  expect_equal(rownames(A), c("impacts:1", "impacts:2", "impacts:3"))
  expect_equal(unname(A[1, "d2"]), 0)
  expect_equal(unname(A[1, "d1"]), 1)
  # Zeroing the adjusted-for entry of the crude row (the shortcut) is wrong.
  shortcut <- c(0.01866481, 0.04010287, 0.05998282)
  expect_gt(max(abs(res$adjusted$adjusted - shortcut)), 1e-3)
  # The full model gives the same answer.
  expect_equal(adjust_additive(m, joint = fit_joint(m))$adjusted$adjusted, truth, tolerance = 1e-8)
  # The raw value matches a weighted least-squares fit on the joint
  # distribution (simulate_raw_impacts).
  sim <- simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids),
                              adjusted_for = c("d2", NA, NA))
  expect_equal(sim$value, raw, tolerance = 1e-9)
})

test_that("treating an adjusted coefficient as crude gives the wrong answer", {
  raw <- c(0.018664808346, 0.054064376914, 0.066681750858)
  res <- deconflate(t2_model(raw))
  expect_gt(max(abs(res$adjusted$adjusted - truth)), 1e-3)
})

test_that("fully adjusted coefficients need no adjustment", {
  raw <- c(0.021, 0.039, 0.062)
  res <- deconflate(t2_model(raw, adjusted_for = "all"))
  expect_equal(res$adjusted$adjusted, raw, tolerance = 1e-14)
  expect_equal(unname(res$conflation$A), diag(3))
  # Listing every other disease is the same as "all".
  res2 <- deconflate(t2_model(raw, adjusted_for = c("d2; d3", "d1; d3", "d1;d2")))
  expect_equal(res2$adjusted$adjusted, raw, tolerance = 1e-14)
})

test_that("mixed estimands recover the truth", {
  raw <- c(0.02, 0.054064376914, 0.066681750858)
  m <- t2_model(raw, adjusted_for = c("all", NA, "d1"))
  expect_equal(deconflate(m)$adjusted$adjusted, truth, tolerance = 1e-9)
  sim <- simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids),
                              adjusted_for = c("all", NA, "d1"))
  expect_equal(sim$value, raw, tolerance = 1e-9)
})

test_that("the estimand of earlier versions is replaced by adjusted_for", {
  expect_error(cm_impacts(t2_ids, truth, estimand = "adjusted_linear"), "no longer used",
               class = "deconflate_unsupported")
  expect_error(simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids), estimand = "crude"),
               "no longer used")
  # A name that is not a disease is a factor outside the model: the row is
  # used as crude, with a note.
  m <- t2_model(c(2.5, 5, 7.5), adjusted_for = c("parity", NA, NA))
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$adjusted$adjusted, deconflate(example_supplement())$adjusted$adjusted)
  expect_match(res$notes, "row 1 is adjusted only for factors outside the model \\(parity\\); used as crude",
               all = FALSE)
  m2 <- t2_model(c(2.5, 5, 7.5), adjusted_for = c("d2; parity", NA, NA))
  r2 <- deconflate(m2, n_draws = 0)
  expect_equal(r2$adjusted$adjusted, deconflate(t2_model(c(2.5, 5, 7.5), adjusted_for = c("d2", NA, NA)))$adjusted$adjusted)
  expect_match(r2$notes, "is ignored", all = FALSE)
})

test_that("a non-identifiable adjustment set is rejected", {
  # b is present exactly when a is: a regression on a and b cannot separate
  # them.
  pop <- cm_population(cm_diseases(c("a", "b"), c(0.2, 0.2)),
                       cm_associations("a", "b", 1 - 1e-12, measure = "prevalence"))
  adj <- cm_model(pop, cm_impacts(c("a", "b"), c(1, 2), adjusted_for = c("b", NA)))
  expect_error(deconflate(adj), class = "deconflate_singular")
  crude <- cm_model(pop, cm_impacts(c("a", "b"), c(1, 2)))
  expect_error(deconflate(crude), class = "deconflate_singular")
})

test_that("zero impacts give a zero aggregate with undefined shares", {
  res <- deconflate(t2_model(c(0, 0, 0)))
  expect_equal(res$adjusted$adjusted, c(0, 0, 0))
  expect_equal(res$totals$adjusted_total, 0)
  expect_true(all(is.na(res$contributions$share)))
  expect_true(all(is.na(res$adjusted$change)))
  expect_equal(contribution_table(res)$contribution, c(0, 0, 0))
})

test_that("cancelling signed impacts are handled without dividing by the total", {
  pop <- cm_population(cm_diseases(c("a", "b"), c(0.3, 0.3)), cm_associations("a", "b", 1))
  res <- deconflate(cm_model(pop, cm_impacts(c("a", "b"), c(5, -5), units = "kg")))
  expect_equal(res$adjusted$adjusted, c(5, -5))
  expect_equal(res$totals$adjusted_total, 0)
  expect_true(all(is.na(res$contributions$share)))
  ct <- contribution_table(res)
  expect_equal(ct$contribution, c(1.5, -1.5))
  expect_true(all(is.na(ct$share)))
})

test_that("jointly infeasible pairs are rejected although A is well conditioned", {
  m <- cm_model(t2_bad_pop(), cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(deconflate(m), class = "deconflate_infeasible")
  expect_error(deconflate(m), "conflicting triples: a-b-c")
  # Unchecked, the infeasible inputs flip the signs of b and c.
  res <- adjust_additive(m, pt = pair_tables(m), warn = FALSE)
  expect_equal(res$diagnostics$sign_changes, "b, c")
  expect_equal(res$diagnostics$condition_number, 6.075710520908865, tolerance = 1e-8)
  expect_warning(j <- fit_joint(m, max_iter = 300), class = "deconflate_nonconvergence")
  expect_false(j$converged)
  expect_error(deconflate(m, joint = j), class = "deconflate_infeasible")
  expect_error(validate_joint(j, m), class = "deconflate_nonconvergence")
})

test_that("a stale joint distribution is rejected; an impact change is not", {
  m <- example_supplement()
  j <- fit_joint(m)
  expect_silent(deconflate(m, joint = j, n_draws = 0))
  m2 <- set_association(m, "d1", "d2", 3)
  expect_error(deconflate(m2, joint = j), "association values", class = "deconflate_error")
  m3 <- m
  m3$diseases$prob[1] <- 0.11
  expect_error(deconflate(m3, joint = j), "disease probabilities", class = "deconflate_error")
  m4 <- cm_model(t2_pop(cm_associations("d1", "d2", 2, measure = "OR_ratio", among = "d3")), m$impacts)
  expect_error(deconflate(m4, joint = j), "different association rows", class = "deconflate_error")
  m6 <- m
  m6$associations <- m$associations[-2, ]
  expect_error(deconflate(m6, joint = j), "different association rows")
  m5 <- t2_model(c(1, 2, 3))
  expect_equal(deconflate(m5, joint = j)$adjusted$adjusted,
               adjust_additive(m5, joint = j)$adjusted$adjusted, tolerance = 1e-8)
  # The pooling is part of the signature.
  expect_error(deconflate(m, joint = j, pool = "fixed"), "pooling")
})

test_that("interactions: the model recovers the truth", {
  # True impacts 0.02, 0.04, 0.06; interactions d1:d2 = 0.01, d2:d3 = 0.015,
  # entered as rows adjusted for all (the interaction itself).
  raw <- c(0.027146191901, 0.061369036759, 0.071166773469)
  m <- t2_model(raw, interactions = c("d1;d2" = 0.01, "d2;d3" = 0.015))
  expect_false(shortcut_ok(m))
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "full model")
  expect_equal(res$adjusted$disease, c(t2_ids, "d1;d2", "d2;d3"))
  expect_equal(res$adjusted$kind, c(rep("disease", 3), rep("interaction", 2)))
  expect_equal(res$adjusted$adjusted, c(truth, 0.01, 0.015), tolerance = 1e-8)
  expect_equal(res$totals$adjusted_total, 0.021095698976771753, tolerance = 1e-9)
  expect_equal(res$contributions$contribution, c(0.002122396963, 0.006547849488, 0.012425452526),
               tolerance = 1e-9)
  expect_equal(sum(res$contributions$contribution), res$totals$adjusted_total, tolerance = 1e-14)
  expect_lt(res$diagnostics$max_reconstruction_residual, 1e-10)
  # raw = A theta, row by row.
  expect_equal(as.vector(res$conflation$A %*% res$adjusted$adjusted), c(raw, 0.01, 0.015),
               tolerance = 1e-9)
  # The interaction interface of earlier versions gives a message.
  expect_error(cm_model(t2_pop(), m$impacts, interactions = data.frame()), "rows of the impact table",
               class = "deconflate_unsupported")
})

test_that("interactions with adjusted estimates are recovered", {
  int <- c("d1;d2" = 0.01, "d2;d3" = 0.015)
  adj <- c("d2", NA, "all")
  raw <- simulate_raw_impacts(t2_pop(), stats::setNames(truth, t2_ids),
                              interactions = int, adjusted_for = adj)
  expect_equal(raw$disease, c(t2_ids, "d1;d2", "d2;d3"))
  res <- deconflate(cm_model(t2_pop(), raw))
  expect_equal(res$method, "full model")
  expect_equal(res$adjusted$adjusted, c(truth, int), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("an OR_ratio row sets the stratum odds ratios and keeps the pair tables", {
  pop <- t2_pop(cm_associations("d1", "d2", 2.5, measure = "OR_ratio", among = "d3"))
  j <- fit_joint(pop)
  expect_true(j$converged)
  p <- j$prob
  x <- j$cells
  sor <- function(s) {
    A <- x[, 3] == s
    pa <- sum(p[A])
    joint_to_or(sum(p[A & x[, 1] == 1 & x[, 2] == 1]) / pa, sum(p[A & x[, 1] == 1]) / pa,
                sum(p[A & x[, 2] == 1]) / pa)
  }
  expect_equal(sor(1) / sor(0), 2.5, tolerance = 1e-6)
  # Pairwise tables and marginals are still matched.
  j0 <- fit_joint(t2_pop())
  J <- crossprod(x * p, x)
  J0 <- crossprod(j0$cells * j0$prob, j0$cells)
  expect_equal(J, J0, tolerance = 1e-9)
  # Without interactions the additive result does not depend on it.
  m <- cm_model(pop, cm_impacts(t2_ids, c(2.5, 5, 7.5)))
  res <- deconflate(m, n_draws = 0)
  expect_equal(res$method, "full model")
  expect_equal(res$adjusted$adjusted, deconflate(example_supplement())$adjusted$adjusted, tolerance = 1e-8)
  # Among "no d3" the ratio is inverted.
  pn <- t2_pop(cm_associations("d1", "d2", 1 / 2.5, measure = "OR_ratio", among = "no d3"))
  expect_equal(fit_joint(pn)$prob, p, tolerance = 1e-9)
  # The three-way interface of earlier versions gives a message.
  expect_error(cm_population(t2_pop()$diseases, three_way = data.frame()), "OR_ratio",
               class = "deconflate_unsupported")
})

test_that("metadata is kept and impacts are put in disease order", {
  imp <- cm_impacts(c("d3", "d1", "d2"), c(7.5, 2.5, 5), label = "yield", units = "%",
                    note = "test")
  m <- cm_model(t2_pop(), imp)
  # Rows keep their order and numbers; results are in disease order.
  expect_equal(m$impacts$disease, c("d3", "d1", "d2"))
  expect_equal(m$impacts$row, 1:3)
  expect_equal(attr(m$impacts, "label"), "yield")
  expect_equal(attr(m$impacts, "units"), "%")
  expect_equal(attr(m$impacts, "kind"), "additive")
  res <- deconflate(m)
  expect_equal(res$label, "yield")
  expect_equal(res$units, "%")
  expect_output(print(res), "yield")
  expect_output(print(m), "cm_model")
  expect_output(print(m), "additive")
  expect_equal(res$adjusted$disease, t2_ids)
  # Impacts replaced after cm_model() are used as they are; cm_model()
  # checks them.
  m2 <- example_supplement()
  m2$impacts <- imp
  expect_equal(deconflate(m2)$adjusted$adjusted, res$adjusted$adjusted)
  expect_error(cm_model(m2, cm_impacts(c("d1", "d2"), c(2.5, 5))), "No impact for")
})

test_that("one population serves several outcomes, one impact table each", {
  pop <- t2_pop()
  yield <- cm_model(pop, cm_impacts(t2_ids, c(2.5, 5, 7.5), label = "yield", units = "%"))
  welfare <- cm_model(pop, cm_impacts(t2_ids, c(1, 0, 3), label = "welfare", units = "score"))
  expect_identical(yield$diseases, welfare$diseases)
  expect_identical(yield$associations, welfare$associations)
  ry <- deconflate(yield)
  rw <- deconflate(welfare)
  expect_equal(ry$adjusted$adjusted, deconflate(example_supplement())$adjusted$adjusted)
  expect_equal(rw$units, "score")
  expect_equal(rw$label, "welfare")
  # One joint fit can serve both.
  j <- fit_joint(pop)
  expect_equal(deconflate(welfare, joint = j)$adjusted$adjusted, rw$adjusted$adjusted, tolerance = 1e-8)
  expect_identical(deconflate(yield, joint = j)$joint, j)
})

test_that("disease ids are validated", {
  expect_error(cm_diseases(c("a:b", "c"), c(0.1, 0.2)), "must not contain")
  expect_error(cm_diseases(c("a;b", "c"), c(0.1, 0.2)), "must not contain")
  expect_error(cm_diseases(c("all", "c"), c(0.1, 0.2)), "reserved")
  expect_error(cm_diseases(c("others", "c"), c(0.1, 0.2)), "reserved")
  expect_error(cm_diseases(c("no SCK", "c"), c(0.1, 0.2)), "start with 'no '")
  expect_error(cm_diseases(c("a", "c"), c(0.1, 0.2), type = "prevalence"), "now called `measure`")
  expect_error(cm_diseases(c("a", "a"), c(0.1, 0.2)), "unique")
  expect_error(cm_diseases(c("a", "b"), c(0.1, 1.2)), class = "deconflate_infeasible")
})

test_that("associations adjusted for factors outside the model are used as crude", {
  d <- cm_diseases(c("a", "b"), c(0.2, 0.3))
  a <- cm_associations("a", "b", 2, adjusted_for = "parity")
  pop <- cm_population(d, a)
  expect_equal(pair_tables(pop)$status, "specified")
  expect_equal(pair_tables(pop)$or, 2)
  m <- cm_model(pop, cm_impacts(c("a", "b"), c(1, 2)))
  expect_match(deconflate(m, n_draws = 0)$notes, "Association row 1 is adjusted only for factors outside the model",
               all = FALSE)
  expect_error(cm_population(d, a, adjusted_associations = "use_as_marginal"), "no longer used",
               class = "deconflate_unsupported")
})

test_that("retired association measures are rejected with what to do instead", {
  for (meas in c("independent", "unknown", "table", "cond_prob")) {
    expect_error(cm_associations("a", "b", 1, measure = meas), class = "deconflate_unsupported")
  }
  # An unrelated pair is an odds ratio of 1; an unknown pair has no row.
  pop <- cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                       cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)))
  expect_equal(pair_tables(pop)$status, c("specified", "unknown", "specified"))
  expect_equal(pair_tables(t2_pop())$status, rep("specified", 3))
})
