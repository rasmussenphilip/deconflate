# fit_joint(): warm starts of the exact backend, and the sampled backend
# (Monte Carlo moment matching, sampling and raking) compared with exact
# enumeration.

t2_ids <- c("d1", "d2", "d3")
t2_pop <- function(three_way = NULL) {
  cm_population(cm_diseases(t2_ids, c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)),
                three_way = three_way)
}
t2_bad_pop <- function() {
  cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
}

sampled <- function(model, ...) {
  fit_joint(model, backend = "sampled", n_samples = 20000, n_chains = 500, seed = 11, ...)
}

test_that("a warm start gives the same fit", {
  j0 <- fit_joint(t2_pop())
  # A small change of one association: the start is close to the new fit.
  pop1 <- set_association(t2_pop(), "d1", "d2", 2.1)
  cold1 <- fit_joint(pop1)
  warm1 <- fit_joint(pop1, start = j0)
  expect_equal(warm1$prob, cold1$prob, tolerance = 1e-8)
  expect_lte(warm1$iterations, cold1$iterations)
  # Probabilities and associations changed: the same fit.
  pop2 <- set_association(t2_pop(), "d1", "d2", 2.5)
  pop2$diseases$value[3] <- pop2$diseases$prob[3] <- 0.22
  cold <- fit_joint(pop2)
  warm <- fit_joint(pop2, start = j0)
  expect_true(warm$converged)
  expect_equal(warm$prob, cold$prob, tolerance = 1e-8)
  expect_lt(warm$max_residual, 1e-10)
  expect_equal(warm$targets, cold$targets)
  # A start of another family is not used: different constrained pairs, a
  # different three-way term, other diseases, or the sampled backend.
  pop3 <- pop2
  pop3$associations <- pop3$associations[-2, ]
  expect_identical(fit_joint(pop3, start = j0)$prob, fit_joint(pop3)$prob)
  pop4 <- set_three_way(pop2, "d1", "d2", "d3", 2)
  expect_identical(fit_joint(pop4, start = j0)$prob, fit_joint(pop4)$prob)
  expect_false(warm_start_ok(j0, c("d1", "d2", "x"), pair_tables(pop2), NULL))
  js0 <- j0
  js0$backend <- "sampled"
  expect_false(warm_start_ok(js0, t2_ids, pair_tables(pop2), NULL))
  expect_false(warm_start_ok(NULL, t2_ids, pair_tables(pop2), NULL))
  expect_true(warm_start_ok(j0, t2_ids, pair_tables(pop2), NULL))
  # A start with the same three-way term is used.
  j4 <- fit_joint(pop4)
  pop5 <- set_association(pop4, "d2", "d3", 2.5)
  expect_equal(fit_joint(pop5, start = j4)$prob, fit_joint(pop5)$prob, tolerance = 1e-8)
})

test_that("the sampled backend matches the targets and the exact distribution (3 diseases)", {
  pop <- t2_pop()
  js <- sampled(pop)
  je <- fit_joint(pop)
  expect_s3_class(js, "cm_joint")
  expect_equal(js$backend, "sampled")
  expect_equal(je$backend, "exact")
  expect_true(js$converged)
  # Raking matches marginals and pairwise tables exactly on the sample.
  expect_lt(js$max_residual, 1e-9)
  expect_equal(unname(colSums(js$cells * js$prob)), pop$diseases$prob, tolerance = 1e-9)
  # The probability of every combination is close to the exact one.
  key <- function(x) apply(x, 1, paste, collapse = "")
  pe <- stats::setNames(je$prob, key(je$cells))
  ps <- stats::setNames(js$prob, key(js$cells))
  expect_lt(max(abs(pe[names(ps)] - ps)), 0.01)
  d <- js$diagnostics
  expect_named(d$summary, c("n_samples", "n_chains", "n_unique", "fit_iterations", "fit_converged",
                            "residual_sample", "residual_calibrated", "max_rhat", "min_ess",
                            "max_mcse"))
  expect_true(all(c("moment", "target", "sampled", "residual", "mcse", "rhat", "ess") %in%
                    names(d$moments)))
  expect_equal(nrow(d$moments), 3 + 3)
  expect_lt(d$summary$max_rhat, 1.1)
  expect_output(print(js), "sampled backend")
})

test_that("the global method gives the same answers with a sampled joint", {
  pop <- t2_pop()
  js <- sampled(pop)
  # Without interactions the result depends on the (raked, exact) pairs only.
  m <- cm_model(pop, cm_impacts(t2_ids, c(2.5, 5, 7.5)))
  expect_equal(deconflate(m, method = "global", joint = js)$adjusted$adjusted,
               deconflate(m)$adjusted$adjusted, tolerance = 1e-7)
  # With interactions it depends on triples: close to the exact result.
  int <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(0.01, 0.015))
  raw <- c(0.027146191901, 0.061369036759, 0.071166773469)
  mi <- cm_model(pop, cm_impacts(t2_ids, raw), int)
  rs <- deconflate(mi, joint = js)
  expect_equal(rs$method, "global")
  expect_equal(rs$adjusted$adjusted, c(0.02, 0.04, 0.06), tolerance = 0.02)
  # The event model works with a sampled joint too.
  ev <- cm_model(pop, cm_impacts(t2_ids, c(1.5, 2.0, 1.3), measure = "HR", estimand = "snapshot_crude"))
  re <- deconflate(ev, event_model = TRUE, overall_risk = 0.25, joint = js, n_draws = 0)
  expect_equal(re$adjusted$adjusted, c(1.3831529713, 1.8909315489, 1.1439870467), tolerance = 0.01)
  expect_equal(sum(re$attributable$by_disease$attributable), re$attributable$summary$attributable,
               tolerance = 1e-10)
})

test_that("without raking the sample residual is reported", {
  js <- sampled(t2_pop(), calibrate = FALSE)
  expect_true(is.na(js$diagnostics$summary$residual_calibrated))
  expect_equal(js$max_residual, js$diagnostics$summary$residual_sample)
  expect_gt(js$max_residual, 0)
  expect_lt(js$max_residual, 0.02)
})

test_that("a sampled joint is validated against the model like an exact one", {
  js <- sampled(t2_pop())
  m2 <- set_association(example_supplement(), "d1", "d2", 3)
  expect_error(deconflate(m2, method = "global", joint = js), class = "deconflate_error")
})

test_that("an unresolved sampled fit is reported, not called infeasible", {
  bad <- t2_bad_pop()
  expect_warning(js <- fit_joint(bad, backend = "sampled", n_samples = 2000, n_chains = 100,
                                 fit_iter = 20, seed = 1),
                 class = "deconflate_nonconvergence")
  expect_false(js$converged)
  m <- cm_model(bad, cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(deconflate(m, method = "global", joint = js), class = "deconflate_nonconvergence")
})

test_that("raking stops at once when no sampled combination can carry a target", {
  pop <- t2_pop()
  pt <- pair_tables(pop)
  cells <- disease_cells(3)
  cells <- cells[!(cells[, 1] == 1 & cells[, 2] == 1), , drop = FALSE]
  fit <- ipf_cells(cells, rep(1 / nrow(cells), nrow(cells)), pop$diseases$prob, pt, t2_ids,
                   1e-10, 10000L)
  expect_false(fit$converged)
  expect_equal(fit$iterations, 0L)
  expect_gt(fit$residual, 0.01)
})

test_that("more than 20 diseases use the sampled backend automatically", {
  n <- 24
  ids <- sprintf("x%02d", seq_len(n))
  pop <- cm_population(cm_diseases(ids, seq(0.05, 0.3, length.out = n)),
                       cm_associations(ids[-n], ids[-1], 2))
  m <- cm_model(pop, cm_impacts(ids, seq(1, 5, length.out = n)))
  expect_error(fit_joint(pop), "backend")
  # Unknown pairs need the global method, here with the sampled backend.
  plan <- plan_method(m, "auto", FALSE)
  expect_equal(plan$method, "global")
  expect_equal(plan$fit_args$backend, "sampled")
  expect_match(plan$notes, "sampled backend", all = FALSE)
  expect_match(plan$notes, "253 pairs without an association", all = FALSE)
  # An explicit backend is kept.
  expect_null(plan_method(m, "auto", FALSE, dots = list(backend = "exact"))$fit_args$backend)
  # With every pair given, the simultaneous method needs no joint at all.
  mi <- cm_model(with_independent_pairs(pop), m$impacts)
  pl <- plan_method(mi, "auto", FALSE)
  expect_equal(pl$method, "simultaneous")
  expect_length(pl$notes, 0)
})

test_that("the sampled backend handles more diseases than enumeration allows", {
  skip_on_cran()
  n <- 24
  ids <- sprintf("x%02d", seq_len(n))
  pop <- with_independent_pairs(cm_population(
    cm_diseases(ids, seq(0.05, 0.3, length.out = n)),
    cm_associations(ids[-n], ids[-1], 2)))
  js <- fit_joint(pop, backend = "sampled", n_samples = 20000, n_chains = 500, seed = 3)
  expect_true(js$converged)
  m <- cm_model(pop, cm_impacts(ids, seq(1, 5, length.out = n)))
  expect_equal(deconflate(m, method = "global", joint = js)$adjusted$adjusted,
               deconflate(m)$adjusted$adjusted, tolerance = 1e-6)
})

test_that("sampled and exact backends agree for the 12-disease global dairy population", {
  skip_on_cran()
  pop <- as_population(example_global_dairy(inputs = "tables"))
  je <- fit_joint(pop)
  js <- fit_joint(pop, backend = "sampled", n_samples = 50000, n_chains = 1000, seed = 7)
  expect_true(js$converged)
  # Distribution of the number of diseases per cow.
  nd <- function(j) vapply(0:5, function(k) sum(j$prob[rowSums(j$cells) == k]), numeric(1))
  expect_lt(max(abs(nd(je) - nd(js))), 0.01)
  # Global yield adjustment with one interaction.
  y <- example_global_dairy("yield", inputs = "tables")
  y$interactions <- cm_interactions("CM", "SCM", 1)
  # The exact solution changes the sign of some global dairy yield impacts
  # (see vignette("reproducing-published")), so warnings are off.
  re <- deconflate(y, method = "global", joint = je, warn = FALSE, n_draws = 0)$adjusted$adjusted
  rs <- deconflate(y, method = "global", joint = js, warn = FALSE, n_draws = 0)$adjusted$adjusted
  expect_lt(max(abs(re - rs)), 0.05 * max(abs(re)))
})
