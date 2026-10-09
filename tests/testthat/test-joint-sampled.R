# Sampled backend of fit_joint(): Monte Carlo moment matching, sampling and
# raking, compared with exact enumeration.

sampled <- function(model, ...) {
  fit_joint(model, backend = "sampled", n_samples = 20000, n_chains = 500, seed = 11, ...)
}

test_that("the sampled backend matches the targets and the exact distribution (3 diseases)", {
  pop <- supp_population()
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
  pop <- supp_population()
  js <- sampled(pop)
  # Without interactions the result depends on the (raked, exact) pairs only.
  m <- cm_model(pop, cm_impacts(ids3, c(2.5, 5, 7.5)))
  expect_equal(deconflate(m, method = "global", joint = js)$adjusted$adjusted,
               deconflate(m)$adjusted$adjusted, tolerance = 1e-7)
  # With interactions it depends on triples: close to the exact result.
  int <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(0.01, 0.015))
  raw <- c(0.027146191901, 0.061369036759, 0.071166773469)
  mi <- cm_model(pop, cm_impacts(ids3, raw), int)
  rs <- deconflate(mi, method = "global", joint = js)
  expect_equal(rs$adjusted$adjusted, c(0.02, 0.04, 0.06), tolerance = 0.02)
})

test_that("without raking the sample residual is reported", {
  js <- sampled(supp_population(), calibrate = FALSE)
  expect_true(is.na(js$diagnostics$summary$residual_calibrated))
  expect_equal(js$max_residual, js$diagnostics$summary$residual_sample)
  expect_gt(js$max_residual, 0)
  expect_lt(js$max_residual, 0.02)
})

test_that("a sampled joint is validated against the model like an exact one", {
  js <- sampled(supp_population())
  m2 <- set_association(example_supplement(), "d1", "d2", 3)
  expect_error(deconflate(m2, method = "global", joint = js), class = "deconflate_error")
})

test_that("an unresolved sampled fit is reported, not called infeasible", {
  bad <- bad_population()
  expect_warning(js <- fit_joint(bad, backend = "sampled", n_samples = 2000, n_chains = 100,
                                 fit_iter = 20, seed = 1),
                 class = "deconflate_nonconvergence")
  expect_false(js$converged)
  m <- cm_model(bad, cm_impacts(c("a", "b", "c"), c(1, 2, 3)))
  expect_error(deconflate(m, method = "global", joint = js), class = "deconflate_nonconvergence")
})

test_that("raking stops at once when no sampled combination can carry a target", {
  pop <- supp_population()
  pt <- pair_tables(pop)
  cells <- disease_cells(3)
  cells <- cells[!(cells[, 1] == 1 & cells[, 2] == 1), , drop = FALSE]
  fit <- ipf_cells(cells, rep(1 / nrow(cells), nrow(cells)), pop$diseases$prob, pt, ids3,
                   1e-10, 10000L)
  expect_false(fit$converged)
  expect_equal(fit$iterations, 0L)
  expect_gt(fit$residual, 0.01)
})

test_that("the sampled backend handles more diseases than enumeration allows", {
  skip_on_cran()
  n <- 24
  ids <- sprintf("x%02d", seq_len(n))
  pop <- cm_population(
    cm_diseases(ids, seq(0.05, 0.3, length.out = n)),
    cm_associations(ids[-n], ids[-1], 2))
  expect_error(fit_joint(pop), "backend")
  js <- fit_joint(pop, backend = "sampled", n_samples = 20000, n_chains = 500, seed = 3)
  expect_true(js$converged)
  m <- cm_model(pop, cm_impacts(ids, seq(1, 5, length.out = n)))
  expect_equal(deconflate(m, method = "global", joint = js)$adjusted$adjusted,
               deconflate(m)$adjusted$adjusted, tolerance = 1e-6)
})

test_that("sampled and exact backends agree for the 12-disease global dairy population", {
  skip_on_cran()
  pop <- example_global_dairy(inputs = "tables")$population
  je <- fit_joint(pop)
  js <- fit_joint(pop, backend = "sampled", n_samples = 50000, n_chains = 1000, seed = 7)
  expect_true(js$converged)
  # Distribution of the number of diseases per cow.
  nd <- function(j) vapply(0:5, function(k) sum(j$prob[rowSums(j$cells) == k]), numeric(1))
  expect_lt(max(abs(nd(je) - nd(js))), 0.01)
  # Global yield adjustment with one interaction.
  y <- example_global_dairy(inputs = "tables")$models$yield
  y$interactions <- cm_interactions("CM", "SCM", 1)
  # The exact solution changes the sign of some global dairy yield impacts
  # (see vignette("reproducing-published")), so warnings are off.
  re <- deconflate(y, method = "global", joint = je, warn = FALSE)$adjusted$adjusted
  rs <- deconflate(y, method = "global", joint = js, warn = FALSE)$adjusted$adjusted
  expect_lt(max(abs(re - rs)), 0.05 * max(abs(re)))
})
