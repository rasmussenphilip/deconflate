# Pairwise tables (R/pairwise.R) and association measures. Reference values
# computed independently in Python (python_reference.py; reference_v040.py,
# section 6).

t2_supp <- function() {
  cm_population(cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
                cm_associations(c("d1", "d1", "d2"), c("d2", "d3", "d3"), c(2, 1, 3)))
}

# a, b, c with OR a:b = 3 and b:c = 2; the pair a:c has no association.
t2_unknown_model <- function() {
  cm_model(cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25)),
           cm_impacts(c("a", "b", "c"), c(3, 4, 5)),
           associations = cm_associations(c("a", "b"), c("b", "c"), c(3, 2)))
}

test_that("odds ratio quadratic reproduces the Supplementary File", {
  p11 <- or_to_joint(2, 0.10, 0.15)
  expect_equal(p11 / 0.15, 0.1631959501785622, tolerance = 1e-10)  # P(1 | 2)
  expect_equal(p11 / 0.10, 0.24479392526784327, tolerance = 1e-10) # P(2 | 1)
  expect_equal(or_to_joint(1, 0.10, 0.15), 0.10 * 0.15)
  expect_error(or_to_joint(-1, 0.1, 0.15), class = "deconflate_infeasible")
})

test_that("odds ratio conversion round-trips, including OR < 1 and OR near 1", {
  grid <- expand.grid(or = c(0.2, 0.4, 0.999999, 1 + 1e-9, 1.5, 3, 9.7),
                      p1 = c(0.02, 0.3, 0.6), p2 = c(0.05, 0.38, 0.8))
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    p11 <- or_to_joint(g$or, g$p1, g$p2)
    expect_true(p11 >= max(0, g$p1 + g$p2 - 1) && p11 <= min(g$p1, g$p2))
    expect_equal(joint_to_or(p11, g$p1, g$p2), g$or, tolerance = 1e-6)
  }
})

test_that("all association measures map to the same 2x2 table", {
  p1 <- 0.12; p2 <- 0.3
  p11 <- or_to_joint(2.5, p1, p2)
  pd1_d2 <- p11 / p2
  pd1_nd2 <- (p1 - p11) / (1 - p2)
  phi <- (p11 - p1 * p2) / sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  expect_equal(association_to_joint("OR", 2.5, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("RR", pd1_d2 / pd1_nd2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("RD", pd1_d2 - pd1_nd2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("cond_prob", pd1_d2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("phi", phi, p1, p2), p11, tolerance = 1e-12)
  # Independence is an odds ratio of 1; the retired measures are not converted.
  expect_equal(association_to_joint("OR", 1, p1, p2), p1 * p2)
  expect_error(association_to_joint("independent", NA, p1, p2), "Cannot convert")
  expect_error(association_to_joint("table", 2, p1, p2), "Cannot convert")
})

test_that("incompatible measures are flagged as infeasible", {
  expect_error(association_to_joint("cond_prob", 0.9, 0.1, 0.5),
               class = "deconflate_infeasible")
  expect_error(association_to_joint("RD", 0.9, 0.1, 0.5),
               class = "deconflate_infeasible")
  expect_error(cm_associations("a", "b", 0, measure = "OR"), class = "deconflate_infeasible")
  expect_error(cm_associations("a", "b", -2, measure = "RR"), class = "deconflate_infeasible")
})

test_that("retired measures and count columns are rejected with what to do instead", {
  expect_error(cm_associations("a", "b", 1, measure = "independent"), class = "deconflate_unsupported")
  expect_error(cm_associations("a", "b", 1, measure = "independent"), "odds ratio of 1")
  expect_error(cm_associations("a", "b", 1, measure = "unknown"), class = "deconflate_unsupported")
  expect_error(cm_associations("a", "b", 1, measure = "unknown"), "leave the pair out")
  expect_error(cm_associations("a", "b", 2, measure = "table"), class = "deconflate_unsupported")
  expect_error(cm_associations("a", "b", 2, measure = "table"), "n11 \\* n00 / \\(n10 \\* n01\\)")
  # A retired measure in any row is reported.
  expect_error(cm_associations(c("a", "a"), c("b", "c"), c(2, 1), measure = c("OR", "independent")),
               class = "deconflate_unsupported")
  # Counts are entered as their odds ratio; there are no count arguments.
  expect_error(cm_associations("a", "b", measure = "OR", n11 = 20, n10 = 30, n01 = 40, n00 = 110),
               "unused argument")
  a <- cm_associations("a", "b", (20 * 110) / (30 * 40))
  expect_equal(a$measure, "OR")
  expect_equal(a$value, 20 * 110 / (30 * 40))
  expect_false("corrected" %in% names(a))
  # A value is required, and cannot be missing.
  expect_error(cm_associations("a", "b"), "value")
  expect_error(cm_associations(c("a", "a"), c("b", "c"), c(2, NA)), "missing")
  expect_error(cm_associations("a", "b", Inf), "finite")
  expect_error(cm_associations(c("a", "b"), c("b", "a"), c(2, 3)), "Duplicate")
  expect_error(cm_associations("a", "a", 2), "two different diseases")
})

test_that("excess matrix matches the Supplementary File", {
  E <- excess_matrix(example_supplement())
  # E[k, i] = P(k | i) - P(k | not i)
  expect_equal(E["d2", "d1"], 0.10532658, tolerance = 1e-7)  # ep_2,1 = 0.11
  expect_equal(E["d1", "d2"], 0.07434818, tolerance = 1e-7)  # ep_1,2 = 0.07
  expect_equal(E["d3", "d2"], 0.20962356, tolerance = 1e-7)  # ep_3,2 = 0.21
  expect_equal(E["d2", "d3"], 0.16704377, tolerance = 1e-7)  # ep_2,3 = 0.17
  expect_equal(E["d1", "d3"], 0, tolerance = 1e-12)
})

test_that("pair tables report specified and unknown pairs", {
  pt <- pair_tables(t2_supp())
  expect_equal(names(pt), c("disease1", "disease2", "measure", "status", "p11", "or",
                            "ep_2_given_1", "ep_1_given_2"))
  expect_equal(pt$status, rep("specified", 3))
  expect_equal(pt$measure, rep("OR", 3))
  expect_equal(pt$or, c(2, 1, 3), tolerance = 1e-10)
  pu <- pair_tables(t2_unknown_model())
  expect_equal(paste(pu$disease1, pu$disease2), c("a b", "a c", "b c"))
  expect_equal(pu$status, c("specified", "unknown", "specified"))
  expect_true(is.na(pu$measure[2]))
  expect_true(is.na(pu$p11[2]))
  expect_true(is.na(pu$or[2]))
  expect_true(all(is.finite(pu$p11[-2])))
  # A pair given an odds ratio of 1 is specified (independent), not unknown.
  p1 <- pair_tables(set_association(t2_unknown_model(), "a", "c", 1))
  expect_equal(p1$status, rep("specified", 3))
  expect_equal(p1$p11[2], 0.2 * 0.25)
  # Without associations every pair is unknown.
  p0 <- pair_tables(cm_population(cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))))
  expect_equal(p0$status, rep("unknown", 3))
  expect_true(all(is.na(p0$measure)))
  # Directional measures keep their orientation.
  pr <- pair_tables(cm_population(cm_diseases(c("a", "b"), c(0.2, 0.3)),
                                  cm_associations("b", "a", 0.5, measure = "cond_prob")))
  expect_equal(pr$p11, 0.5 * 0.2)
  expect_equal(nrow(pair_tables(cm_population(cm_diseases("a", 0.2)))), 0)
})

test_that("unknown pairs are rejected by pairwise methods and filled in by the global method", {
  m <- t2_unknown_model()
  expect_equal(global_reasons(m), "1 pair without an association (unknown)")
  expect_error(adjust_impacts(m, method = "simultaneous"), class = "deconflate_unknown_pairs")
  expect_error(adjust_impacts(m, method = "simultaneous"), "odds ratio of 1")
  expect_error(adjust_impacts(m, method = "published"), class = "deconflate_unknown_pairs")
  expect_error(excess_matrix(m), class = "deconflate_unknown_pairs")
  # deconflate() uses the global method, with a note giving the reason.
  res <- deconflate(m, n_draws = 0)
  expect_s3_class(res, "cm_result")
  expect_equal(res$method, "global")
  expect_equal(res$notes, "The global method was used because of 1 pair without an association (unknown).")
  # method = "simultaneous" behaves the same, and says so.
  rs <- deconflate(m, method = "simultaneous", n_draws = 0)
  expect_equal(rs$method, "global")
  expect_match(rs$notes, "instead of the simultaneous method")
  expect_equal(rs$adjusted$adjusted, res$adjusted$adjusted)
  expect_length(deconflate(m, method = "global", n_draws = 0)$notes, 0)
  # Python reference (reference_v040.py, section 6).
  expect_equal(res$adjusted$adjusted, c(2.1, 2.9872624936, 4.4812715344), tolerance = 1e-8)
  expect_equal(res$totals$adjusted_total, 2.4364966316833225, tolerance = 1e-9)
  expect_equal(res$unknown_pairs$disease1, "a")
  expect_equal(res$unknown_pairs$disease2, "c")
  expect_equal(res$unknown_pairs$fitted_or, 1.1945293054612587, tolerance = 1e-7)
  expect_output(print(res), "Unknown pairs: 1")
  # Unconstrained a:c is implied by the maximum-entropy fit, so it differs
  # from stating independence (an odds ratio of 1).
  ind <- deconflate(set_association(m, "a", "c", 1), n_draws = 0)
  expect_equal(ind$method, "simultaneous")
  expect_equal(nrow(ind$unknown_pairs), 0)
  expect_equal(ind$adjusted$adjusted, c(2.2632638445, 2.9469446222, 4.5487022146), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(res$adjusted$adjusted, ind$adjusted$adjusted)))
  # compare_methods() runs the methods as asked: the pairwise ones fail.
  cmp <- compare_methods(m)
  expect_equal(cmp$methods, "global")
  expect_setequal(names(cmp$failed), c("published", "simultaneous"))
  expect_match(cmp$failed[["simultaneous"]], "No association for a-c")
})

test_that("deconflate() needs at least one association", {
  m <- cm_model(cm_diseases(c("a", "b"), c(0.2, 0.3)), cm_impacts(c("a", "b"), c(1, 2)))
  expect_null(m$associations)
  expect_error(deconflate(m), class = "deconflate_unsupported")
  expect_error(deconflate(m), "screen_associations")
})
