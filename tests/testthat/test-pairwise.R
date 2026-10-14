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

test_that("odds ratio quadratic matches the Supplementary File", {
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
  expect_equal(association_to_joint("prevalence", pd1_d2, p1, p2), p11, tolerance = 1e-12)
  expect_equal(association_to_joint("phi", phi, p1, p2), p11, tolerance = 1e-12)
  # Independence is an odds ratio of 1; the retired measures are not converted.
  expect_equal(association_to_joint("OR", 1, p1, p2), p1 * p2)
  expect_error(association_to_joint("independent", NA, p1, p2), "Cannot convert")
  expect_error(association_to_joint("table", 2, p1, p2), "Cannot convert")
  expect_error(association_to_joint("cond_prob", pd1_d2, p1, p2), "Cannot convert")
})

test_that("incompatible measures are flagged as infeasible", {
  expect_error(association_to_joint("prevalence", 0.9, 0.1, 0.5),
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
               "Unknown arguments: n11")
  expect_error(cm_associations("a", "b", 0.3, measure = "cond_prob"), "now 'prevalence'",
               class = "deconflate_unsupported")
  expect_error(cm_associations("a", "b", 2, adjusted = TRUE), "adjusted_for", class = "deconflate_unsupported")
  expect_error(cm_associations(disease1 = "a", disease2 = "b", value = 2), "associated_with",
               class = "deconflate_unsupported")
  a <- cm_associations("a", "b", (20 * 110) / (30 * 40))
  expect_equal(a$measure, "OR")
  expect_equal(a$value, 20 * 110 / (30 * 40))
  expect_false("corrected" %in% names(a))
  # A value is required, and cannot be missing.
  expect_error(cm_associations("a", "b"), "value")
  expect_error(cm_associations(c("a", "a"), c("b", "c"), c(2, NA)), "missing")
  expect_error(cm_associations("a", "b", Inf), "finite")
  # Several rows about one pair are allowed (they are pooled); a row about
  # one disease on both sides is not.
  expect_s3_class(cm_associations(c("a", "b"), c("b", "a"), c(2, 3)), "cm_associations")
  expect_error(cm_population(cm_diseases("a", 0.2), cm_associations("a", "a", 2)), "share a")
  expect_error(cm_associations("a", "b", 1.2, measure = "prevalence"), "between 0 and 1")
  expect_error(cm_associations("a", "b", 0, measure = "OR_ratio"), "positive")
})

test_that("excess probabilities match the Supplementary File", {
  pt <- pair_tables(example_supplement())
  rownames(pt) <- paste(pt$disease1, pt$disease2, sep = ":")
  # ep_2_given_1 = P(disease2 | disease1) - P(disease2 | not disease1)
  expect_equal(pt["d1:d2", "ep_2_given_1"], 0.10532658, tolerance = 1e-7)  # ep_2,1 = 0.11
  expect_equal(pt["d1:d2", "ep_1_given_2"], 0.07434818, tolerance = 1e-7)  # ep_1,2 = 0.07
  expect_equal(pt["d2:d3", "ep_2_given_1"], 0.20962356, tolerance = 1e-7)  # ep_3,2 = 0.21
  expect_equal(pt["d2:d3", "ep_1_given_2"], 0.16704377, tolerance = 1e-7)  # ep_2,3 = 0.17
  expect_equal(pt["d1:d3", "ep_1_given_2"], 0, tolerance = 1e-12)
})

test_that("pair tables report specified and unknown pairs", {
  pt <- pair_tables(t2_supp())
  expect_equal(names(pt), c("disease1", "disease2", "measure", "status", "rows", "p11", "or",
                            "ep_2_given_1", "ep_1_given_2"))
  expect_equal(pt$rows, c("1", "2", "3"))
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
                                  cm_associations("b", "a", 0.5, measure = "prevalence")))
  expect_equal(pr$p11, 0.5 * 0.2)
  expect_equal(nrow(pair_tables(cm_population(cm_diseases("a", 0.2)))), 0)
  # A pair that only other kinds of rows mention.
  po <- pair_tables(cm_population(cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25)),
                                  cm_associations("a", "b", 2, among = "c")))
  expect_equal(po$status, c("other rows", "unknown", "unknown"))
  expect_true(is.na(po$p11[1]))
})

test_that("unknown pairs rule out the pairwise shortcut and are filled in by the model", {
  m <- t2_unknown_model()
  expect_false(shortcut_ok(m))
  # deconflate() fits the full model, which fills the pair in.
  res <- deconflate(m, n_draws = 0)
  expect_s3_class(res, "cm_result")
  expect_equal(res$method, "full model")
  expect_length(res$notes, 0)
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
  expect_equal(ind$method, "pairwise shortcut")
  expect_equal(nrow(ind$unknown_pairs), 0)
  expect_equal(ind$adjusted$adjusted, c(2.2632638445, 2.9469446222, 4.5487022146), tolerance = 1e-8)
  expect_false(isTRUE(all.equal(res$adjusted$adjusted, ind$adjusted$adjusted)))
})

test_that("deconflate() needs at least one association", {
  m <- cm_model(cm_diseases(c("a", "b"), c(0.2, 0.3)), cm_impacts(c("a", "b"), c(1, 2)))
  expect_null(m$associations)
  expect_error(deconflate(m), class = "deconflate_unsupported")
  expect_error(deconflate(m), "screen_associations")
})
