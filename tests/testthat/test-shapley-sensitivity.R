test_that("cell-wise Shapley equals the closed form for pairwise interactions", {
  m <- example_supplement()
  ints <- cm_interactions("d1", "d2", 0.01, outcome = "yield")
  truth <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)
  m2 <- cm_model(m$diseases, m$associations,
                 simulate_raw_impacts(m, truth, ints, outcome = "yield"), interactions = ints)
  res <- deconflate(m2, method = "global")
  closed <- attribute_burden(res)
  adj <- stats::setNames(res$adjusted$adjusted, res$adjusted$disease)
  sh <- shapley_by_cell(res$joint, loss_additive(adj, list(list(diseases = c("d1", "d2"), value = 0.01))))
  expect_equal(sh$shapley, closed$total, tolerance = 1e-7)
  expect_equal(sum(sh$shapley), attr(sh, "total"), tolerance = 1e-12)
})

test_that("higher-order additive terms are split equally", {
  j <- fit_joint(example_supplement())
  main <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)
  sh <- shapley_by_cell(j, loss_additive(main, list(list(diseases = c("d1", "d2", "d3"), value = 0.03))))
  p123 <- sum(j$prob[rowSums(j$cells) == 3])
  P <- c(0.10, 0.15, 0.20)
  expect_equal(sh$shapley, unname(main * P + 0.03 * p123 / 3), tolerance = 1e-8)
})

test_that("multiplicative losses match the Python reference", {
  j <- fit_joint(example_supplement())
  sh <- shapley_by_cell(j, loss_multiplicative(c(d1 = 0.02, d2 = 0.04, d3 = 0.06)))
  expect_equal(sh$shapley, c(0.0019783456, 0.0059222732, 0.0119200649), tolerance = 1e-7)
  expect_equal(attr(sh, "total"), 0.01982068370187362, tolerance = 1e-9)
})

test_that("association screening reports influential pairs", {
  sc <- screen_associations(example_supplement(), or_values = c(0.5, 2), multipliers = c(0.5, 2))
  expect_s3_class(sc, "cm_screen")
  expect_equal(nrow(sc), 6)
  base <- attr(sc, "baseline_total")
  expect_equal(base, 0.02109228876453694, tolerance = 1e-9)
  # d1:d3 is specified (OR = 1), so it is scaled: 1 x 2 = 2.
  up <- sc[sc$pair == "d1:d3" & sc$scenario == "OR x 2", ]
  expect_equal(up$total, 0.02012881081303463, tolerance = 1e-9)
  expect_lt(up$rel_change, 0)
  expect_output(print(sc), "cm_screen")
})

test_that("interaction screening and one-at-a-time sensitivity run", {
  si <- screen_interactions(example_supplement(), outcome = "yield", values = c(-0.005, 0.005))
  expect_equal(nrow(si), 6)
  expect_true(all(is.finite(si$total)))
  oat <- sensitivity_oat(example_supplement(), variation = 0.2)
  expect_s3_class(oat, "cm_oat")
  expect_true(all(oat$swing >= 0))
  imp <- oat[grepl("^impact:", oat$input), ]
  expect_true(all(imp$total_high > imp$total_low))
})

test_that("scenario comparison and model editing", {
  base <- example_supplement()
  alt <- set_association(base, "d1", "d3", 2)
  expect_equal(alt$associations$value[pair_key(alt$associations$disease1, alt$associations$disease2) == "d1|d3"], 2)
  cs <- compare_scenarios(base = base, alt = alt)
  expect_equal(cs$totals$rel_to_first[1], 0)
  expect_equal(cs$totals$total, c(0.02109228876453694, 0.02012881081303463), tolerance = 1e-9)
  withint <- set_interaction(base, "d1", "d2", 0.01, "yield")
  expect_equal(nrow(withint$interactions), 1)
})
