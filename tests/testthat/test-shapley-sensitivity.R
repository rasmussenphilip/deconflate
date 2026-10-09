# Shapley attribution (R/shapley.R) and sensitivity / scenario tools
# (R/sensitivity.R). Supplement impacts are in percent (2.5, 5, 7.5);
# reference values from inst/validation/reference_v020_tests.py and
# reference_2022_shapley.py (multiplicative Shapley values).

supp_truth <- c(d1 = 2, d2 = 4, d3 = 6)
supp_base_total <- 2.109228876453694       # simultaneous aggregate, percent
supp_int_total <- 2.0890878497857703       # global, interaction d1:d2 = 1
supp_int_tw_total <- c(`0.5` = 2.0913976662985534, `2` = 2.086745618267496)

test_that("additive contributions add up to the aggregate", {
  m <- example_supplement()
  res <- deconflate(m, method = "global")
  expect_equal(sum(res$contributions$total), res$totals$adjusted_total)
  expect_equal(res$totals$adjusted_total, supp_base_total, tolerance = 1e-8)
  expect_equal(attribute_burden(res), res$contributions)
  expect_error(attribute_burden(m), "deconflate")
  res_i <- deconflate(set_interaction(m, "d1", "d2", 1), method = "global")
  expect_equal(sum(res_i$contributions$total), res_i$totals$adjusted_total)
  expect_equal(res_i$totals$adjusted_total, supp_int_total, tolerance = 1e-8)
  P12 <- res_i$joint_pairs["d1", "d2"]
  expect_equal(res_i$totals$interaction_total, P12)
  expect_equal(res_i$contributions$interaction, c(0.5 * P12, 0.5 * P12, 0))
  expect_equal(res_i$totals$adjusted_total,
               sum(c(0.10, 0.15, 0.20) * res_i$adjusted$adjusted) + P12, tolerance = 1e-12)
})

test_that("cell-wise Shapley equals the closed form for pairwise interactions", {
  pop <- supp_population()
  ints <- cm_interactions("d1", "d2", 1)
  m2 <- cm_model(pop, simulate_raw_impacts(pop, supp_truth, ints), ints)
  res <- deconflate(m2, method = "global")
  expect_equal(res$adjusted$adjusted, unname(supp_truth), tolerance = 1e-8)
  closed <- attribute_burden(res)
  adj <- stats::setNames(res$adjusted$adjusted, res$adjusted$disease)
  expect_silent(sh <- shapley_by_cell(res$joint,
                                      loss_additive(adj, list(list(diseases = c("d1", "d2"), value = 1)))))
  expect_equal(sh$disease, ids3)
  expect_equal(sh$shapley, closed$total, tolerance = 1e-7)
  expect_equal(sum(sh$shapley), attr(sh, "total"), tolerance = 1e-12)
  expect_equal(attr(sh, "total"), res$totals$adjusted_total, tolerance = 1e-8)
  expect_equal(attr(sh, "skipped_mass"), 0)
  expect_equal(attr(sh, "skipped_loss"), 0)
  expect_equal(sum(sh$share), 1)
})

test_that("higher-order additive terms are split equally", {
  j <- fit_joint(supp_population())
  main <- c(d1 = 2, d2 = 4, d3 = 6)
  sh <- shapley_by_cell(j, loss_additive(main, list(list(diseases = c("d1", "d2", "d3"), value = 3))))
  p123 <- sum(j$prob[rowSums(j$cells) == 3])
  expect_equal(sh$shapley, unname(main * c(0.10, 0.15, 0.20) + 3 * p123 / 3), tolerance = 1e-8)
  expect_error(loss_additive(c(2, 4)), "named")
  expect_error(loss_additive(main, list(list(diseases = "d1", value = 1))), "two or more")
  expect_error(loss_multiplicative(c(0.1, 0.2)), "named")
  expect_error(shapley_by_cell(supp_population(), loss_additive(main)), "fit_joint")
})

test_that("multiplicative losses match the Python reference", {
  j <- fit_joint(supp_population())
  sh <- shapley_by_cell(j, loss_multiplicative(c(d1 = 0.02, d2 = 0.04, d3 = 0.06)))
  expect_equal(sh$shapley, c(0.0019783456, 0.0059222732, 0.0119200649), tolerance = 1e-7)
  expect_equal(attr(sh, "total"), 0.01982068370187362, tolerance = 1e-9)
  expect_equal(attr(sh, "allocated"), sum(sh$shapley))
})

test_that("skipped combinations are reported, not silently dropped", {
  j <- fit_joint(supp_population())
  f <- loss_multiplicative(c(d1 = 0.02, d2 = 0.04, d3 = 0.06))
  expect_warning(sh <- shapley_by_cell(j, f, max_present = 2), class = "deconflate_incomplete_allocation")
  p123 <- j$prob[rowSums(j$cells) == 3]
  expect_equal(attr(sh, "skipped_mass"), p123)
  expect_equal(attr(sh, "skipped_loss"), p123 * f(c(d1 = 1, d2 = 1, d3 = 1)))
  expect_gt(attr(sh, "skipped_loss"), 0)
  expect_equal(attr(sh, "allocated"), sum(sh$shapley))
  expect_equal(attr(sh, "allocated") + attr(sh, "skipped_loss"), attr(sh, "total"), tolerance = 1e-12)
  # The total still covers every combination.
  expect_equal(attr(sh, "total"), 0.01982068370187362, tolerance = 1e-9)
  expect_equal(sum(sh$share), 1)
})

test_that("association screening reports influential pairs", {
  sc <- screen_associations(example_supplement(), or_values = c(0.5, 2), multipliers = c(0.5, 2))
  expect_s3_class(sc, "cm_screen")
  expect_equal(nrow(sc), 6)
  expect_true(all(c("pair", "status", "scenario", "total", "change", "rel_change", "max_rank_shift",
                    "rank_corr", "failed") %in% names(sc)))
  expect_true(all(is.na(sc$failed)))
  expect_equal(attr(sc, "baseline_total"), supp_base_total, tolerance = 1e-10)
  # d1:d3 is specified (OR = 1), so it is scaled: 1 x 2 = 2.
  up <- sc[sc$pair == "d1:d3" & sc$scenario == "OR x 2", ]
  expect_equal(up$total, 2.0128810813034628, tolerance = 1e-10)
  expect_lt(up$rel_change, 0)
  expect_equal(up$change, up$total - supp_base_total, tolerance = 1e-10)
  dn <- sc[sc$pair == "d2:d3" & sc$scenario == "OR x 0.5", ]
  expect_equal(dn$total, 2.306649129926649, tolerance = 1e-10)
  # Sorted by the absolute relative change.
  expect_equal(sc$pair[1], "d2:d3")
  expect_true(all(diff(abs(sc$rel_change)) <= 0))
  expect_output(print(sc), "cm_screen")
  # Restricting to pairs (in either order).
  sp <- screen_associations(example_supplement(), pairs = "d3:d1")
  expect_equal(unique(sp$pair), "d1:d3")
  expect_equal(nrow(sp), 2)
  expect_equal(attr(sp, "baseline_total"), supp_base_total, tolerance = 1e-10)
  expect_equal(sp$total[sp$scenario == "OR x 2"], 2.0128810813034628, tolerance = 1e-10)
  # The metric is always the adjusted aggregate: there is no valuation.
  expect_false("valuation" %in% names(formals(screen_associations)))
  expect_error(screen_associations(example_supplement(), pairs = "d1:d3",
                                   valuation = list(observed = 10000)),
               "unused argument")
})

test_that("association screening keeps failed scenarios with the reason", {
  m <- cm_model(cm_population(cm_diseases(ids3, c(0.10, 0.15, 0.20)),
                              cm_associations(c("d1", "d2"), c("d2", "d3"), c(0.3, 3),
                                              measure = c("cond_prob", "OR"))),
                cm_impacts(ids3, c(2.5, 5, 7.5)))
  sc <- screen_associations(m, multipliers = c(-1, 0.5, 4), pairs = c("d1:d2", "d2:d3"))
  expect_equal(nrow(sc), 6)
  key <- paste(sc$pair, sc$scenario)
  failed <- stats::setNames(sc$failed, key)
  # P(d1 | d2) = -0.3 or 1.2 is impossible; a negative odds ratio is invalid.
  expect_equal(unname(failed[c("d1:d2 cond_prob x -1", "d1:d2 cond_prob x 4", "d2:d3 OR x -1")]),
               c("infeasible", "infeasible", "invalid value"))
  ok <- c("d1:d2 cond_prob x 0.5", "d2:d3 OR x 0.5", "d2:d3 OR x 4")
  expect_true(all(is.na(failed[ok])))
  expect_true(all(is.finite(sc$total[key %in% ok])))
  expect_true(all(is.na(sc$total[!is.na(sc$failed)])))
  # Failed scenarios come last.
  expect_true(all(!is.na(sc$failed[4:6])))
})

test_that("screens report the adjusted aggregate, whatever its size", {
  # Impacts 45 times the supplement's (percent): an aggregate of 94.9%, and
  # 103.8% when the d2:d3 odds ratio is halved (reference_v020_tests.py
  # screen values, scaled). The aggregate is linear in the impacts and is
  # not converted into anything else, so every scenario is kept.
  m <- supp_model(45 * c(2.5, 5, 7.5))
  sc <- screen_associations(m)
  expect_equal(nrow(sc), 6)
  expect_true(all(is.na(sc$failed)))
  expect_true(all(is.finite(sc$total)))
  expect_equal(attr(sc, "baseline_total"), 45 * supp_base_total, tolerance = 1e-9)
  big <- sc$pair == "d2:d3" & sc$scenario == "OR x 0.5"
  expect_equal(sum(big), 1)
  expect_equal(sc$total[big], 45 * 2.306649129926649, tolerance = 1e-9)
  expect_gt(sc$total[big], 100)
  # The same scenario in a scenario comparison.
  cs <- compare_scenarios(base = m, weak = set_association(m, "d2", "d3", 1.5))
  expect_true(all(is.na(cs$totals$failed)))
  expect_equal(cs$totals$total, 45 * c(supp_base_total, 2.306649129926649), tolerance = 1e-9)
  expect_false("valuation" %in% names(formals(compare_scenarios)))
})

test_that("screens without scenarios return an empty screen", {
  two <- cm_model(cm_population(cm_diseases(c("a", "b"), c(0.1, 0.2)), cm_associations("a", "b", 2)),
                  cm_impacts(c("a", "b"), c(1, 2)))
  st <- screen_three_way(two)
  expect_s3_class(st, "cm_screen")
  expect_equal(nrow(st), 0)
  expect_true(all(c("pair", "total", "failed") %in% names(st)))
  expect_equal(nrow(screen_interactions(two, values = 1)), 1)
  expect_error(screen_interactions(two, values = 1, pairs = "a:c"), "Unknown or malformed pairs")
  expect_error(screen_associations(two, pairs = "a:a"), "Unknown or malformed pairs")
  oat <- sensitivity_oat(supp_model(c(0, 0, 0)), inputs = "impact")
  expect_s3_class(oat, "cm_oat")
  expect_equal(nrow(oat), 0)
})

test_that("interaction screening runs every pair and value", {
  si <- screen_interactions(example_supplement(), values = c(-0.5, 0.5))
  expect_s3_class(si, "cm_screen")
  expect_equal(nrow(si), 6)
  expect_true(all(is.finite(si$total)))
  expect_true(all(is.na(si$failed)))
  expect_equal(unique(si$status), "interaction")
  expect_equal(attr(si, "baseline_total"), supp_base_total, tolerance = 1e-8)
  s1 <- screen_interactions(example_supplement(), values = 1, pairs = "d1:d2")
  expect_equal(nrow(s1), 1)
  expect_equal(s1$scenario, "delta = 1")
  expect_equal(s1$total, supp_int_total, tolerance = 1e-8)
})

test_that("three-way terms change global results with interactions only", {
  base <- example_supplement()
  # Without interactions, the global result depends on the pairs alone.
  st <- screen_three_way(base)
  expect_s3_class(st, "cm_screen")
  expect_equal(nrow(st), 2)
  expect_equal(st$pair, rep("d1:d2:d3", 2))
  expect_equal(unique(st$status), "three-way")
  expect_true(all(is.na(st$failed)))
  expect_lt(max(abs(st$rel_change)), 1e-7)
  tw <- set_three_way(base, "d1", "d2", "d3", 2)
  expect_equal(deconflate(tw, method = "global")$adjusted$adjusted,
               deconflate(base, method = "simultaneous")$adjusted$adjusted, tolerance = 1e-7)
  # With an interaction, the probability of the triple matters.
  int <- set_interaction(base, "d1", "d2", 1)
  sti <- screen_three_way(int, ratios = c(0.5, 2), triples = list(c("d1", "d2", "d3")))
  expect_equal(attr(sti, "baseline_total"), supp_int_total, tolerance = 1e-8)
  expect_true(all(abs(sti$rel_change) > 1e-4))
  tot <- stats::setNames(sti$total, sti$scenario)
  expect_equal(unname(tot[c("ratio = 0.5", "ratio = 2")]), unname(supp_int_tw_total), tolerance = 1e-7)
  res_tw <- deconflate(set_three_way(int, "d1", "d2", "d3", 2), method = "global")
  expect_equal(res_tw$totals$adjusted_total, supp_int_tw_total[["2"]], tolerance = 1e-7)
})

test_that("set_three_way, set_interaction and set_association edit models", {
  base <- example_supplement()
  tw <- set_three_way(base, "d1", "d2", "d3", 2)
  expect_s3_class(tw$three_way, "cm_three_way")
  tw2 <- set_three_way(tw, "d3", "d1", "d2", 0.5)
  expect_equal(nrow(tw2$three_way), 1)
  expect_equal(tw2$three_way$ratio, 0.5)
  expect_error(set_three_way(base, "d1", "d2", "zz", 2), "Unknown diseases")
  expect_error(set_three_way(base, "d1", "d2", "d3", -1), "positive")
  pop_tw <- set_three_way(supp_population(), "d1", "d2", "d3", 2)
  expect_s3_class(pop_tw, "cm_population")

  wi <- set_interaction(base, "d1", "d2", 1)
  expect_equal(nrow(wi$interactions), 1)
  wi2 <- set_interaction(wi, "d2", "d1", 2)
  expect_equal(nrow(wi2$interactions), 1)
  expect_equal(wi2$interactions$value, 2)
  wi3 <- set_interaction(wi2, "d2", "d3", 0.5)
  expect_equal(nrow(wi3$interactions), 2)
  expect_error(set_interaction(supp_population(), "d1", "d2", 1), "cm_model")
  expect_error(set_interaction(base, "d1", "zz", 1), "Unknown diseases")
  expect_error(set_association(base, "zz", "d1", 2), "Unknown diseases")

  alt <- set_association(base, "d3", "d1", 2)
  a <- alt$associations
  expect_equal(nrow(a), 3)
  expect_equal(a$value[pair_key(a$disease1, a$disease2) == "d1|d3"], 2)
  new_pop <- set_association(cm_population(cm_diseases(ids3, c(0.1, 0.15, 0.2))), "d1", "d2", 2)
  expect_equal(nrow(new_pop$associations), 1)
})

test_that("one-at-a-time sensitivity varies every input", {
  oat <- sensitivity_oat(example_supplement(), variation = 0.2)
  expect_s3_class(oat, "cm_oat")
  expect_equal(names(oat), c("input", "value", "total_low", "total_high", "swing", "rel_swing"))
  expect_equal(nrow(oat), 9)
  expect_setequal(oat$input, c("prob:d1", "prob:d2", "prob:d3", "assoc:d1:d2", "assoc:d1:d3",
                               "assoc:d2:d3", "impact:d1", "impact:d2", "impact:d3"))
  expect_true(all(oat$swing >= 0))
  expect_true(all(diff(oat$swing) <= 0))
  expect_equal(attr(oat, "baseline_total"), supp_base_total, tolerance = 1e-10)
  imp <- oat[grepl("^impact:", oat$input), ]
  expect_true(all(imp$total_high > imp$total_low))
  # The aggregate is linear in the impacts: low and high are symmetric.
  expect_equal((imp$total_low + imp$total_high) / 2, rep(supp_base_total, 3), tolerance = 1e-10)
  ov <- sensitivity_oat(example_supplement(), inputs = "impact")
  expect_equal(nrow(ov), 3)
  expect_equal(attr(ov, "baseline_total"), supp_base_total, tolerance = 1e-10)
  expect_error(sensitivity_oat(example_supplement(), inputs = "impact",
                               valuation = list(observed = 10000)),
               "unused argument")
})

test_that("scenario comparison tabulates totals and rankings", {
  base <- example_supplement()
  alt <- set_association(base, "d1", "d3", 2)
  cs <- compare_scenarios(base = base, alt = alt)
  expect_equal(cs$totals$scenario, c("base", "alt"))
  expect_equal(cs$totals$rel_to_first[1], 0)
  expect_equal(cs$totals$total, c(supp_base_total, 2.0128810813034628), tolerance = 1e-10)
  expect_equal(nrow(cs$by_disease), 6)
  expect_equal(names(cs$by_disease), c("scenario", "disease", "contribution", "rank"))
  # Models with interactions or three-way terms use the global method.
  int <- set_interaction(base, "d1", "d2", 1)
  cs2 <- compare_scenarios(int = int, tw = set_three_way(int, "d1", "d2", "d3", 2))
  expect_equal(cs2$totals$total, c(supp_int_total, supp_int_tw_total[["2"]]), tolerance = 1e-7)
  expect_error(compare_scenarios(base, alt), "named")
  # Every scenario must be a model (e.g. an old `valuation` list is rejected).
  expect_error(compare_scenarios(base = base, valuation = list(observed = 10000)), "cm_model")
  # A scenario that cannot be run is kept with NA.
  bad <- set_association(base, "d1", "d2", 1.5, measure = "cond_prob")
  cs3 <- compare_scenarios(base = base, bad = bad)
  expect_true(is.na(cs3$totals$total[2]))
  expect_equal(unique(cs3$by_disease$scenario), "base")
})
