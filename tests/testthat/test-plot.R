# Plots and their data (R/plot.R): every plot draws, and plot_data() gives
# the data each plot draws, taken from the object without new computation.

supp_draws <- function() {
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impacts:1" = dist_normal(2.5, 0.5), "impacts:3" = dist_normal(7.5, 1)))
  deconflate(m, n_draws = 40, seed = 1)
}

test_that("every plot of a result draws and returns its data", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  res <- supp_draws()
  for (type in c("impacts", "contributions", "ranks", "rows")) {
    out <- withVisible(plot(res, type = type))
    expect_false(out$visible, info = type)
    expect_identical(out$value, plot_data(res, type), info = type)
  }
  ev <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = dist_beta(25, 75), n_draws = 20, seed = 1)
  for (type in c("impacts", "contributions", "ranks", "rows")) expect_invisible(plot(ev, type = type))
  # Interactions, several rows per term and no draws.
  w <- deconflate(read_five("impacts_welfare.csv"), n_draws = 0)
  for (type in c("impacts", "contributions", "rows")) expect_invisible(plot(w, type = type))
  # Event impacts with ratio and risk-difference rows: two panels.
  cu <- deconflate(read_five("culling.csv"), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_invisible(plot(cu, type = "rows"))
  expect_invisible(plot(res, main = "My title"))
  expect_error(plot(res, type = "bars"), class = "deconflate_invalid_input")
  expect_error(plot(w, type = "ranks"), "no draws", class = "deconflate_invalid_input")
})

test_that("plot_data() of a result: impacts, contributions, ranks and rows", {
  res <- supp_draws()
  d <- plot_data(res)
  expect_equal(names(d), c("term", "kind", "measure", "raw", "adjusted", "lower", "upper", "same_scale"))
  expect_equal(d$adjusted, res$adjusted$adjusted)
  expect_equal(d$upper, res$adjusted$upper)
  expect_true(all(d$same_scale) && all(is.na(d$measure)))
  ct <- plot_data(res, "contributions")
  expect_equal(names(ct), c("disease", "raw_contribution", "contribution", "lower", "upper", "main",
                            "interaction", "share"))
  expect_equal(ct$contribution, res$contributions$contribution)
  expect_equal(sum(ct$raw_contribution), res$totals$raw_total)
  rk <- plot_data(res, "ranks")
  expect_equal(names(rk), c("disease", "rank", "share", "mean_rank"))
  expect_equal(nrow(rk), 9L)
  expect_equal(as.vector(tapply(rk$share, rk$disease, sum)), c(1, 1, 1))
  expect_equal(rk$share[rk$disease == "d1" & rk$rank == 2], res$ranks$rank_2[res$ranks$disease == "d1"])
  rw <- plot_data(res, "rows")
  expect_equal(names(rw), c("row", "term", "among", "adjusted_for", "measure", "study", "value", "lower", "upper",
                            "fitted", "std_residual", "disagrees", "log_scale"))
  # 95% range of a row with a distribution (normal, sd 0.5); none for an exact row.
  expect_equal(c(rw$lower[1], rw$upper[1]), stats::qnorm(c(0.025, 0.975), 2.5, 0.5))
  expect_true(is.na(rw$lower[2]))
  # Event impacts: raw values on another scale are marked; ratio rows go on
  # a log axis; the range is that of each row's distribution.
  cu <- deconflate(read_five("culling.csv"), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  d <- plot_data(cu)
  expect_equal(d$same_scale, d$measure %in% c("HR", "rate_ratio"))
  rw <- plot_data(cu, "rows")
  e <- cu$estimates[cu$estimates$table == "impacts", ]
  expect_equal(c(rw$lower[1], rw$upper[1]), c(1.40, 2.16), tolerance = 0.01)  # lognormal_ci(1.74, 1.40, 2.16)
  expect_equal(rw$log_scale, e$measure != "RD")
  expect_equal(rw$upper[5], cu$model$distributions[["impacts:5"]]$q(0.975))
  expect_true(is.na(rw$upper[2]))
  expect_error(plot_data(example_supplement()), class = "deconflate_invalid_input")
})

test_that("sensitivity and threshold plots draw, and their data match the objects", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  res <- deconflate(example_supplement(), n_draws = 0)
  s <- cm_sensitivity(res)
  expect_identical(withVisible(plot(s, n = 3))$value, plot_data(s, n = 3))
  ds <- plot_data(s)
  expect_equal(names(ds), c("input", "label", "kind", "low", "high", "total_low", "total_high", "central", "swing"))
  expect_equal(ds$swing, s$inputs$swing[is.finite(s$inputs$swing)])
  expect_true(all(ds$central == res$totals$adjusted_total))
  # One search: the curves are the compared quantity itself.
  th <- cm_threshold(res, "pair:d2:d3", c(0.5, 20), conclusion = "total", target = 2.8, n_grid = 21)
  cv <- plot_data(th)
  expect_equal(names(cv), c("input", "label", "value", "item", "quantity", "level", "status"))
  expect_equal(cv$quantity - cv$level, unname(th$values[, 1]))
  expect_equal(cv$label[1], "d2-d3 OR")
  pts <- plot_data(th, "thresholds")
  expect_equal(pts$threshold, th$thresholds$threshold)
  expect_equal(pts$level, rep(2.8, nrow(pts)))
  expect_invisible(plot(th))
  rk <- cm_threshold(res, "pair:d1:d2", c(0.2, 20), conclusion = "rank")
  expect_invisible(plot(rk))
  expect_invisible(plot(rk, items = "d1 vs d2"))
  expect_error(plot(rk, items = "d9"), "Unknown items", class = "deconflate_invalid_input")
  # Sign of event impacts: hazard ratios against 1.
  ev <- deconflate(supp_hr_model(), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  sg <- cm_threshold(ev, "pair:d1:d2", c(0.5, 10), conclusion = "sign", n_grid = 11)
  cs <- plot_data(sg)
  expect_equal(cs$quantity, exp(as.vector(sg$values)))
  expect_true(all(cs$level == 1))
  expect_invisible(plot(sg))
  # A set of inputs: one panel each.
  set <- cm_threshold(res, conclusion = "rank", diseases = c("d1", "d2"), top = 2, sensitivity = s, n_grid = 11)
  cv <- plot_data(set)
  expect_equal(sort(unique(cv$input)), sort(set$inputs$input))
  expect_equal(cv$label[match(set$inputs$input, cv$input)], set$inputs$label)
  expect_equal(nrow(plot_data(set, "thresholds")), nrow(set$thresholds))
  expect_invisible(plot(set))
})

test_that("probabilistic threshold plots draw", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  m <- cm_model(example_supplement(), example_supplement()$impacts,
                distributions = list("impacts:1" = dist_normal(2.5, 0.5)))
  p <- cm_threshold_prob(m, "pair:d1:d2", c(0.5, 10), conclusion = "rank", diseases = c("d1", "d2"),
                         n_draws = 20, n_grid = 5, seed = 1)
  cv <- plot_data(p)
  expect_equal(names(cv), c("input", "label", "value", "item", "probability", "n"))
  expect_equal(cv$probability, p$curves$probability)
  expect_identical(plot_data(p, "crossings"), p$crossings)
  expect_identical(withVisible(plot(p))$value, cv)
  expect_true(is.logical(p$inputs$log_scale))
})

test_that("titles use the outcome when there is one", {
  res <- deconflate(example_supplement(), n_draws = 0)
  expect_equal(res$outcome, "yield")
  expect_equal(with_outcome("Adjusted impacts", res$outcome), "Adjusted impacts: yield")
  expect_equal(with_outcome("Adjusted impacts", NULL), "Adjusted impacts")
  expect_equal(with_units("Impact", "%"), "Impact (%)")
  expect_equal(with_units("Impact", NULL), "Impact")
  res$outcome <- NULL
  res$units <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  for (type in c("impacts", "contributions", "rows")) expect_invisible(plot(res, type = type))
})
