# Groups of diseases (deconflate(groups = )): sums of the diseases'
# contributions, reported with the diseases staying separate in the model.

test_that("group contributions are sums of the diseases' contributions and add up to the total", {
  w <- read_five("impacts_welfare.csv")
  groups <- list(uterine = c("MET", "RP"), udder = "MAS")
  res <- deconflate(w, n_draws = 40, seed = 1, groups = groups)
  ref <- deconflate(w, n_draws = 40, seed = 1)
  # The model and the per-disease results are unchanged.
  expect_equal(res$adjusted, ref$adjusted)
  expect_equal(res$contributions, ref$contributions)
  g <- res$groups
  expect_equal(names(g), c("group", "diseases", "raw_contribution", "contribution", "lower", "upper", "main",
                           "interaction", "share"))
  # Groups first, then each disease not in a group.
  expect_equal(g$group, c("uterine", "udder", "LAM", "SCK"))
  expect_equal(g$diseases, c("MET, RP", "MAS", "LAM", "SCK"))
  ct <- res$contributions
  part <- function(col, ids) sum(ct[[col]][ct$disease %in% ids])
  for (col in c("contribution", "raw_contribution", "main", "interaction", "share")) {
    expect_equal(g[[col]][1], part(col, c("MET", "RP")), info = col)
  }
  expect_equal(sum(g$contribution), res$totals$adjusted_total)
  expect_equal(sum(g$share), 1)
  # The MET;RP interaction is wholly in the group: both halves are there.
  expect_equal(g$interaction[1], sum(ct$interaction[ct$disease %in% c("MET", "RP")]))
  # A disease alone has its own interval; a group's comes from the sum in
  # each draw (not from adding the intervals).
  expect_equal(g$lower[g$group == "LAM"], ct$lower[ct$disease == "LAM"])
  s <- res$draws$values[, "contribution:MET"] + res$draws$values[, "contribution:RP"]
  ss <- summarise_draws(list(values = cbind(uterine = s), block = res$draws$block, n_blocks = res$draws$n_blocks))
  expect_equal(c(g$lower[1], g$upper[1]), c(ss$lower, ss$upper))
  # Ranks of the groups over the draws.
  gr <- res$group_ranks
  expect_equal(names(gr)[1:2], c("group", "mean_rank"))
  expect_equal(sort(gr$group), sort(g$group))
  expect_equal(unname(rowSums(as.matrix(gr[, grep("^rank_", names(gr))]))), rep(1, 4))
  # contribution_table() names each disease's group.
  expect_equal(contribution_table(res)$group, c(NA, "udder", "uterine", NA, "uterine"))
  expect_true(all(is.na(contribution_table(ref)$group)))
  expect_null(ref$groups)
})

test_that("groups of event impacts add up to the attributable risk", {
  cu <- read_five("culling.csv")
  res <- deconflate(cu, event_model = TRUE, overall_risk = 0.25, n_draws = 0, groups = list(uterine = c("MET", "RP")))
  g <- res$groups
  ct <- res$contributions
  expect_equal(g$contribution[1], sum(ct$contribution[ct$disease %in% c("MET", "RP")]))
  expect_equal(sum(g$contribution), sum(ct$contribution))
  expect_true(all(is.na(g$lower)))
  expect_null(res$group_ranks)
})

test_that("groups are checked", {
  w <- read_five("impacts_welfare.csv")
  err <- function(groups, pattern) {
    expect_error(deconflate(w, n_draws = 0, groups = groups), pattern, class = "deconflate_invalid_input")
  }
  err(list(a = "XX"), "not in the model: XX")
  err(list(LAM = "MAS"), "cannot have the name of a disease")
  err(list(a = "MET", b = c("RP", "MET")), "one group only")
  err(list("MET"), "needs a name")
  err("MET", "named list")
  err(list(a = c("MET", "MET;RP")), "not in the model")
  err(list(a = c("MET", "RP"), a = "LAM"), "unique")
  # A named character vector is not taken as groups (c(a = c("x", "y"))
  # would silently become groups a1 and a2).
  err(c(uterine = "MET"), "named list")
})

test_that("groups are printed, summarised and plotted", {
  w <- read_five("impacts_welfare.csv")
  res <- deconflate(w, n_draws = 30, seed = 2, groups = list(uterine = c("MET", "RP")))
  out <- capture.output(print(res))
  expect_true(any(grepl("Groups (sums of the diseases' contributions)", out, fixed = TRUE)))
  expect_true(any(grepl("uterine", out)))
  expect_true(all(nchar(out) <= 80))
  so <- capture.output(print(summary(res)))
  expect_true(any(grepl("Contributions by group", so, fixed = TRUE)))
  expect_true(any(grepl("Rank of the groups over the draws", so, fixed = TRUE)))
  expect_true(all(nchar(so) <= 80))
  expect_equal(summary(res)$groups, res$groups)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(withVisible(plot(res, type = "groups"))$value, res$groups)
  gr <- plot_data(res, "group_ranks")
  expect_equal(names(gr), c("group", "rank", "share", "mean_rank"))
  expect_invisible(plot(res, type = "group_ranks"))
  ref <- deconflate(w, n_draws = 0)
  expect_error(plot(ref, type = "groups"), "no groups", class = "deconflate_invalid_input")
  expect_error(plot_data(ref, "group_ranks"), "no groups", class = "deconflate_invalid_input")
})
