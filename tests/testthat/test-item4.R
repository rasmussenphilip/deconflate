# v0.4 item 4: one model, rows about combinations and subgroups, several
# rows per item. Self-consistency checks: values computed from a known
# distribution of disease combinations (and known impacts) are entered as
# rows, and the model must recover what produced them. Reference values of
# the five-disease example (inst/validation/reference_v041.py) are in
# test-import.R and test-event.R.

ids5 <- c("LAM", "MAS", "MET", "SCK", "RP")

# A known log-linear distribution of the five diseases, with pairwise terms
# and (optionally) a three-way term LAM:MAS:SCK.
known_joint <- function(three = 0.5) {
  cells <- disease_cells(5)
  colnames(cells) <- ids5
  lam <- list("LAM;MAS" = 0.6, "LAM;SCK" = 0.7, "MET;RP" = 1.8, "MET;SCK" = 0.4, "MAS;SCK" = 0.3,
              "MAS;MET" = 0.5, "LAM;MET" = 0.2, "RP;SCK" = 0.35, "MAS;RP" = 0.4, "LAM;MAS;SCK" = three)
  main <- c(LAM = -1.3, MAS = -1.2, MET = -2.3, SCK = -0.8, RP = -3.0)
  lp <- as.vector(cells %*% main[ids5])
  for (k in names(lam)) lp <- lp + lam[[k]] * all_present(cells, split_ids(k))
  p <- exp(lp - max(lp))
  ctx <- new.env()
  ctx$cells <- cells
  ctx$prob <- p / sum(p)
  ctx$masks <- list()
  ctx
}

# The value a row would have under the known distribution.
true_value <- function(ctx, disease, associated_with = NA, among = NA, measure = "OR", adjusted_for = NA) {
  a <- cm_associations(disease, associated_with, 0.5, measure = measure, among = among, adjusted_for = adjusted_for)
  row_value(association_spec(a, ids5), ctx)
}

# The nine pairs of the known distribution (RP:LAM is independent and left
# unknown) as crude odds ratios.
known_pairs <- function(ctx) {
  pr <- rbind(c("LAM", "MAS"), c("LAM", "SCK"), c("MET", "RP"), c("MET", "SCK"), c("MAS", "SCK"),
              c("MAS", "MET"), c("LAM", "MET"), c("RP", "SCK"), c("MAS", "RP"))
  data.frame(disease = pr[, 1], associated_with = pr[, 2],
             value = vapply(seq_len(nrow(pr)), function(i) true_value(ctx, pr[i, 1], pr[i, 2]), numeric(1)),
             stringsAsFactors = FALSE)
}

known_diseases <- function(ctx) cm_diseases(ids5, unname(colSums(ctx$cells * ctx$prob)))

test_that("every kind of association row recovers a known distribution", {
  skip_on_cran()
  ctx <- known_joint()
  pp <- known_pairs(ctx)
  dis <- known_diseases(ctx)
  # The pairs alone do not determine the three-way term.
  f0 <- fit_joint(cm_population(dis, cm_associations(pp$disease, pp$associated_with, pp$value)))
  expect_gt(max(abs(f0$prob - ctx$prob)), 1e-3)
  kinds <- list(
    list("LAM", "MAS;SCK", NA, "OR"),          # a combination as comparison
    list("LAM", "MAS", "SCK", "OR"),           # within a subgroup
    list("LAM", "MAS", "no SCK", "OR"),        # within its complement
    list("LAM", "MAS", "SCK", "OR_ratio"),
    list("LAM", "MAS", "no SCK", "OR_ratio"),
    list("LAM", "MAS", "SCK", "RR"),
    list("LAM", "MAS;SCK", NA, "RD"),
    list("LAM", "MAS;SCK", NA, "prevalence"),  # P(LAM | MAS and SCK)
    list("LAM;MAS;SCK", NA, NA, "prevalence"), # P(all three)
    list("MAS;SCK", "LAM", NA, "OR"))          # a combination as outcome
  for (k in kinds) {
    v <- true_value(ctx, k[[1]], k[[2]], k[[3]], k[[4]])
    a <- cm_associations(c(pp$disease, k[[1]]), c(pp$associated_with, k[[2]]), c(pp$value, v),
                         measure = c(rep("OR", nrow(pp)), k[[4]]), among = c(rep(NA, nrow(pp)), k[[3]]))
    f <- fit_joint(cm_population(dis, a))
    expect_true(f$converged, label = paste(unlist(k), collapse = " "))
    expect_lt(max(abs(f$prob - ctx$prob)), 1e-8, label = paste(unlist(k), collapse = " "))
  }
})

test_that("an association adjusted for all is a log-linear term, recovered exactly", {
  skip_on_cran()
  ctx <- known_joint(three = 0)
  # In a pairwise log-linear model, the odds ratio adjusted for every other
  # disease is exp(lambda).
  expect_equal(true_value(ctx, "LAM", "SCK", adjusted_for = "all"), exp(0.7), tolerance = 1e-10)
  pp <- known_pairs(ctx)
  adj <- pp$disease == "LAM" & pp$associated_with == "SCK"
  a <- cm_associations(pp$disease, pp$associated_with, ifelse(adj, exp(0.7), pp$value),
                       adjusted_for = ifelse(adj, "all", NA))
  f <- fit_joint(cm_population(known_diseases(ctx), a))
  expect_lt(max(abs(f$prob - ctx$prob)), 1e-8)
  # Adjusted for some diseases (LAM ~ SCK given MAS).
  v <- true_value(ctx, "LAM", "SCK", adjusted_for = "MAS")
  a2 <- cm_associations(pp$disease, pp$associated_with, ifelse(adj, v, pp$value),
                        adjusted_for = ifelse(adj, "MAS", NA))
  f2 <- fit_joint(cm_population(known_diseases(ctx), a2))
  expect_lt(max(abs(f2$prob - ctx$prob)), 1e-8)
  # Only odds ratios between single diseases can be adjusted.
  expect_error(cm_population(known_diseases(ctx),
                             cm_associations("LAM", "MAS;SCK", 2, adjusted_for = "MET")), "adjusted")
  expect_error(cm_population(known_diseases(ctx),
                             cm_associations("LAM", "SCK", 0.1, measure = "RD", adjusted_for = "MAS")), "adjusted")
})

test_that("adjustment for factors outside the model is noted and used as crude", {
  ctx <- known_joint(three = 0)
  pp <- known_pairs(ctx)
  hit <- pp$disease == "LAM" & pp$associated_with == "MAS"
  a <- cm_associations(pp$disease, pp$associated_with, pp$value,
                       adjusted_for = ifelse(hit, "parity;MAS;herd", NA))
  pop <- cm_population(known_diseases(ctx), a)
  # adjusted_for keeps the model's diseases; the row's own diseases drop out.
  notes <- adjustment_notes(association_specs(pop), "association")
  expect_match(paste(notes, collapse = " "), "parity")
  expect_match(paste(notes, collapse = " "), "crude")
  f <- fit_joint(pop)
  expect_lt(max(abs(f$prob - ctx$prob)), 1e-8)
})

# The five diseases with associations of several measures, for impact rows.
impact_population <- function() {
  dis <- cm_diseases(ids5, c(0.25, 0.30, 0.10, 0.35, 0.06),
                     measure = c("prevalence", "incidence_rate", "prevalence", "prevalence", "prevalence"))
  assoc <- cm_associations(c("LAM", "LAM", "MET", "MET", "MAS", "MAS", "LAM", "RP", "LAM", "MAS"),
                           c("MAS", "SCK", "RP", "SCK", "SCK", "MET", "MET", "SCK", "RP", "RP"),
                           c(1.8, 2.0, 6.0, 1.5, 0.35, 0.08, 0.05, 30 * 610 / (40 * 320), 1.0, 1.5),
                           measure = c("OR", "OR", "OR", "RR", "prevalence", "RD", "phi", "OR", "OR", "OR"))
  cm_population(dis, assoc)
}
b5 <- c(LAM = 4.8, MAS = 3.0, MET = 5.0, SCK = 1.5, RP = 2.5)

test_that("combination impacts (interactions) are recovered from crude and adjusted rows", {
  skip_on_cran()
  pop <- impact_population()
  d <- c("LAM;MAS" = 2.0, "LAM;MAS;SCK" = -1.0)
  # Crude rows for the diseases and LAM;MAS, a row for the triple adjusted
  # for all (its own value).
  rows <- data.frame(disease = c(ids5, "LAM;MAS", "LAM;MAS;SCK"), adjusted_for = c(rep(NA, 6), "all"))
  raw <- simulate_raw_impacts(pop, b5, d, rows = rows)
  expect_equal(raw$value[7], -1)
  r <- deconflate(cm_model(pop, raw), n_draws = 0)
  expect_equal(r$diagnostics$method, "full model")
  expect_equal(r$adjusted$disease, c(ids5, "LAM;MAS", "LAM;MAS;SCK"))
  expect_equal(r$adjusted$kind, c(rep("disease", 5), rep("interaction", 2)))
  expect_equal(r$adjusted$adjusted, unname(c(b5, d)), tolerance = 1e-8)
  # Contributions: interactions shared equally among their diseases, and
  # adding up to the aggregate.
  expect_equal(sum(r$contributions$total), r$totals$adjusted_total, tolerance = 1e-10)
  j <- fit_joint(pop)
  loss <- function(x) sum(b5[names(x)] * x) + 2 * x[["LAM"]] * x[["MAS"]] - x[["LAM"]] * x[["MAS"]] * x[["SCK"]]
  sh <- shapley_by_cell(j, loss)
  expect_equal(sh$shapley, r$contributions$total, tolerance = 1e-7)

  # A row within a subgroup tells LAM;MAS apart without a crude combination row.
  rows2 <- data.frame(disease = c(ids5, "LAM;MAS", "LAM", "LAM;MAS;SCK"),
                      adjusted_for = c(rep(NA, 5), "all", NA, NA),
                      among = c(rep(NA, 6), "MAS", NA))
  raw2 <- simulate_raw_impacts(pop, b5, d, rows = rows2)
  r2 <- deconflate(cm_model(pop, raw2), n_draws = 0)
  expect_equal(r2$adjusted$adjusted[1:7], unname(c(b5, d)), tolerance = 1e-8)
  expect_equal(r2$adjusted$rows[1], 2L)   # LAM: crude and among MAS

  # "Only vs healthy": LAM among animals with none of the other diseases.
  raw3 <- simulate_raw_impacts(pop, b5, c("LAM;MAS" = 2),
                               rows = data.frame(disease = c(ids5, "LAM", "LAM;MAS"),
                                                 among = c(rep(NA, 5), "no others", NA)))
  # With none of the others, LAM's crude difference is its own impact.
  expect_equal(raw3$value[6], 4.8, tolerance = 1e-10)
  r3 <- deconflate(cm_model(pop, raw3[-1, ]), n_draws = 0)
  expect_equal(r3$adjusted$adjusted, unname(c(b5, 2)), tolerance = 1e-8)
})

test_that("rows that do not tell terms apart are an error naming them", {
  pop <- impact_population()
  imp <- cm_impacts(c("MAS", "MET", "SCK", "RP", "LAM", "LAM;MAS"), c(1, 1, 1, 1, 1, 1),
                    adjusted_for = c(NA, NA, NA, NA, "all", "all"))
  # LAM adjusted for all and LAM;MAS adjusted for all: fine.
  expect_silent(deconflate(cm_model(pop, imp), n_draws = 0, warn = FALSE))
  # Two rows for the same comparison, both exact and different: a conflict.
  imp2 <- cm_impacts(c(ids5, "LAM;MAS", "LAM;MAS"), c(1, 1, 1, 1, 1, 1, 1.2))
  expect_error(deconflate(cm_model(pop, imp2), n_draws = 0), "LAM;MAS")
  # A combination row cannot stand in for the diseases' own rows.
  expect_error(cm_model(pop, cm_impacts(c("MAS", "MET", "SCK", "RP", "LAM;MAS"), rep(1, 5))), "LAM")
})

test_that("event combinations are recovered", {
  skip_on_cran()
  pop <- impact_population()
  j <- fit_joint(pop)
  beta <- log(c(LAM = 1.6, MAS = 1.5, MET = 1.1, SCK = 1.2, RP = 1.4))
  gam <- log(1.3)
  cells <- j$cells
  mult <- exp(as.vector(cells %*% beta[ids5]) + gam * all_present(cells, c("LAM", "MAS")))
  hr <- function(g) {
    G <- all_present(cells, g)
    N <- none_present(cells, g)
    (sum(j$prob[G] * mult[G]) / sum(j$prob[G])) / (sum(j$prob[N] * mult[N]) / sum(j$prob[N]))
  }
  vals <- c(vapply(ids5, hr, 1), hr(c("LAM", "MAS")))
  ev <- cm_model(pop, cm_impacts(c(ids5, "LAM;MAS"), vals, measure = "HR", label = "culling"))
  re <- deconflate(ev, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(re$adjusted$adjusted, unname(exp(c(beta, gam))), tolerance = 1e-7)
  expect_equal(re$adjusted$kind, c(rep("disease", 5), "interaction"))
  # HR adjusted for all is the term itself.
  ev2 <- cm_model(pop, cm_impacts(c(ids5, "LAM;MAS"), c(vals[1:5], 1.3), measure = "HR",
                                  adjusted_for = c(rep(NA, 5), "all")))
  re2 <- deconflate(ev2, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  expect_equal(re2$adjusted$adjusted, unname(exp(c(beta, gam))), tolerance = 1e-7)
})

# The supplement with two LAM-MAS-like rows for d1:d2 (each with a
# confidence interval) and Woolf's standard error for the second.
pooled_model <- function() {
  m <- example_supplement()
  or2 <- 60 * 650 / (140 * 150)
  woolf <- sqrt(1 / 60 + 1 / 140 + 1 / 150 + 1 / 650)
  a <- cm_associations(c("d1", "d1", "d1", "d2"), c("d2", "d2", "d3", "d3"), c(1.8, or2, 1, 3))
  cm_model(m$diseases, m$impacts, associations = a,
           distributions = list("associations:1" = dist_lognormal_ci(1.8, 1.3, 2.5),
                                "associations:2" = dist_lognormal_ci(or2, or2 * exp(-1.959964 * woolf),
                                                                     or2 * exp(1.959964 * woolf))))
}

test_that("several association rows are pooled, random effects by default", {
  m <- pooled_model()
  pt <- pair_tables(m)
  expect_equal(pt$measure[pt$disease1 == "d1" & pt$disease2 == "d2"], "pooled")
  r <- deconflate(m, n_draws = 0)
  h <- r$heterogeneity
  expect_equal(nrow(h), 1)
  expect_true(all(c("Q", "I2", "tau2", "pooled") %in% names(h)))
  # By hand: DerSimonian-Laird on the log odds ratios.
  y <- log(c(1.8, 60 * 650 / (140 * 150)))
  se <- c((log(2.5) - log(1.3)) / (2 * 1.959964), sqrt(1 / 60 + 1 / 140 + 1 / 150 + 1 / 650))
  w <- 1 / se^2
  fe <- sum(w * y) / sum(w)
  Q <- sum(w * (y - fe)^2)
  tau2 <- max(0, (Q - 1) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (se^2 + tau2)
  expect_equal(h$Q, Q, tolerance = 1e-4)
  expect_equal(h$tau2, tau2, tolerance = 1e-4)
  expect_equal(h$pooled, exp(sum(ws * y) / sum(ws)), tolerance = 1e-4)
  # The fitted d1:d2 table has the pooled odds ratio.
  p11 <- r$joint_pairs["d1", "d2"]
  expect_equal(unname(p11 * (1 - 0.25 + p11) / ((0.10 - p11) * (0.15 - p11))), h$pooled, tolerance = 1e-6)
  # Fixed effect.
  rf <- deconflate(m, n_draws = 0, pool = "fixed")
  expect_equal(rf$heterogeneity$pooled, exp(fe), tolerance = 1e-4)
  # Here Q < df, so tau2 = 0 and the draws need no between-study term.
  expect_equal(tau2, 0)
  rd <- deconflate(m, n_draws = 50, seed = 1)
  expect_false(any(grepl("^between:", names(rd$draws$specs))))
  # Rows that disagree: tau2 > 0, and the draws add the between-study
  # spread (one standard normal per pooled group).
  mh <- m
  mh$associations$value[2] <- 4
  mh$distributions[["associations:2"]] <- dist_lognormal_ci(4, 3.4, 4.7)
  rh <- deconflate(mh, n_draws = 50, seed = 1)
  expect_gt(rh$heterogeneity$tau2, 0)
  expect_gt(rh$heterogeneity$I2, 0.5)
  expect_equal(grep("^between:", names(rh$draws$specs), value = TRUE), "between:associations:1,2")
  # The between-study spread widens the interval of d1.
  ci <- function(x) diff(unlist(x$adjusted[1, c("lower", "upper")]))
  rh0 <- deconflate(mh, n_draws = 50, seed = 1, pool = "fixed")
  expect_gt(ci(rh), ci(rh0))
})

test_that("several impact rows are combined by weighted least squares", {
  m <- example_supplement()
  im <- cm_impacts(c("d1", "d1", "d2", "d3"), c(2.5, 1.9, 5, 7.5), adjusted_for = c(NA, "d2", NA, NA))
  m3 <- cm_model(m$diseases, im, associations = m$associations,
                 distributions = list("impacts:1" = dist_normal(2.5, 0.5), "impacts:2" = dist_normal(1.9, 0.4)))
  r <- deconflate(m3, n_draws = 0)
  expect_equal(r$adjusted$rows, c(2L, 1L, 1L))
  expect_true(is.na(r$adjusted$raw[1]))
  # Raw total with the diseases unrelated: both d1 rows are d1's own impact
  # there, fitted together (weights 1 / SD^2).
  expect_equal(r$totals$raw_total, 0.10 * (2.5 / 0.5^2 + 1.9 / 0.4^2) / (1 / 0.5^2 + 1 / 0.4^2) + 0.15 * 5 + 0.20 * 7.5,
               tolerance = 1e-10)
  e <- r$estimates[r$estimates$table == "impacts", ]
  expect_equal(nrow(e), 4)
  # The two d1 rows disagree a little: the fit lies between them, and the
  # exact rows (d2, d3) are matched.
  expect_true(all(abs(e$fitted[3:4] - e$value[3:4]) < 1e-8))
  expect_true(e$fitted[1] < 2.5 && e$fitted[2] > 1.9)
  expect_equal(e$std_residual[1:2], (e$value[1:2] - e$fitted[1:2]) / c(0.5, 0.4), tolerance = 1e-8)
  h <- r$heterogeneity[r$heterogeneity$table == "impacts", ]
  expect_equal(nrow(h), 1)
  expect_equal(h$Q, sum(e$std_residual[1:2]^2), tolerance = 1e-8)
  expect_equal(h$df, 1)
  # With the weights swapped the fit moves to the other row.
  m4 <- m3
  m4$distributions <- list("impacts:1" = dist_normal(2.5, 0.4), "impacts:2" = dist_normal(1.9, 0.5))
  r4 <- deconflate(m4, n_draws = 0)
  expect_false(isTRUE(all.equal(r4$adjusted$adjusted[1], r$adjusted$adjusted[1])))
})

test_that("functions and arguments of earlier versions say what to do instead", {
  m <- example_supplement()
  expect_error(deconflate(m, method = "global"), "method")
  expect_error(deconflate(m, feasibility = "none"), "feasibility")
  expect_error(cm_impacts("d1", 1, estimand = "snapshot_crude"), "adjusted_for")
  expect_error(cm_associations("d1", "d2", 0.3, measure = "cond_prob"), "prevalence")
  expect_error(cm_model(m$diseases, m$impacts, associations = m$associations, interactions = data.frame()),
               "combination")
  expect_error(cm_population(m$diseases, m$associations, three_way = data.frame()), "OR_ratio")
})

test_that("a hazard ratio adjusted for all within a subgroup includes the subgroup's terms", {
  dis <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
  pop <- cm_population(dis, cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(2, 1.5, 1.8)))
  # a among c, adjusted for all: b_a + delta_ac (additive), beta_a + gamma_ac (event).
  imA <- cm_impacts(c("a", "a;c", "b", "c"), c(3, 2, 1, 1), among = c("c", NA, NA, NA),
                    adjusted_for = c("all", "all", NA, NA))
  rA <- deconflate(cm_model(pop, imA), n_draws = 0, warn = FALSE)
  expect_equal(rA$adjusted$adjusted[c(1, 4)], c(1, 2), tolerance = 1e-8)
  imE <- cm_impacts(c("a", "a;c", "b", "c"), c(3, 2, 1.2, 1.2), among = c("c", NA, NA, NA),
                    adjusted_for = c("all", "all", NA, NA), measure = "HR")
  rE <- deconflate(cm_model(pop, imE), event_model = TRUE, overall_risk = 0.2, n_draws = 0, warn = FALSE)
  expect_equal(rE$adjusted$adjusted[c(1, 4)], c(1.5, 2), tolerance = 1e-8)
})

test_that("exact association rows that disagree stop every tool", {
  dis <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
  a2 <- cm_associations(c("a", "a", "a", "b"), c("b", "b", "c", "c"), c(2, 3, 1.5, 1.8))
  m2 <- cm_model(cm_population(dis, a2), cm_impacts(c("a", "b", "c"), c(2, 3, 4)))
  msg <- "rows 1 and 2 describe the same quantity"
  expect_error(deconflate(m2, n_draws = 0), msg, class = "deconflate_infeasible")
  expect_error(cm_sensitivity(m2), msg, class = "deconflate_infeasible")
  expect_error(pair_tables(m2), msg)
})

test_that("a pair row adjusted only for outside factors is the pair's crude row everywhere", {
  dis <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
  a <- cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(3, 1.5, 1.8), adjusted_for = c("parity", NA, NA))
  m <- cm_model(cm_population(dis, a), cm_impacts(c("a", "b", "c"), c(2, 3, 4)))
  expect_equal(canonical_key("assoc:a:b", m), "associations:1")
  sc <- cm_sensitivity(m, keys = "associations:1")
  expect_gt(sc$inputs$swing / abs(sc$central$total), 0.01)
  th <- cm_threshold(m, "pair:a:b", c(1.5, 6), conclusion = "total", target = 1.8, n_grid = 5)
  expect_equal(th$baseline, 3)
  expect_true(diff(range(th$scan$total)) > 0.1)
  s <- set_association(m, "b", "a", 2)
  expect_equal(nrow(s$associations), 3)
  expect_true(is.na(s$associations$adjusted_for[1]))
  # An OR_ratio row adjusted only for outside factors is used as it is.
  expect_silent(cm_population(dis, cm_associations("a", "b", 1.5, measure = "OR_ratio", among = "c",
                                                   adjusted_for = "parity")))
})

test_that("OR_ratio rows are one quantity whichever disease is the subgroup", {
  dis <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
  im <- cm_impacts(c("a", "b", "c", "a;b"), c(2, 3, 4, 1), adjusted_for = c(NA, NA, NA, "all"))
  mk <- function(d, w, among, v) {
    cm_model(cm_population(dis, cm_associations(c("a", "a", "b", d), c("b", "c", "c", w), c(2, 1.5, 1.8, v),
                                                measure = c("OR", "OR", "OR", "OR_ratio"),
                                                among = c(NA, NA, NA, among))), im)
  }
  r1 <- deconflate(mk("a", "b", "c", 1.5), n_draws = 0)
  r2 <- deconflate(mk("a", "b", "no c", 1 / 1.5), n_draws = 0)
  r3 <- deconflate(mk("c", "a", "b", 1.5), n_draws = 0)
  expect_equal(r2$joint$prob, r1$joint$prob, tolerance = 1e-9)
  expect_equal(r3$joint$prob, r1$joint$prob, tolerance = 1e-9)
  # The fitted value is reported in the row's own form.
  e2 <- r2$estimates
  expect_equal(e2$fitted[e2$table == "associations" & e2$row == 4], 1 / 1.5, tolerance = 1e-8)
  # Two uncertain rows in different forms are pooled (log scale).
  a <- cm_associations(c("a", "a", "b", "a", "a"), c("b", "c", "c", "b", "b"), c(2, 1.5, 1.8, 1.5, 0.5),
                       measure = c("OR", "OR", "OR", "OR_ratio", "OR_ratio"), among = c(NA, NA, NA, "c", "no c"))
  m <- cm_model(cm_population(dis, a), im,
                distributions = list("associations:4" = dist_lognormal_ci(1.5, 1, 2.25),
                                     "associations:5" = dist_lognormal_ci(0.5, 1 / 3, 0.75)))
  h <- deconflate(m, n_draws = 0)$heterogeneity
  expect_equal(h$rows, "4, 5")
  expect_equal(h$pooled, sqrt(1.5 * 2), tolerance = 1e-6)
  # The threshold baseline of a row among "no S" is its ratio among S.
  th <- cm_threshold(mk("a", "b", "no c", 1 / 1.5), "or_ratio:a:b|c", c(0.5, 3), conclusion = "change",
                     target = 0.5, n_grid = 3)
  expect_equal(th$baseline, 1.5)
})

test_that("an old interactions file passed as an impact table is recognised", {
  p <- cm_check_inputs(diseases = data.frame(disease = c("a", "b", "c"), value = c(0.2, 0.3, 0.25)),
                       impacts = data.frame(disease1 = "a", disease2 = "b", value = 2))
  expect_match(paste(p$problem, collapse = " "), "interactions file")
})
