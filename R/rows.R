# Rows of the input tables ------------------------------------------------------
#
# Each association row and each impact row is parsed into a "spec": the
# diseases it is about (G = `disease`, H = `associated_with`), the subgroup
# (`among`: present `apos`, absent `aneg`) and the adjustment set. See
# ?cm_associations and ?cm_impacts for what each row means.

# Members of a disease or combination ("a;b"), checked against the ids.
parse_members <- function(x, ids, where, col, required = TRUE) {
  m <- split_ids(x)
  if (!length(m)) {
    if (required) cm_abort(sprintf("%s: `%s` is missing.", where, col), class = "deconflate_invalid_input")
    return(character(0))
  }
  unk <- setdiff(m, ids)
  if (length(unk)) cm_abort(sprintf("%s: unknown disease%s in `%s`: %s.", where,
                                    if (length(unk) > 1L) "s" else "", col, paste(unk, collapse = ", ")), class = "deconflate_invalid_input")
  if (anyDuplicated(m)) cm_abort(sprintf("%s: `%s` names a disease twice.", where, col), class = "deconflate_invalid_input")
  ids[ids %in% m]   # in the order of the diseases table
}

# A subgroup: diseases present (`pos`) and absent (`neg`). "no others" makes
# every disease not named in the row (`named`) absent.
parse_among <- function(x, ids, named, where) {
  items <- split_ids(x)
  pos <- character(0)
  neg <- character(0)
  others <- FALSE
  for (it in items) {
    if (tolower(gsub("\\s+", " ", it)) == "no others") {
      others <- TRUE
    } else if (grepl("^no\\s+", it, ignore.case = TRUE)) {
      neg <- c(neg, trimws(sub("^no\\s+", "", it, ignore.case = TRUE)))
    } else {
      pos <- c(pos, it)
    }
  }
  unk <- setdiff(c(pos, neg), ids)
  if (length(unk)) cm_abort(sprintf("%s: unknown disease%s in `among`: %s.", where,
                                    if (length(unk) > 1L) "s" else "", paste(unk, collapse = ", ")), class = "deconflate_invalid_input")
  if (anyDuplicated(c(pos, neg))) cm_abort(sprintf("%s: `among` names a disease twice.", where), class = "deconflate_invalid_input")
  clash <- intersect(c(pos, neg), named)
  if (length(clash)) {
    cm_abort(sprintf("%s: `among` names %s, which the row is about.", where, paste(clash, collapse = ", ")), class = "deconflate_invalid_input")
  }
  if (others) neg <- c(neg, setdiff(ids, c(named, pos, neg)))
  list(pos = ids[ids %in% pos], neg = ids[ids %in% neg], others = others,
       single = length(items) == 1L && !others)
}

# An adjustment set: `all`, diseases of the model, and other factors (names
# that are not diseases of the model).
parse_adjusted <- function(x, ids, exclude, where) {
  s <- split_ids(x)
  if (!length(s)) return(list(all = FALSE, diseases = character(0), other = character(0)))
  if (length(s) == 1L && tolower(s) == "all") {
    return(list(all = TRUE, diseases = setdiff(ids, exclude), other = character(0)))
  }
  if (any(tolower(s) == "all")) cm_abort(sprintf("%s: `adjusted_for` is either 'all' or a list of names.", where), class = "deconflate_invalid_input")
  # The row's own diseases (and its subgroup's, which are constant within
  # it) carry no adjustment.
  d <- setdiff(intersect(s, ids), exclude)
  list(all = FALSE, diseases = ids[ids %in% d], other = setdiff(s, ids))
}

# Label of a combination in disease order ("a;b").
term_label <- function(members) paste(members, collapse = ";")

# Parse one association row (a one-row data frame) against the disease ids.
association_spec <- function(a, ids) {
  where <- sprintf("associations, row %s", a$row %||% "?")
  g <- parse_members(a$disease, ids, where, "disease")
  h <- parse_members(a$associated_with, ids, where, "associated_with", required = FALSE)
  if (length(intersect(g, h))) cm_abort(sprintf("%s: `disease` and `associated_with` share %s.", where,
                                                paste(intersect(g, h), collapse = ", ")), class = "deconflate_invalid_input")
  am <- parse_among(a$among, ids, c(g, h), where)
  adj <- parse_adjusted(a$adjusted_for, ids, c(g, h, am$pos, am$neg), where)
  m <- a$measure
  v <- a$value
  if (m != "prevalence" && !length(h)) {
    cm_abort(sprintf("%s: measure %s needs `associated_with` (the comparison).", where, m), class = "deconflate_invalid_input")
  }
  if (m == "prevalence" && !length(h) && length(g) == 1L && !length(c(am$pos, am$neg))) {
    cm_abort(sprintf("%s: the prevalence of one disease in the whole population belongs in the diseases table.", where), class = "deconflate_invalid_input")
  }
  adjusted <- length(adj$diseases) > 0L
  if (adjusted) {
    if (m != "OR" || length(g) != 1L || length(h) != 1L) {
      cm_abort(sprintf("%s: an adjusted association must be an odds ratio between two single diseases (measure OR).", where),
               class = "deconflate_unsupported", table = "associations", row = a$row)
    }
  }
  if (m == "OR_ratio") {
    if (length(g) != 1L || length(h) != 1L) {
      cm_abort(sprintf("%s: OR_ratio needs one disease in `disease` and one in `associated_with`.", where), class = "deconflate_invalid_input")
    }
    if (!am$single) {
      cm_abort(sprintf("%s: OR_ratio needs one disease in `among` (present, e.g. 'SCK', or absent, 'no SCK').", where), class = "deconflate_invalid_input")
    }
    if (length(adj$diseases) || adj$all) {
      cm_abort(sprintf("%s: OR_ratio rows cannot be adjusted for other diseases.", where), class = "deconflate_unsupported",
               table = "associations", row = a$row)
    }
  }
  # The OR_ratio is the same quantity whichever of its three diseases is the
  # subgroup, and among "no S" it is the inverse of that among "S": rows are
  # stored in one form (the three diseases in model order, the last as the
  # subgroup, present), so that rows about the same triple are pooled.
  invert <- FALSE
  if (m == "OR_ratio") {
    trip <- ids[ids %in% c(g, h, am$pos, am$neg)]
    invert <- length(am$neg) > 0L
    g <- trip[1]
    h <- trip[2]
    am$pos <- trip[3]
    am$neg <- character(0)
  }
  pos <- unique(c(g, h, am$pos))
  pos <- ids[ids %in% pos]
  kind <- if (m == "OR_ratio") "ratio" else if (m == "prevalence" && !length(h)) "prevalence" else "table"
  list(row = a$row, where = where, g = g, h = h, apos = am$pos, aneg = am$neg, adj = adj$diseases,
       adj_all = adj$all, adj_other = adj$other, adjusted = adjusted, measure = m, value = v,
       kind = kind, study = a$study %||% NA_character_, invert = invert,
       # The event whose probability the row fixes, and the identity of its
       # 2x2 table (rows of one table are pooled on the log odds ratio scale).
       event = paste(paste(pos, collapse = ","), paste(am$neg, collapse = ","), sep = "|"),
       table = if (kind == "table") {
         gh <- if (length(g) == 1L && length(h) == 1L) sort(c(g, h)) else c(term_label(g), "~", term_label(h))
         paste(paste(gh, collapse = ","), paste(am$pos, collapse = ","), paste(am$neg, collapse = ","), sep = "|")
       } else NA_character_,
       pair = kind == "table" && length(g) == 1L && length(h) == 1L && !length(c(am$pos, am$neg)) && !adjusted)
}

# Parse one impact row against the disease ids.
impact_spec <- function(r, ids) {
  where <- sprintf("impacts, row %s", r$row %||% "?")
  g <- parse_members(r$disease, ids, where, "disease")
  am <- parse_among(r$among, ids, g, where)
  adj <- parse_adjusted(r$adjusted_for, ids, c(g, am$pos, am$neg), where)
  list(row = r$row, where = where, g = g, term = term_label(g), apos = am$pos, aneg = am$neg,
       adj = adj$diseases, adj_all = adj$all, adj_other = adj$other,
       adjusted = length(adj$diseases) > 0L || adj$all,
       measure = if (is.null(r$measure)) NA_character_ else r$measure, value = r$value,
       study = r$study %||% NA_character_)
}

association_specs <- function(model) {
  a <- model$associations
  if (is.null(a) || !nrow(a)) return(list())
  ids <- model$diseases$disease
  lapply(seq_len(nrow(a)), function(r) association_spec(a[r, , drop = FALSE], ids))
}

impact_specs <- function(model) {
  im <- model$impacts
  ids <- model$diseases$disease
  lapply(seq_len(nrow(im)), function(r) impact_spec(im[r, , drop = FALSE], ids))
}

# Notes about adjustment sets that name factors outside the model.
adjustment_notes <- function(specs, table) {
  other <- Filter(function(s) length(s$adj_other) > 0L, specs)
  if (!length(other)) return(character(0))
  crude <- vapply(other, function(s) !s$adjusted, logical(1))
  out <- character(0)
  if (any(crude)) {
    out <- c(out, sprintf("%s row%s %s %s adjusted only for factors outside the model (%s); used as crude.",
                          tools_cap(table), if (sum(crude) > 1L) "s" else "",
                          paste(vapply(other[crude], `[[`, numeric(1), "row"), collapse = ", "),
                          if (sum(crude) > 1L) "are" else "is",
                          paste(unique(unlist(lapply(other[crude], `[[`, "adj_other"))), collapse = ", ")))
  }
  if (any(!crude)) {
    out <- c(out, sprintf("%s row%s %s: the adjustment for factors outside the model (%s) is ignored; the adjustment for the model's diseases is used.",
                          tools_cap(table), if (sum(!crude) > 1L) "s" else "",
                          paste(vapply(other[!crude], `[[`, numeric(1), "row"), collapse = ", "),
                          paste(unique(unlist(lapply(other[!crude], `[[`, "adj_other"))), collapse = ", ")))
  }
  out
}

tools_cap <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

# Distribution of a row, if any (key "associations:<row>" or "impacts:<row>").
row_dist <- function(model, table, row) model$distributions[[paste0(table, ":", row)]]

# Standard error of a row on a scale z(value): half the distance between the
# transformed 15.9% and 84.1% quantiles of the row's distribution (exact for a
# normal distribution on that scale); 0 for a row without a distribution
# (exact).
row_se <- function(d, zfun) {
  if (is.null(d) || identical(d$type, "fixed")) return(0)
  se <- tryCatch(abs(zfun(d$q(stats::pnorm(1))) - zfun(d$q(stats::pnorm(-1)))) / 2,
                 error = function(e) NA_real_)
  if (!is.finite(se) || se <= 0) {
    se <- tryCatch(abs(zfun(d$q(stats::pnorm(0.5))) - zfun(d$q(stats::pnorm(-0.5)))),
                   error = function(e) NA_real_)
  }
  if (!is.finite(se) || se <= 0) NA_real_ else se
}

# Pool values z (one quantity) with standard errors se (0 = exact). Exact
# rows win (the first one; the others are checked after the fit). `frozen`:
# weights and a between-study SD fixed at the central fit (for draws), with
# the draw's between-study deviate `eps`.
pool_values <- function(z, se, pool = "random", frozen = NULL, eps = 0) {
  k <- length(z)
  exact <- !is.na(se) & se == 0
  info <- list(k = k, estimate = NA_real_, Q = NA_real_, df = NA_integer_, I2 = NA_real_,
               tau2 = NA_real_, weights = rep(NA_real_, k), exact = any(exact), extra_sd = 0,
               se = NA_real_)
  if (any(exact)) {
    info$estimate <- z[which(exact)[1]]
    info$weights <- as.numeric(seq_len(k) == which(exact)[1])
    info$se <- 0
    return(info)
  }
  if (k == 1L) {
    info$estimate <- z
    info$weights <- 1
    info$se <- se
    return(info)
  }
  if (anyNA(se)) {
    # Without usable SDs every row counts equally.
    se <- rep(1, k)
  }
  if (!is.null(frozen)) {
    w <- frozen$weights
    info$estimate <- sum(w * z) / sum(w) + frozen$extra_sd * eps
    info[c("Q", "df", "I2", "tau2", "weights", "extra_sd", "se")] <-
      frozen[c("Q", "df", "I2", "tau2", "weights", "extra_sd", "se")]
    return(info)
  }
  w <- 1 / se^2
  zf <- sum(w * z) / sum(w)
  Q <- sum(w * (z - zf)^2)
  df <- k - 1L
  tau2 <- if (pool == "random") max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w))) else 0
  ws <- 1 / (se^2 + tau2)
  info$estimate <- sum(ws * z) / sum(ws)
  info$Q <- Q
  info$df <- df
  info$I2 <- if (Q > 0) max(0, (Q - df) / Q) else 0
  info$tau2 <- tau2
  info$weights <- ws
  info$se <- sqrt(1 / sum(ws))
  # Draws pool drawn values with these weights; the between-study part of
  # the variance (1 / sum(ws) minus that of the weighted mean) is added.
  info$extra_sd <- sqrt(max(0, 1 / sum(ws) - sum(ws^2 * se^2) / sum(ws)^2))
  info
}

# Probabilities a row's conversion needs, on a joint distribution (`ctx` with
# cells and their probabilities) or, for crude pair rows only, from the
# disease probabilities alone (`ctx$P`).
table_probs <- function(s, ctx) {
  if (is.null(ctx$cells)) {
    P <- ctx$P
    return(list(pR = 1, pG = P[[s$g]], pH = if (length(s$h)) P[[s$h]] else NA_real_, pA = 1,
                pHA = if (length(s$h)) P[[s$h]] else NA_real_))
  }
  m <- spec_masks(s, ctx)
  p <- ctx$prob
  pR <- sum(p[m$R])
  list(pR = pR, pG = sum(p[m$GR]) / pR, pH = if (length(s$h)) sum(p[m$HR]) / pR else NA_real_,
       pA = sum(p[m$A]), pHA = if (length(s$h)) sum(p[m$HA]) else NA_real_)
}

# Masks of a row on the cells: the subgroup A, the 2x2 table's population R
# (A, and all or none of H), G and H within R, H within A and the event E.
spec_masks <- function(s, ctx) {
  key <- paste(s$event, term_label(s$g), term_label(s$h), sep = "#")
  hit <- ctx$masks[[key]]
  if (!is.null(hit)) return(hit)
  cells <- ctx$cells
  allm <- function(v) if (length(v)) rowSums(cells[, v, drop = FALSE]) == length(v) else rep(TRUE, nrow(cells))
  nonem <- function(v) if (length(v)) rowSums(cells[, v, drop = FALSE]) == 0 else rep(TRUE, nrow(cells))
  A <- allm(s$apos) & nonem(s$aneg)
  aH <- allm(s$h)
  R <- if (length(s$h) > 1L) A & (aH | nonem(s$h)) else A
  aG <- allm(s$g)
  out <- list(A = A, R = R, GR = aG & R, HR = aH & R, HA = aH & A, E = which(aG & aH & A))
  ctx$masks[[key]] <- out
  out
}

# Log odds ratio of a 2x2 table from P(G and H), P(G), P(H).
log_or <- function(p11, pG, pH) log(p11) + log(1 - pG - pH + p11) - log(pG - p11) - log(pH - p11)

# The event probability a 2x2 row implies, given its table's probabilities.
table_joint <- function(measure, v, pG, pH) association_to_joint(measure, v, pG, pH)

# Logistic regression coefficient of H in the regression of G on H and the
# adjustment diseases within A, on the joint distribution (an adjusted odds
# ratio), and the crude log odds ratio of the same table.
adjusted_shift <- function(s, ctx) {
  m <- spec_masks(s, ctx)
  cells <- ctx$cells
  A <- m$A
  w <- ctx$prob[A]
  y <- cells[A, s$g]
  X <- cbind(1, cells[A, c(s$h, s$adj), drop = FALSE])
  # Columns constant within the subgroup carry no information.
  keep <- c(TRUE, apply(X[, -1, drop = FALSE], 2, function(x) stats::var(x) > 0))
  if (!keep[2]) {
    cm_abort(sprintf("%s: the comparison is constant within the subgroup.", s$where), class = "deconflate_singular",
             table = "associations", row = s$row)
  }
  X <- X[, keep, drop = FALSE]
  beta <- rep(0, ncol(X))
  for (it in seq_len(100L)) {
    eta <- as.vector(X %*% beta)
    mu <- stats::plogis(eta)
    g <- crossprod(X, w * (y - mu))
    H <- crossprod(X * (w * mu * (1 - mu)), X)
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step)) cm_abort(sprintf("%s: the adjusted odds ratio is not identifiable in this population.", s$where),
                                class = "deconflate_singular", table = "associations", row = s$row)
    beta <- beta + as.vector(step)
    if (max(abs(step)) < 1e-12) break
  }
  tp <- table_probs(s, ctx)
  p11 <- sum(ctx$prob[m$E]) / tp$pR
  log_or(p11, tp$pG, tp$pH) - beta[2]
}

# Solve the OR_ratio row for P(G, H, S): the odds ratio of G and H among S
# divided by that among not-S equals v (rows are stored with S present; see
# association_spec()).
ratio_target <- function(s, v, ctx) {
  S <- s$apos
  cells <- ctx$cells
  p <- ctx$prob
  allm <- function(x) rowSums(cells[, x, drop = FALSE]) == length(x)
  pS <- sum(p[allm(S)])
  pgS <- sum(p[allm(c(s$g, S))])
  phS <- sum(p[allm(c(s$h, S))])
  pgT <- sum(p[allm(s$g)])
  phT <- sum(p[allm(s$h)])
  pghT <- sum(p[allm(c(s$g, s$h))])
  a <- pgT - pgS
  b <- phT - phS
  mm <- 1 - pS
  f <- function(t) {
    log_or(t / pS, pgS / pS, phS / pS) - log_or((pghT - t) / mm, a / mm, b / mm) - log(v)
  }
  lo <- max(0, pgS + phS - pS, pghT - a, pghT - b)
  hi <- min(pgS, phS, pghT, pghT - (a + b - mm))
  eps <- 1e-12 * max(hi - lo, 1e-300)
  infeasible <- function(e = NULL) {
    cm_abort(sprintf("%s: no joint distribution has this odds ratio ratio with the other inputs.", s$where),
             class = "deconflate_infeasible", table = "associations", row = s$row)
  }
  if (!(hi - lo > 2 * eps)) infeasible()
  # (No root in the interval: the ratio is beyond what the other inputs allow.)
  tryCatch(stats::uniroot(f, c(lo + eps, hi - eps), tol = 1e-15, maxiter = 1000L)$root, error = infeasible)
}

# The value of an OR_ratio row in its stored form (among S present).
ratio_value <- function(s, v) if (isTRUE(s$invert)) 1 / v else v

# A row's event probability (target) for a value v, at the current fit.
row_target <- function(s, v, ctx, shift = NULL) {
  if (s$kind == "ratio") return(ratio_target(s, ratio_value(s, v), ctx))
  tp <- table_probs(s, ctx)
  if (s$kind == "prevalence") return(v * tp$pA)
  if (s$adjusted) {
    or <- exp(log(v) + shift)
    return(or_to_joint(or, tp$pG, tp$pH) * tp$pR)
  }
  table_joint(s$measure, v, tp$pG, tp$pH) * tp$pR
}

# A row's value on its group's pooling scale, for a value v at the current
# fit: "logOR" (the log odds ratio of the group's 2x2 table), "logit_prev"
# (the logit of a prevalence), "log_ratio" (log OR_ratio) or "logit_t" (the
# logit of the event probability, when the rows of a group differ in kind).
row_z <- function(s, v, ctx, scale, shift = NULL) {
  switch(scale,
    logOR = {
      if (s$adjusted) return(log(v) + shift)
      if (s$measure == "OR") return(log(v))
      tp <- table_probs(s, ctx)
      log_or(table_joint(s$measure, v, tp$pG, tp$pH), tp$pG, tp$pH)
    },
    logit_prev = stats::qlogis(v),
    log_ratio = log(ratio_value(s, v)),
    logit_t = stats::qlogis(row_target(s, v, ctx, shift)))
}

# From a pooled value on the group's scale back to the event probability.
z_target <- function(s, z, ctx, scale, shift = NULL) {
  switch(scale,
    logOR = {
      tp <- table_probs(s, ctx)
      or_to_joint(exp(z), tp$pG, tp$pH) * tp$pR
    },
    logit_prev = stats::plogis(z) * table_probs(s, ctx)$pA,
    log_ratio = ratio_target(s, exp(z), ctx),
    logit_t = stats::plogis(z))
}

# Group association rows by the event they fix. Each group has a pooling
# scale: the log odds ratio of a 2x2 table shared by all its rows, the logit
# of a prevalence, the log OR_ratio, or (rows of different kinds) the logit
# of the event probability.
association_groups <- function(specs) {
  if (!length(specs)) return(list())
  ev <- vapply(specs, `[[`, character(1), "event")
  keys <- unique(ev)
  lapply(keys, function(k) {
    ss <- specs[ev == k]
    kinds <- unique(vapply(ss, `[[`, character(1), "kind"))
    tabs <- unique(vapply(ss, function(s) s$table %||% NA_character_, character(1)))
    scale <- if (length(kinds) == 1L && kinds == "table" && length(tabs) == 1L) "logOR" else
      if (length(kinds) == 1L && kinds == "prevalence") "logit_prev" else
        if (length(kinds) == 1L && kinds == "ratio" && length(unique(vapply(ss, function(s) paste(c(s$g, s$h, s$apos, s$aneg), collapse = ","), character(1)))) == 1L) "log_ratio" else
          "logit_t"
    list(key = k, specs = ss, rows = vapply(ss, `[[`, numeric(1), "row"), scale = scale,
         pair = all(vapply(ss, `[[`, logical(1), "pair")), label = group_label(ss[[1]]))
  })
}

# Name of the between-study deviate of a pooled group in the draws:
# "associations:<rows>" (e.g. "associations:1,2").
between_key <- function(rows) paste0("associations:", paste(sort(rows), collapse = ","))

group_label <- function(s) {
  lab <- term_label(s$g)
  if (length(s$h)) lab <- paste(lab, "~", term_label(s$h))
  am <- c(s$apos, if (length(s$aneg)) paste("no", s$aneg))
  if (length(am)) lab <- sprintf("%s (among %s)", lab, paste(am, collapse = "; "))
  if (s$kind == "ratio") lab <- paste(lab, "OR ratio")
  lab
}

# Targets of every group at the current fit: each row's value on the group's
# scale, pooled. `values`: the rows' values (by row), `frozen`: pooling
# weights fixed at the central fit, `between`: between-study deviates of a
# draw (by group key).
group_targets <- function(groups, ctx, model, pool = "random", frozen = NULL, between = NULL) {
  out <- vector("list", length(groups))
  for (j in seq_along(groups)) {
    gr <- groups[[j]]
    z <- se <- numeric(length(gr$specs))
    for (k in seq_along(gr$specs)) {
      s <- gr$specs[[k]]
      shift <- if (s$adjusted) adjusted_shift(s, ctx)
      zf <- function(v) row_z(s, v, ctx, gr$scale, shift)
      z[k] <- zf(s$value)
      d <- row_dist(model, "associations", s$row)
      # SDs matter only for pooling several rows; one row needs only to
      # know whether it is exact.
      se[k] <- if (!is.null(frozen)) 0 else if (length(gr$specs) > 1L) row_se(d, zf) else
        if (is.null(d) || identical(d$type, "fixed")) 0 else NA_real_
      gr$specs[[k]]$shift <- shift
    }
    info <- pool_values(z, if (is.null(frozen)) se else frozen[[gr$key]]$se_rows,
                        pool = pool, frozen = if (!is.null(frozen) && !frozen[[gr$key]]$exact && length(z) > 1L) frozen[[gr$key]],
                        eps = between[[between_key(gr$rows)]] %||% 0)
    info$z <- z
    info$se_rows <- if (is.null(frozen)) se else frozen[[gr$key]]$se_rows
    s1 <- gr$specs[[if (info$exact) which(info$se_rows == 0)[1] else 1L]]
    info$target <- z_target(s1, info$estimate, ctx, gr$scale, s1$shift)
    out[[j]] <- info
  }
  out
}

# Value of a row's measure on a joint distribution (for diagnostics).
row_value <- function(s, ctx) {
  if (s$kind == "ratio") {
    S <- c(s$apos, s$aneg)
    one <- function(present) {
      s2 <- s
      s2$kind <- "table"
      s2$apos <- if (present) S else character(0)
      s2$aneg <- if (present) character(0) else S
      s2$event <- paste0(s$event, if (present) "+S" else "-S")
      tp <- table_probs(s2, ctx)
      exp(log_or(sum(ctx$prob[spec_masks(s2, ctx)$E]) / tp$pR, tp$pG, tp$pH))
    }
    return(ratio_value(s, one(TRUE) / one(FALSE)))
  }
  tp <- table_probs(s, ctx)
  m <- spec_masks(s, ctx)
  pE <- sum(ctx$prob[m$E])
  if (s$kind == "prevalence") return(pE / tp$pA)
  if (s$measure == "prevalence") return(pE / tp$pHA)
  q <- pE / tp$pR
  if (s$adjusted) return(exp(log_or(q, tp$pG, tp$pH) - adjusted_shift(s, ctx)))
  switch(s$measure,
    OR = exp(log_or(q, tp$pG, tp$pH)),
    RR = (q / tp$pH) / ((tp$pG - q) / (1 - tp$pH)),
    RD = q / tp$pH - (tp$pG - q) / (1 - tp$pH),
    phi = (q - tp$pG * tp$pH) / sqrt(tp$pG * (1 - tp$pG) * tp$pH * (1 - tp$pH)))
}

# Rake the probabilities of the cells to the event targets and the disease
# probabilities: each event's cells are scaled to its target and the other
# cells to the rest (iterative proportional fitting on events).
rake_sweep <- function(prob, events, targets, margins, P) {
  for (k in seq_along(events)) {
    m <- events[[k]]
    t <- targets[k]
    s <- sum(prob[m])
    if (s <= 0 || s >= 1) next
    f0 <- (1 - t) / (1 - s)
    prob <- prob * f0
    prob[m] <- prob[m] * (t / s) / f0
  }
  for (i in seq_along(margins)) {
    m <- margins[[i]]
    s <- sum(prob[m])
    if (s <= 0 || s >= 1) next
    f0 <- (1 - P[i]) / (1 - s)
    prob <- prob * f0
    prob[m] <- prob[m] * (P[i] / s) / f0
  }
  prob
}

rake_residual <- function(prob, events, targets, margins, P) {
  r1 <- if (length(events)) max(abs(vapply(events, function(m) sum(prob[m]), numeric(1)) - targets)) else 0
  r2 <- max(abs(vapply(margins, function(m) sum(prob[m]), numeric(1)) - P))
  max(r1, r2)
}

# Maximum-entropy joint distribution with the disease probabilities and the
# association rows: raking on events, with the targets recomputed at the
# current fit until they no longer move (a fixed point; for crude pair rows
# the targets do not depend on the fit).
fit_events <- function(model, cells, prob, groups, tol = 1e-10, max_iter = 10000L, pool = "random",
                       frozen = NULL, between = NULL) {
  ids <- model$diseases$disease
  P <- model$diseases$prob
  ctx <- new.env(parent = emptyenv())
  ctx$cells <- cells
  ctx$masks <- list()
  ctx$prob <- prob
  margins <- lapply(seq_along(ids), function(i) which(cells[, i] == 1))
  events <- lapply(groups, function(gr) spec_masks(gr$specs[[1]], ctx)$E)
  fixed <- all(vapply(groups, function(gr) gr$pair && gr$scale == "logOR", logical(1)))
  sweeps <- 0L
  info <- NULL
  targets <- numeric(0)
  converged <- FALSE
  resid <- Inf
  # Crude pair rows need only the disease probabilities.
  tctx <- if (fixed) list(P = stats::setNames(P, ids)) else ctx
  repeat {
    ctx$prob <- prob
    info <- tryCatch(group_targets(groups, tctx, model, pool, frozen, between),
                     deconflate_infeasible = function(e) e)
    if (inherits(info, "condition")) {
      # A conversion is impossible at the current fit: if it is at the
      # start, the inputs are infeasible; later, the fit has gone astray.
      stop(info)
    }
    new_t <- vapply(info, `[[`, numeric(1), "target")
    moved <- if (length(targets)) max(abs(new_t - targets), 0) else Inf
    targets <- new_t
    if (moved < tol && resid < tol) {
      converged <- TRUE
      break
    }
    if (sweeps >= max_iter) break
    inner <- if (fixed) max_iter - sweeps else min(10L, max_iter - sweeps)
    for (it in seq_len(inner)) {
      prob <- rake_sweep(prob, events, targets, margins, P)
      sweeps <- sweeps + 1L
      resid <- rake_residual(prob, events, targets, margins, P)
      if (!is.finite(resid) || resid < tol / 10) break
    }
    if (!is.finite(resid)) break
    if (fixed && resid < tol) {
      ctx$prob <- prob
      converged <- TRUE
      break
    }
  }
  ctx$prob <- prob
  if (converged) check_exact_rows(groups, info)
  list(prob = prob, converged = converged && is.finite(resid) && resid < tol, iterations = sweeps,
       residual = resid, groups = groups, pooling = info, ctx = ctx)
}

# Exact rows (no distribution) of one group must agree with each other (at
# the fit, for conversions that depend on it).
check_exact_rows <- function(groups, pooling) {
  for (j in seq_along(groups)) {
    gr <- groups[[j]]
    info <- pooling[[j]]
    exact <- which(info$se_rows == 0)
    if (length(exact) > 1L) {
      off <- abs(info$z[exact] - info$estimate) > 1e-6 * max(1, abs(info$estimate))
      if (any(off)) {
        cm_abort(sprintf("Association rows %s describe the same quantity (%s) with different values and no uncertainty, so they cannot all hold. Give them distributions (they are then pooled), or keep one.",
                         paste(gr$rows[exact], collapse = " and "), gr$label), class = "deconflate_infeasible",
                 table = "associations", rows = gr$rows[exact])
      }
    }
  }
  invisible(TRUE)
}

# Check that exact rows of one group agree with each other at the fit, and
# return the diagnostics of every association row.
association_diagnostics <- function(fit, model) {
  groups <- fit$groups
  ctx <- fit$ctx
  rows <- list()
  check_exact_rows(groups, fit$pooling)
  for (j in seq_along(groups)) {
    gr <- groups[[j]]
    info <- fit$pooling[[j]]
    for (k in seq_along(gr$specs)) {
      s <- gr$specs[[k]]
      fitted <- tryCatch(row_value(s, ctx), error = function(e) NA_real_)
      se <- info$se_rows[k]
      rows[[length(rows) + 1L]] <- data.frame(
        table = "associations", row = s$row, item = gr$label, measure = s$measure,
        study = s$study, value = s$value, fitted = fitted,
        weight = if (length(gr$specs) > 1L) info$weights[k] / sum(info$weights) else 1,
        std_residual = if (isTRUE(se > 0) && length(gr$specs) > 1L) (info$z[k] - info$estimate) / se else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out <- out[order(out$row), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Heterogeneity of groups with several rows.
association_heterogeneity <- function(groups, pooling, pool) {
  keep <- which(vapply(groups, function(g) length(g$specs) > 1L, logical(1)))
  if (!length(keep)) return(NULL)
  out <- do.call(rbind, lapply(keep, function(j) {
    gr <- groups[[j]]
    info <- pooling[[j]]
    est <- info$estimate
    data.frame(item = gr$label, rows = paste(gr$rows, collapse = ", "), k = info$k,
               scale = switch(gr$scale, logOR = "log odds ratio", logit_prev = "logit prevalence",
                              log_ratio = "log OR ratio", logit_t = "logit probability"),
               pooled = switch(gr$scale, logOR = , log_ratio = exp(est), logit_prev = , logit_t = stats::plogis(est)),
               method = if (info$exact) "exact row" else pool,
               Q = info$Q, df = info$df,
               p = if (is.finite(info$Q)) stats::pchisq(info$Q, info$df, lower.tail = FALSE) else NA_real_,
               I2 = info$I2, tau2 = info$tau2, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# Frozen pooling (for draws): the weights, SDs and between-study SD of every
# group at the central fit, by group key.
freeze_pooling <- function(groups, pooling) {
  if (!length(groups)) return(NULL)
  out <- lapply(seq_along(groups), function(j) {
    c(pooling[[j]][c("weights", "se_rows", "extra_sd", "exact", "Q", "df", "I2", "tau2", "se")],
      list(rows = groups[[j]]$rows))
  })
  names(out) <- vapply(groups, `[[`, character(1), "key")
  out
}
