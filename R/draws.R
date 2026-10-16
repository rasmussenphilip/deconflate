# Uncertainty: draws of the uncertain inputs -----------------------------------
#
# deconflate() propagates the uncertainty of the inputs
# that have a distribution (the model's `distributions`, and the overall or
# baseline risk of an event model) by drawing input sets, adjusting each, and summarising
# the results over the accepted draws.

# The distributions to draw from: the model's, the overall or baseline risk ("risk"),
# and a standard normal between-study deviate ("between:associations:<rows>",
# e.g. "between:associations:1,2") for each pooled group of association rows
# with a random-effects term.
draw_specs <- function(model, risk_dist = NULL, frozen = NULL) {
  specs <- model$distributions %||% list()
  if (!is.null(risk_dist)) specs[["risk"]] <- risk_dist
  for (g in names(frozen)) {
    if (isTRUE(frozen[[g]]$extra_sd > 0)) specs[[paste0("between:", between_key(frozen[[g]]$rows))]] <- dist_normal(0, 1)
  }
  specs
}

# Set inputs of a model to the given values (named by canonical key, see
# canonical_key()). Impossible values (a probability outside (0, 1), a
# non-positive ratio, ...) stop with class deconflate_infeasible.
set_inputs <- function(model, values) {
  for (k in names(values)) {
    x <- values[[k]]
    if (!is.finite(x)) cm_abort(sprintf("Draw %s is not finite.", k), class = "deconflate_infeasible")
    if (startsWith(k, "between:")) {
      model$between[[substring(k, 9L)]] <- x
      next
    }
    if (!grepl("^(prob|associations|impacts):", k)) k <- canonical_key(k, model)
    parts <- strsplit(k, ":", fixed = TRUE)[[1]]
    invalid <- function(what) {
      cm_abort(sprintf("%s = %g is not a valid %s.", k, x, what), class = "deconflate_infeasible")
    }
    switch(parts[1],
      prob = {
        r <- match(parts[2], model$diseases$disease)
        if (is.na(r)) cm_abort(sprintf("Unknown disease in '%s'.", k), class = "deconflate_invalid_input")
        prob <- if (model$diseases$measure[r] == "incidence_rate") 1 - exp(-x) else x
        if (x < 0 || prob <= 0 || prob >= 1) invalid("disease value (it gives a probability outside (0, 1))")
        model$diseases$value[r] <- x
        model$diseases$prob[r] <- prob
      },
      associations = {
        a <- model$associations
        r <- if (is.null(a)) NA_integer_ else match(as.integer(parts[2]), a$row)
        if (is.na(r)) cm_abort(sprintf("No association row for '%s'.", k), class = "deconflate_invalid_input")
        meas <- a$measure[r]
        ok <- switch(meas, OR = , RR = , OR_ratio = x > 0, prevalence = x > 0 && x < 1,
                     RD = , phi = x >= -1 && x <= 1, TRUE)
        if (!ok) invalid(meas)
        model$associations$value[r] <- x
      },
      impacts = {
        im <- model$impacts
        r <- if (is.null(im)) NA_integer_ else match(as.integer(parts[2]), im$row)
        if (is.na(r)) cm_abort(sprintf("No impact row for '%s'.", k), class = "deconflate_invalid_input")
        meas <- im$measure[r]
        if (!is.na(meas)) {
          if (meas == "RD" && (x <= -1 || x >= 1)) invalid("risk difference")
          if (meas != "RD" && x <= 0) invalid(meas)
        }
        model$impacts$value[r] <- x
      },
      cm_abort(sprintf("Unknown input key '%s'.", k), class = "deconflate_invalid_input")
    )
  }
  model
}

# Check the arguments that control the draws.
check_draw_args <- function(n_draws, seed = NULL, sampling = "random", lhs_replicates = 10L) {
  whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  if (!whole(n_draws) || n_draws < 0) cm_abort("`n_draws` must be a single whole number, 0 or more.", class = "deconflate_invalid_input")
  if (n_draws == 1) cm_abort("`n_draws` must be 0 (point estimates only) or at least 2.", class = "deconflate_invalid_input")
  if (!is.null(seed) && (!whole(seed) || abs(seed) > .Machine$integer.max)) {
    cm_abort("`seed` must be a single whole number (as for set.seed()).", class = "deconflate_invalid_input")
  }
  if (!whole(lhs_replicates) || lhs_replicates < 2) cm_abort("`lhs_replicates` must be a whole number, at least 2.", class = "deconflate_invalid_input")
  if (identical(sampling, "lhs") && n_draws > 0 && n_draws < 4) {
    cm_abort("Latin hypercube sampling needs at least 4 draws (2 blocks of 2).", class = "deconflate_invalid_input")
  }
  invisible(TRUE)
}

# Uniforms for the draws: simple random, or Latin hypercube in replicate
# blocks (the Monte Carlo error is then estimated from the block means).
draw_uniforms <- function(n_draws, keys, sampling, lhs_replicates) {
  k <- length(keys)
  if (sampling == "random") {
    U <- matrix(stats::runif(n_draws * k), n_draws, k, dimnames = list(NULL, keys))
    return(list(U = U, block = NULL))
  }
  R <- min(as.integer(lhs_replicates), n_draws %/% 2L)
  if (R < 2L) cm_abort("Latin hypercube sampling needs at least 4 draws (2 blocks of 2).", class = "deconflate_invalid_input")
  block <- sort(rep(seq_len(R), length.out = n_draws))
  U <- matrix(NA_real_, n_draws, k, dimnames = list(NULL, keys))
  for (b in seq_len(R)) {
    idx <- which(block == b)
    nb <- length(idx)
    for (kk in keys) U[idx, kk] <- (sample.int(nb) - stats::runif(nb)) / nb
  }
  list(U = U, block = block)
}

# Run the draws. `evaluate(model, risk)` adjusts one drawn model and returns a
# named numeric vector of quantities (the same names for every draw). Draws
# that fail are rejected and counted by type.
run_draws <- function(model, specs, n_draws, seed, sampling, lhs_replicates, evaluate) {
  keys <- names(specs)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  set.seed(seed)
  uu <- draw_uniforms(n_draws, keys, sampling, lhs_replicates)
  U <- uu$U
  vals <- vector("list", n_draws)
  params <- matrix(NA_real_, n_draws, length(keys), dimnames = list(NULL, keys))
  rej_draw <- integer(0)
  rej_type <- character(0)
  rej_reason <- character(0)
  for (d in seq_len(n_draws)) {
    x <- vapply(keys, function(k) specs[[k]]$q(U[d, k]), numeric(1))
    params[d, ] <- x
    res <- tryCatch(
      withCallingHandlers({
        m <- set_inputs(model, x[keys != "risk"])
        risk <- if ("risk" %in% keys) x[["risk"]] else NULL
        if (!is.null(risk) && !(risk > 0 && risk < 1)) {
          cm_abort(sprintf("risk = %g is not a valid risk (between 0 and 1).", risk), class = "deconflate_infeasible")
        }
        q <- evaluate(m, risk)
        if (!all(is.finite(q))) cm_abort("Non-finite results.", class = "deconflate_nonfinite")
        q
      }, deconflate_nonconvergence = function(w) {
        if (inherits(w, "warning")) invokeRestart("muffleWarning")
      }),
      deconflate_error = function(e) e)
    if (inherits(res, "condition")) {
      rej_draw <- c(rej_draw, d)
      rej_type <- c(rej_type, condition_type(res))
      rej_reason <- c(rej_reason, conditionMessage(res))
    } else {
      vals[[d]] <- res
    }
  }
  keep <- !vapply(vals, is.null, logical(1))
  values <- if (any(keep)) do.call(rbind, vals[keep]) else NULL
  rej <- data.frame(draw = rej_draw, type = rej_type, reason = rej_reason, stringsAsFactors = FALSE)
  list(values = values, params = params[keep, , drop = FALSE], uniforms = U[keep, , drop = FALSE],
       draw = which(keep),
       block = if (is.null(uu$block)) NULL else uu$block[keep],
       n_blocks = if (is.null(uu$block)) NULL else length(unique(uu$block)),
       n_draws = n_draws, n_rejected = sum(!keep), rejections = rej, seed = seed,
       sampling = sampling, specs = specs)
}

# Summary of each quantity over the accepted draws: mean, SD, Monte Carlo
# standard error, interval and stability (see draw_stats()).
summarise_draws <- function(dr, level = 0.95) {
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  if (is.null(dr$values) || nrow(dr$values) < 2L) return(NULL)
  qn <- colnames(dr$values)
  rows <- lapply(qn, function(q) {
    g <- data.frame(x = dr$values[, q], block = if (is.null(dr$block)) NA_integer_ else dr$block)
    # The raw total can be missing in some draws (a failed fit with the
    # diseases unrelated); it is summarised over the others.
    g <- g[is.finite(g$x), , drop = FALSE]
    if (nrow(g) < 2L) {
      return(data.frame(quantity = q, mean = NA_real_, sd = NA_real_, mcse = NA_real_, lower = NA_real_,
                        median = NA_real_, upper = NA_real_, stability = NA_character_,
                        stringsAsFactors = FALSE))
    }
    st <- draw_stats(g$x, g$block, probs, n_blocks = dr$n_blocks)
    data.frame(quantity = q, mean = st$mean, sd = st$sd, mcse = st$mcse,
               lower = st$q[1], median = st$q[2], upper = st$q[3], stability = st$stability,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Summary and stability of one quantity over the accepted draws (equal
# weights; `block`: the Latin hypercube replicate of each draw, or NA).
draw_stats <- function(x, block, probs, n_blocks = NULL) {
  n <- length(x)
  w <- rep(1 / n, n)
  mu <- sum(w * x)
  s <- sqrt(sum(w * (x - mu)^2))
  q <- weighted_quantile(x, w, probs)
  if (all(!is.na(block))) {
    # The pooled mean is a ratio estimator, mu = sum_b S_b / sum_b W_b, over
    # the R replicate blocks (S_b: weighted sum, W_b: weight of block b).
    # Its standard error is sqrt(sum_b (S_b - mu W_b)^2 / (R (R - 1))) / mean(W_b).
    # Blocks whose draws were all rejected are replicates with S_b = W_b = 0:
    # they count in R. With equal block weights this is
    # sd(block means) / sqrt(R). Fewer than two blocks with positive weight
    # give no precision information.
    ix_b <- split(seq_along(x), block)
    Wb <- vapply(ix_b, function(ix) sum(w[ix]), numeric(1))
    Sb <- vapply(ix_b, function(ix) sum(w[ix] * x[ix]), numeric(1))
    R <- max(n_blocks %||% 0L, length(ix_b))
    mcse <- if (sum(Wb > 0) >= 2L && R >= 2L) {
      sqrt(sum((Sb - mu * Wb)^2) / (R * (R - 1))) / (sum(Wb) / R)
    } else NA_real_
  } else {
    mcse <- sqrt(sum(w^2 * (x - mu)^2))
  }
  rel <- if (isTRUE(abs(mu) > 1e-12)) mcse / abs(mu) else NA_real_
  cdev <- w * (x - mu)^2
  k <- max(1L, ceiling(0.01 * n))
  tail_share <- if (isTRUE(sum(cdev) > 0)) sum(sort(cdev, decreasing = TRUE)[seq_len(k)]) / sum(cdev) else 0
  stability <- if (n >= 50 && tail_share > 0.6) {
    "heavy_tail"
  } else if (is.na(mcse) || (is.na(rel) && mcse > 0)) {
    "insufficient_info"
  } else if (!is.na(rel) && rel > 0.05) {
    "imprecise"
  } else {
    "ok"
  }
  list(mean = mu, sd = s, mcse = mcse, q = q, stability = stability)
}

# Rows of a result's draw summary for one kind of quantity ("adjusted",
# "contribution", "total", ...), in the order of the result's tables.
draw_rows <- function(x, kind) {
  s <- x$draws$summary
  if (kind == "adjusted") {
    s[match(paste0(kind, ":", x$adjusted$disease), s$quantity), , drop = FALSE]
  } else if (kind == "contribution") {
    s[match(paste0(kind, ":", x$contributions$disease), s$quantity), , drop = FALSE]
  } else {
    s[match(kind, s$quantity), , drop = FALSE]
  }
}

# Quantities of one adjusted result, for the draws.
result_quantities <- function(r) {
  terms <- r$adjusted$disease
  if (inherits(r, "cm_event_result")) {
    ar <- r$attributable
    q <- c(stats::setNames(r$adjusted$adjusted, paste0("adjusted:", terms)),
           stats::setNames(ar$summary$attributable, "total"),
           stats::setNames(ar$summary$attributable_fraction, "fraction"))
    ct <- r$contributions
    if (!is.null(ct)) q <- c(q, stats::setNames(ct$contribution, paste0("contribution:", ct$disease)))
    # Always present when computed (NA in a draw where it fails), so that
    # every draw has the same quantities.
    if (!is.null(ar$summary$raw_attributable)) {
      q <- c(q, stats::setNames(ar$summary$raw_attributable, "raw_total"))
    }
    return(q)
  }
  q <- c(stats::setNames(r$adjusted$adjusted, paste0("adjusted:", terms)),
         stats::setNames(r$contributions$contribution, paste0("contribution:", r$contributions$disease)),
         stats::setNames(r$totals$adjusted_total, "total"))
  if (!is.null(r$totals$raw_total)) q <- c(q, stats::setNames(r$totals$raw_total, "raw_total"))
  q
}

# Draws for a deconflate() result: every draw is adjusted the same way as the
# central result, starting from its joint distribution and solution, with
# the pooling weights and impact-row SDs of the central estimate.
deconflate_draws <- function(model, specs, plan, central, risk, n_draws, seed, sampling,
                             lhs_replicates, dots) {
  joint_start <- if (!is.null(central$joint) && identical(central$joint$backend, "exact")) central$joint
  beta_start <- if (inherits(central, "cm_event_result")) log(central$adjusted$adjusted)
  frozen <- list(associations = central$pooling$frozen, impacts = central$row_se)
  evaluate <- function(m, r) {
    res <- run_point(m, plan, risk = r %||% risk, joint = NULL, warn = FALSE, dots = dots,
                     joint_start = joint_start, start = beta_start, frozen = frozen)
    result_quantities(res)
  }
  dr <- run_draws(model, specs, n_draws, seed, sampling, lhs_replicates, evaluate)
  dr$summary <- summarise_draws(dr)
  dr
}

# Add the interval columns of the draws to a result's tables.
attach_intervals <- function(res) {
  if (is.null(res$draws$summary)) return(res)
  a <- draw_rows(res, "adjusted")
  res$adjusted$lower <- a$lower
  res$adjusted$upper <- a$upper
  t <- draw_rows(res, "total")
  if (!is.null(res$contributions)) {
    cc <- draw_rows(res, "contribution")
    res$contributions$lower <- cc$lower
    res$contributions$upper <- cc$upper
  }
  if (inherits(res, "cm_event_result")) {
    if (!is.null(res$attributable)) {
      res$attributable$summary$attributable_lower <- t$lower
      res$attributable$summary$attributable_upper <- t$upper
    }
  } else {
    res$totals$adjusted_total_lower <- t$lower
    res$totals$adjusted_total_upper <- t$upper
  }
  res
}

# Notes about the draws (rejections, stability) and the result's own notes.
draw_notes <- function(dr, event = FALSE) {
  if (is.null(dr)) return(character(0))
  out <- character(0)
  if (dr$n_rejected == dr$n_draws) {
    return(sprintf("All %d draws were rejected (see $draws$rejections), so there are no intervals.", dr$n_draws))
  }
  if (dr$n_draws - dr$n_rejected < 2L) {
    return(sprintf("Only 1 of %d draws was accepted (see $draws$rejections), too few for intervals.", dr$n_draws))
  }
  share <- dr$n_rejected / dr$n_draws
  if (share > 0.1) {
    tab <- table(dr$rejections$type)
    out <- c(out, sprintf("%.0f%% of the draws were rejected (%s): the input distributions often combine values that cannot hold together. The intervals describe the accepted draws only; see $draws$rejections.",
                          100 * share, paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = ", ")))
  }
  s <- dr$summary
  bad <- s$quantity[!is.na(s$stability) & s$stability != "ok"]
  if (length(bad)) {
    out <- c(out, sprintf("Too few draws for precise intervals of %s: increase n_draws, or use sampling = \"lhs\" (see $draws$summary).",
                          imprecise_text(bad, event)))
  }
  out
}

# The quantities of the draw summary ("adjusted:LAM", "contribution:MET",
# "total", ...) in plain words: "the adjusted impacts of LAM and LAM;MAS
# and the contributions of MET and SCK".
imprecise_text <- function(q, event = FALSE) {
  and_list <- function(x) if (length(x) < 2L) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
  kind <- sub(":.*$", "", q)
  what <- sub("^[^:]*:", "", q)
  parts <- character(0)
  adj <- what[kind == "adjusted"]
  if (length(adj)) {
    parts <- c(parts, sprintf("the adjusted %s of %s",
                              if (event) if (length(adj) == 1L) "hazard ratio" else "hazard ratios"
                              else if (length(adj) == 1L) "impact" else "impacts", and_list(adj)))
  }
  con <- what[kind == "contribution"]
  if (length(con)) {
    parts <- c(parts, sprintf("the contribution%s of %s", if (length(con) == 1L) "" else "s", and_list(con)))
  }
  single <- c(total = if (event) "the attributable risk" else "the adjusted total",
              raw_total = if (event) "the raw attributable risk" else "the raw total",
              fraction = "the attributable fraction")
  parts <- c(parts, unname(single[intersect(names(single), kind)]))
  other <- setdiff(kind, c("adjusted", "contribution", names(single)))
  if (length(other)) parts <- c(parts, q[kind %in% other])
  and_list(parts)
}

# The draws line of print() and summary().
uncertainty_line <- function(dr) {
  if (is.null(dr)) return(character(0))
  sprintf("95%% intervals from %d draws (%d rejected; %s sampling; seed %d).",
          dr$n_draws, dr$n_rejected, if (identical(dr$sampling, "lhs")) "Latin hypercube" else "random", dr$seed)
}

# Notes as a bulleted list, wrapped to the console width.
print_notes <- function(notes, title = "Notes:") {
  if (!length(notes)) return(invisible(NULL))
  cat("\n", title, "\n", sep = "")
  for (nt in notes) {
    lines <- strwrap(nt, width = max(40L, getOption("width", 80L) - 4L))
    cat("* ", paste(lines, collapse = "\n  "), "\n", sep = "")
  }
  invisible(NULL)
}

# How often each disease ranks 1st, 2nd, ... by its contribution (for event
# impacts, its share of the attributable risk) over the accepted draws, and
# its mean rank. NULL without draws.
draw_ranks <- function(res) {
  v <- res$draws$values
  if (is.null(v) || nrow(v) < 2L) return(NULL)
  ids <- res$contributions$disease
  cols <- paste0("contribution:", ids)
  if (is.null(ids) || !all(cols %in% colnames(v))) return(NULL)
  C <- v[, cols, drop = FALSE]
  C <- C[apply(C, 1, function(x) all(is.finite(x))), , drop = FALSE]
  if (nrow(C) < 2L) return(NULL)
  R <- t(apply(C, 1, function(x) rank(-x, ties.method = "first")))
  if (length(ids) == 1L) R <- matrix(R, ncol = 1L)
  k <- length(ids)
  P <- vapply(seq_len(k), function(r) colMeans(R == r), numeric(k))
  P <- matrix(P, nrow = k, dimnames = list(NULL, paste0("rank_", seq_len(k))))
  out <- data.frame(disease = ids, mean_rank = unname(colMeans(R)), P, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
