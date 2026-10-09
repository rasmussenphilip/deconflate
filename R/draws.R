# Uncertainty: draws of the uncertain inputs -----------------------------------
#
# deconflate() and compare_methods() propagate the uncertainty of the inputs
# that have a distribution (the model's `distributions`, and the overall risk
# of an event model) by drawing input sets, adjusting each, and summarising
# the results over the accepted draws.

# The distributions to draw from: the model's, plus the overall risk ("risk").
draw_specs <- function(model, risk_dist = NULL) {
  specs <- model$distributions %||% list()
  if (!is.null(risk_dist)) specs[["risk"]] <- risk_dist
  specs
}

# Set inputs of a model to the given values (named by key, as in
# cm_model(distributions = )). Impossible values (a probability outside
# (0, 1), a non-positive ratio, ...) stop with class deconflate_infeasible.
set_inputs <- function(model, values) {
  for (k in names(values)) {
    x <- values[[k]]
    if (!is.finite(x)) cm_abort(sprintf("Draw %s is not finite.", k), class = "deconflate_infeasible")
    parts <- strsplit(k, ":", fixed = TRUE)[[1]]
    invalid <- function(what) {
      cm_abort(sprintf("%s = %g is not a valid %s.", k, x, what), class = "deconflate_infeasible")
    }
    switch(parts[1],
      prob = {
        r <- match(parts[2], model$diseases$id)
        if (is.na(r)) cm_abort(sprintf("Unknown disease in '%s'.", k))
        type <- model$diseases$type[r]
        prob <- if (type == "incidence_rate") 1 - exp(-x) else x
        if (x < 0 || prob <= 0 || prob >= 1) invalid("disease value (it gives a probability outside (0, 1))")
        model$diseases$value[r] <- x
        model$diseases$prob[r] <- prob
      },
      assoc = {
        a <- model$associations
        r <- if (is.null(a)) NA_integer_ else match(pair_key(parts[2], parts[3]), pair_key(a$disease1, a$disease2))
        if (is.na(r)) cm_abort(sprintf("No association row for '%s'.", k))
        meas <- a$measure[r]
        ok <- switch(meas, OR = , RR = x > 0, cond_prob = x >= 0 && x <= 1,
                     RD = , phi = x >= -1 && x <= 1, TRUE)
        if (!ok) invalid(meas)
        model$associations$value[r] <- x
      },
      impact = {
        im <- model$impacts
        r <- match(parts[2], im$disease)
        if (is.na(r)) cm_abort(sprintf("No impact row for '%s'.", k))
        meas <- if (is.null(im$measure)) NA_character_ else im$measure[r]
        if (!is.na(meas)) {
          if (meas == "RD" && (x <= -1 || x >= 1)) invalid("risk difference")
          if (meas != "RD" && x <= 0) invalid(meas)
        }
        model$impacts$value[r] <- x
      },
      inter = {
        it <- model$interactions
        r <- if (is.null(it)) NA_integer_ else match(pair_key(parts[2], parts[3]), pair_key(it$disease1, it$disease2))
        if (is.na(r)) cm_abort(sprintf("No interaction row for '%s'.", k))
        model$interactions$value[r] <- x
      },
      three = {
        tw <- model$three_way
        tkey <- if (is.null(tw) || !nrow(tw)) character(0) else
          apply(tw[, c("disease1", "disease2", "disease3")], 1, function(z) paste(sort(z), collapse = "|"))
        r <- match(paste(sort(parts[2:4]), collapse = "|"), tkey)
        if (is.na(r)) cm_abort(sprintf("No three-way row for '%s'.", k))
        if (x <= 0) invalid("ratio (it must be positive)")
        model$three_way$ratio[r] <- x
      },
      cm_abort(sprintf("Unknown input key '%s'.", k))
    )
  }
  model
}

# Check the arguments that control the draws (deconflate(), compare_methods()).
check_draw_args <- function(n_draws, seed = NULL, sampling = "random", lhs_replicates = 10L) {
  whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == round(x)
  if (!whole(n_draws) || n_draws < 0) cm_abort("`n_draws` must be a single whole number, 0 or more.")
  if (n_draws == 1) cm_abort("`n_draws` must be 0 (point estimates only) or at least 2.")
  if (!is.null(seed) && (!whole(seed) || abs(seed) > .Machine$integer.max)) {
    cm_abort("`seed` must be a single whole number (as for set.seed()).")
  }
  if (!whole(lhs_replicates) || lhs_replicates < 2) cm_abort("`lhs_replicates` must be a whole number, at least 2.")
  if (identical(sampling, "lhs") && n_draws > 0 && n_draws < 4) {
    cm_abort("Latin hypercube sampling needs at least 4 draws (2 blocks of 2).")
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
  if (R < 2L) cm_abort("Latin hypercube sampling needs at least 4 draws (2 blocks of 2).")
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
          cm_abort(sprintf("risk = %g is not a valid overall risk.", risk), class = "deconflate_infeasible")
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
  list(values = values, params = params[keep, , drop = FALSE], draw = which(keep),
       block = if (is.null(uu$block)) NULL else uu$block[keep],
       n_blocks = if (is.null(uu$block)) NULL else length(unique(uu$block)),
       n_draws = n_draws, n_rejected = sum(!keep), rejections = rej, seed = seed,
       sampling = sampling, specs = specs)
}

# Summary of each quantity over the accepted draws: mean, SD, Monte Carlo
# standard error, interval and stability (see mc_stats()).
summarise_draws <- function(dr, level = 0.95) {
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  if (is.null(dr$values) || nrow(dr$values) < 2L) return(NULL)
  qn <- colnames(dr$values)
  rows <- lapply(qn, function(q) {
    g <- data.frame(x = dr$values[, q], w = 1, disease = q, method = "",
                    block = if (is.null(dr$block)) NA_integer_ else dr$block,
                    raw = NA_real_, adjusted = NA_real_, stringsAsFactors = FALSE)
    st <- mc_stats(g, probs, trim = 0.05, check_pole = FALSE, n_blocks = dr$n_blocks)$row
    data.frame(quantity = q, mean = st$mean, sd = st$sd, mcse = st$mcse,
               lower = st[[paste0("q", probs[1])]], median = st[[paste0("q", probs[2])]],
               upper = st[[paste0("q", probs[3])]], stability = st$stability,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Rows of a result's draw summary for one kind of quantity ("adjusted",
# "contribution", "total", ...), in disease order where it applies.
draw_rows <- function(x, kind, method = NULL) {
  s <- x$draws$summary
  pre <- if (is.null(method)) "" else paste0(method, "|")
  if (kind %in% c("adjusted", "contribution")) {
    ids <- x$adjusted$disease
    s[match(paste0(pre, kind, ":", ids), s$quantity), , drop = FALSE]
  } else {
    s[match(paste0(pre, kind), s$quantity), , drop = FALSE]
  }
}

# Quantities of one adjusted result, for the draws.
result_quantities <- function(r, prefix = "") {
  ids <- r$adjusted$disease
  if (inherits(r, "cm_event_result")) {
    ar <- r$attributable
    q <- c(stats::setNames(r$adjusted$adjusted, paste0(prefix, "adjusted:", ids)),
           stats::setNames(ar$summary$attributable, paste0(prefix, "total")),
           stats::setNames(ar$summary$attributable_fraction, paste0(prefix, "fraction")))
    if (!is.null(ar$by_disease)) {
      q <- c(q, stats::setNames(ar$by_disease$attributable, paste0(prefix, "contribution:", ids)))
    }
    return(q)
  }
  c(stats::setNames(r$adjusted$adjusted, paste0(prefix, "adjusted:", ids)),
    stats::setNames(r$contributions$total, paste0(prefix, "contribution:", ids)),
    stats::setNames(r$totals$adjusted_total, paste0(prefix, "total")),
    stats::setNames(r$totals$raw_sum, paste0(prefix, "raw_sum")))
}

# Draws for a deconflate() result: every draw is adjusted the same way as the
# central result, starting from its joint distribution and solution.
deconflate_draws <- function(model, specs, plan, central, risk, n_draws, seed, sampling,
                             lhs_replicates, feasibility, dots) {
  joint_start <- if (!is.null(central$joint) && identical(central$joint$backend, "exact")) central$joint
  beta_start <- if (inherits(central, "cm_event_result")) log(central$adjusted$adjusted)
  evaluate <- function(m, r) {
    res <- run_point(m, plan, risk = r %||% risk, joint = NULL, warn = FALSE,
                     feasibility = feasibility, dots = dots, joint_start = joint_start,
                     start = beta_start)
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
  if (inherits(res, "cm_event_result")) {
    if (!is.null(res$attributable)) {
      res$attributable$summary$attributable_lower <- t$lower
      res$attributable$summary$attributable_upper <- t$upper
      if (!is.null(res$attributable$by_disease)) {
        cc <- draw_rows(res, "contribution")
        res$attributable$by_disease$lower <- cc$lower
        res$attributable$by_disease$upper <- cc$upper
      }
    }
  } else {
    cc <- draw_rows(res, "contribution")
    res$contributions$lower <- cc$lower
    res$contributions$upper <- cc$upper
    res$totals$adjusted_total_lower <- t$lower
    res$totals$adjusted_total_upper <- t$upper
  }
  res
}

# Notes about the draws (rejections, stability) and the result's own notes.
draw_notes <- function(dr) {
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
  bad <- s[s$stability != "ok", , drop = FALSE]
  if (nrow(bad)) {
    out <- c(out, sprintf("Monte Carlo precision is limited for %d quantit%s (%s): increase n_draws, or use sampling = \"lhs\". See $draws$summary.",
                          nrow(bad), if (nrow(bad) == 1L) "y" else "ies",
                          paste(utils::head(sprintf("%s: %s", bad$quantity, bad$stability), 5), collapse = "; ")))
  }
  out
}

print_draws_and_notes <- function(x) {
  dr <- x$draws
  if (!is.null(dr)) {
    cat(sprintf("Uncertainty: %s intervals from %d draws (%d rejected; %s sampling; seed %d).\n",
                "95%", dr$n_draws, dr$n_rejected,
                if (identical(dr$sampling, "lhs")) "Latin hypercube" else "random", dr$seed))
  }
  notes <- x$notes
  if (length(notes)) {
    cat("\nNotes:\n")
    for (nt in notes) {
      lines <- strwrap(nt, width = max(40L, getOption("width", 80L) - 4L))
      cat("* ", paste(lines, collapse = "\n  "), "\n", sep = "")
    }
  }
  invisible(x)
}
