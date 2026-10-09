# Monte Carlo for the reproduction of Rasmussen et al. (2024) -------------------
#
# Internal. reproduce_rasmussen_2024() draws the 2024 inputs with a batch
# sampler (shared disease and association draws for every analysis; see
# sampler.R) and adjusts each analysis with the published approximation.
# New analyses use deconflate(n_draws = ), whose draws are in draws.R; the
# summary statistics (mc_stats()) are shared.

# Run a batch sampler: `n_draws` draws of the shared population inputs, each
# used by every analysis, adjusted with `method`.
mc_batch <- function(bs, n_draws, method, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  nms <- names(bs$samplers)
  accs <- stats::setNames(lapply(nms, function(nm) new_accumulator(n_draws)), nms)
  for (d in seq_len(n_draws)) {
    pop_vals <- vapply(bs$population_keys, function(k) bs$population_specs[[k]]$r(1), numeric(1))
    for (nm in nms) {
      s <- bs$samplers[[nm]]
      res <- run_draw(function() s(d, values = pop_vals), method, ...)
      accs[[nm]] <- record_draw(accs[[nm]], d, res, method)
    }
  }
  out <- lapply(nms, function(nm) finish_mc(accs[[nm]], n_draws, method))
  names(out) <- nms
  structure(list(analyses = out, n_draws = n_draws, population_keys = bs$population_keys),
            class = "cm_mc_batch")
}

# Run one draw: sample a model and adjust it with every method. Returns the
# results or a classed failure.
run_draw <- function(make_model, method, ...) {
  tryCatch(
    withCallingHandlers({
      model <- make_model()
      rs <- lapply(method, function(m) adjust_impacts(model, method = m, warn = FALSE, ...))
      bad <- !vapply(rs, result_is_finite, logical(1))
      if (any(bad)) {
        cm_abort(sprintf("Non-finite adjusted impacts (method %s).", paste(method[bad], collapse = ", ")),
                 class = "deconflate_nonfinite")
      }
      list(model = model, results = rs)
    }, deconflate_nonconvergence = function(w) {
      if (inherits(w, "warning")) invokeRestart("muffleWarning")
    }),
    deconflate_error = function(e) e
  )
}

new_accumulator <- function(n) {
  list(draws = vector("list", n), totals = vector("list", n),
       rej_draw = integer(0), rej_type = character(0), rej_reason = character(0),
       label = NULL, units = NULL)
}

record_draw <- function(acc, d, res, method) {
  if (inherits(res, "condition")) {
    acc$rej_draw <- c(acc$rej_draw, d)
    acc$rej_type <- c(acc$rej_type, condition_type(res))
    acc$rej_reason <- c(acc$rej_reason, conditionMessage(res))
    return(acc)
  }
  acc$draws[[d]] <- do.call(rbind, lapply(seq_along(method), function(j) {
    r <- res$results[[j]]
    data.frame(draw = d, method = method[j], disease = r$adjusted$disease,
               raw = r$adjusted$raw, adjusted = r$adjusted$adjusted,
               contribution = r$contributions$total, stringsAsFactors = FALSE)
  }))
  acc$totals[[d]] <- do.call(rbind, lapply(seq_along(method), function(j) {
    r <- res$results[[j]]
    data.frame(draw = d, method = method[j], raw_sum = r$totals$raw_sum,
               adjusted_total = r$totals$adjusted_total, stringsAsFactors = FALSE)
  }))
  acc$label <- acc$label %||% res$results[[1]]$label
  acc$units <- acc$units %||% res$results[[1]]$units
  acc
}

finish_mc <- function(acc, n_draws, method) {
  keep <- !vapply(acc$draws, is.null, logical(1))
  structure(list(
    draws = do.call(rbind, acc$draws[keep]),
    totals = do.call(rbind, acc$totals[keep]),
    n_draws = n_draws, n_rejected = length(acc$rej_draw),
    rejections = data.frame(draw = acc$rej_draw, type = acc$rej_type, reason = acc$rej_reason,
                            stringsAsFactors = FALSE),
    method = method, label = acc$label, units = acc$units
  ), class = "cm_mc")
}

# Summary of adjusted impacts by disease and method over the accepted draws
# of one analysis (equal weights).
mc_summary <- function(object, probs = c(0.025, 0.5, 0.975), trim = 0.05) {
  if (is.null(object$draws)) cm_abort("All draws were rejected.")
  d <- object$draws
  d$x <- d$adjusted
  d$w <- 1
  d$block <- NA_integer_
  grp <- list(factor(d$disease, levels = unique(d$disease)),
              factor(d$method, levels = unique(d$method)))
  groups <- split(d, grp, drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(g) mc_stats(g, probs, trim)$row))
  rownames(out) <- NULL
  out
}

# Summary and stability statistics for one group of draws (columns x, w,
# block; and disease, method, raw, adjusted for the labels and the pole check
# of the published approximation).
mc_stats <- function(g, probs, trim, check_pole = TRUE, n_blocks = NULL) {
  w <- g$w / sum(g$w)
  x <- g$x
  mu <- sum(w * x)
  s <- sqrt(sum(w * (x - mu)^2))
  q <- weighted_quantile(x, w, probs)
  lim <- weighted_quantile(x, w, c(trim, 1 - trim))
  inside <- x >= lim[1] & x <= lim[2]
  tmean <- if (any(inside)) sum(w[inside] * x[inside]) / sum(w[inside]) else NA_real_
  if (all(!is.na(g$block))) {
    # The pooled mean is a ratio estimator, mu = sum_b S_b / sum_b W_b, over
    # the R replicate blocks (S_b: weighted sum, W_b: weight of block b).
    # Its standard error is sqrt(sum_b (S_b - mu W_b)^2 / (R (R - 1))) / mean(W_b).
    # Blocks whose draws were all rejected are replicates with S_b = W_b = 0:
    # they count in R. With equal block weights this is
    # sd(block means) / sqrt(R). Fewer than two blocks with positive weight
    # give no precision information.
    ix_b <- split(seq_along(x), g$block)
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
  k <- max(1L, ceiling(0.01 * length(x)))
  tail_share <- if (isTRUE(sum(cdev) > 0)) sum(sort(cdev, decreasing = TRUE)[seq_len(k)]) / sum(cdev) else 0
  # The published approximation is m^2 / (m + c). With c = 0 (no associated
  # impacts) it reduces to m and the singularity is removable. A pole lies
  # inside the sampled inputs when m + c (recovered as m^2 / adjusted) takes
  # both signs and differs in sign from m in some draws.
  pole <- NA_real_
  if (check_pole && identical(g$method[1], "published") && all(!is.na(g$raw))) {
    m <- g$raw
    a <- g$adjusted
    ok <- is.finite(m) & is.finite(a) & abs(m) > 1e-12 & abs(a) > 1e-12
    if (sum(ok) > 1) {
      den <- m[ok]^2 / a[ok]
      both <- any(den > 0) && any(den < 0)
      flip <- mean(sign(den) != sign(m[ok]))
      pole <- if (both && flip > 0) flip else 0
    }
  }
  stability <- if (!is.na(pole) && pole > 0) {
    "possible_pole"
  } else if (length(x) >= 50 && tail_share > 0.6) {
    "heavy_tail"
  } else if (is.na(mcse) || (is.na(rel) && mcse > 0)) {
    "insufficient_info"
  } else if (!is.na(rel) && rel > 0.05) {
    "imprecise"
  } else {
    "ok"
  }
  row <- data.frame(disease = g$disease[1], method = g$method[1],
                    mean = mu, sd = s, mcse = mcse,
                    t(stats::setNames(q, paste0("q", probs))),
                    trimmed_mean = tmean, rel_mcse = rel, tail_share = tail_share,
                    stability = stability, check.names = FALSE, stringsAsFactors = FALSE)
  list(row = row, info = data.frame(pole = pole, n = length(x), stringsAsFactors = FALSE))
}

# Summary of every analysis of a batch run.
mc_batch_summary <- function(object) {
  parts <- lapply(names(object$analyses), function(nm) {
    s <- tryCatch(mc_summary(object$analyses[[nm]]), deconflate_error = function(e) NULL)
    if (is.null(s) || !nrow(s)) return(NULL)
    cbind(analysis = nm, s, stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(data.frame(analysis = character(0), stringsAsFactors = FALSE))
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}
