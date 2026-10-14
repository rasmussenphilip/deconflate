# Diagnosing the inputs: deconflate(diagnose = TRUE) ---------------------
#
# Which inputs drive problematic adjusted impacts, and which rows disagree
# with other rows for the same quantity. Design: section 14 of the v0.4
# design note; reference values: inst/validation/reference_v05_input_check.py.
#
# A disease's adjusted impact a is flagged when it lies on the other side of
# zero from its raw impact r (for event impacts, of 1), or has shrunk by more
# than 90%: c = sign(r) a - 0.1 |r| < 0, on the impact's scale (additive:
# units; event: log adjusted HR, and the row on its own scale; for risk
# difference rows only the sign counts, c = sign(r) a). r is fixed at the
# inputs as entered, so the flag disappears where a reaches 0.1 r (or 0).

check_grid <- c(0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.6, 2.2, 3)
check_offsets <- c(-0.5, -0.3, -0.15, 0, 0.15, 0.3, 0.5)
check_h <- 0.01
check_max_sd <- 3
check_top <- 4L

# Raw reference of each disease term with impact rows: its row, or with
# several rows their inverse-variance weighted mean on the impact's scale
# (an exact row's value when one is exact). For event impacts, `measure` is
# the measure the adjusted result is compared on: the rows' measure when
# they share one (HR, rate ratio, RR, OR or RD: the adjusted HR, or the
# risk ratio, odds ratio or risk difference of an animal with only the
# disease against a disease-free animal); with mixed measures, the hazard
# ratio and the ratio rows (on the log scale).
raw_reference <- function(model, event) {
  specs <- impact_specs(model)
  ids <- model$diseases$disease
  term <- vapply(specs, `[[`, character(1), "term")
  out <- list()
  for (d in intersect(ids, term)) {
    ss <- specs[term == d]
    meas <- vapply(ss, function(s) s$measure %||% NA_character_, character(1))
    on <- if (!event) "units" else if (length(unique(meas)) == 1L) {
      if (meas[1] == "rate_ratio") "HR" else meas[1]
    } else "HR"
    if (event && length(unique(meas)) > 1L) {
      ss <- ss[meas != "RD"]
      meas <- meas[meas != "RD"]
      if (!length(ss)) next
    }
    scale <- function(v, m) if (event) event_scale(v, m) else v
    v <- vapply(seq_along(ss), function(k) scale(ss[[k]]$value, meas[k]), numeric(1))
    se <- vapply(seq_along(ss), function(k) {
      row_se(row_dist(model, "impacts", ss[[k]]$row), function(x) scale(x, meas[k]))
    }, numeric(1))
    se[is.na(se)] <- 1
    r <- if (any(se == 0)) v[which(se == 0)[1]] else sum(v / se^2) / sum(1 / se^2)
    out[[d]] <- list(r = r, on = on, measure = if (length(ss) == 1L) meas[1] else "pooled",
                     rows = vapply(ss, `[[`, numeric(1), "row"))
  }
  out
}

# The adjusted result of each disease in `ref` on the measure it is compared
# on (natural scale: units, HR, RR, OR or RD), from a result's values.
compared_values <- function(vals, ref) {
  vapply(names(ref), function(d) {
    switch(ref[[d]]$on, units = , HR = vals$adj[[d]], RR = vals$rr[[d]], OR = vals$or[[d]], RD = vals$rd[[d]])
  }, numeric(1))
}

# The criterion c for each disease in `ref`, from its compared value (see
# compared_values()).
check_criterion <- function(cmp, ref, event) {
  vapply(names(ref), function(d) {
    rr <- ref[[d]]
    a <- if (event && rr$on != "RD") log(cmp[[d]]) else cmp[[d]]
    sign(rr$r) * a - 0.1 * abs(rr$r)
  }, numeric(1))
}

# Current value of an input (key as in the draws; "risk" for the overall risk).
input_value <- function(model, key, risk) {
  if (key == "risk") return(risk)
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  switch(parts[1],
    prob = model$diseases$value[match(parts[2], model$diseases$disease)],
    associations = model$associations$value[match(as.integer(parts[2]), model$associations$row)],
    impacts = model$impacts$value[match(as.integer(parts[2]), model$impacts$row)])
}

# Measure of an input (to choose its search range).
input_measure <- function(model, key) {
  if (key == "risk") return("risk")
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  switch(parts[1],
    prob = model$diseases$measure[match(parts[2], model$diseases$disease)],
    associations = model$associations$measure[match(as.integer(parts[2]), model$associations$row)],
    impacts = {
      m <- model$impacts$measure
      if (is.null(m)) NA_character_ else m[match(as.integer(parts[2]), model$impacts$row)]
    })
}

# Short labels of inputs, e.g. "SCM prevalence", "CM-SCM OR", "MAS impact
# (adjusted for SCK)"; the row number is added where labels repeat.
input_labels <- function(model, keys, risk_kind = "overall") {
  uk <- unique(keys)
  if (length(uk) < length(keys)) return(input_labels(model, uk, risk_kind)[match(keys, uk)])
  lab <- vapply(keys, function(k) {
    if (k == "risk") return(risk_label(risk_kind))
    parts <- strsplit(k, ":", fixed = TRUE)[[1]]
    qual <- function(among, adj) {
      q <- c(if (!is.na(among) && nzchar(among)) paste("among", among),
             if (!is.na(adj) && nzchar(adj)) paste("adjusted for", adj))
      if (length(q)) sprintf(" (%s)", paste(q, collapse = ", ")) else ""
    }
    switch(parts[1],
      prob = {
        r <- match(parts[2], model$diseases$disease)
        paste(parts[2], gsub("_", " ", model$diseases$measure[r]))
      },
      associations = {
        a <- model$associations
        r <- match(as.integer(parts[2]), a$row)
        paste0(a$disease[r], "-", a$associated_with[r], " ", a$measure[r], qual(a$among[r], a$adjusted_for[r]))
      },
      impacts = {
        im <- model$impacts
        r <- match(as.integer(parts[2]), im$row)
        m <- if (is.null(im$measure) || is.na(im$measure[r])) "impact" else im$measure[r]
        paste0(im$disease[r], " ", m, qual(im$among[r], im$adjusted_for[r]))
      })
  }, character(1), USE.NAMES = FALSE)
  dup <- lab %in% lab[duplicated(lab)]
  rows <- sub("^[a-z]+:", "", keys)
  lab[dup] <- sprintf("%s [row %s]", lab[dup], rows[dup])
  lab
}

# A model without one row (for the disagreeing rows).
drop_row <- function(model, key) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  tab <- parts[1]
  r <- as.integer(parts[2])
  d <- model[[tab]]
  keep <- d$row != r
  model[[tab]] <- structure(d[keep, , drop = FALSE], kind = attr(d, "kind"), label = attr(d, "label"),
                            units = attr(d, "units"))
  model$distributions[[key]] <- NULL
  model
}

# Adjusted values (natural scale) and the aggregate of a result; for event
# results also the risk ratio, odds ratio and risk difference of an animal
# with only the disease against a disease-free animal.
check_values <- function(r) {
  a <- r$adjusted
  if (inherits(r, "cm_event_result")) {
    return(list(adj = stats::setNames(a$adjusted, a$disease), rr = stats::setNames(a$rr_alone, a$disease),
                or = stats::setNames(a$or_alone, a$disease), rd = stats::setNames(a$rd_alone, a$disease),
                agg = r$attributable$summary$attributable %||% NA_real_))
  }
  list(adj = stats::setNames(a$adjusted, a$disease), agg = r$totals$adjusted_total)
}

# The diagnosis of a central result.
diagnose_inputs <- function(model, plan, res, risk, dots) {
  event <- inherits(res, "cm_event_result")
  joint_start <- if (!is.null(res$joint) && identical(res$joint$backend, "exact")) res$joint
  beta_start <- if (event) log(res$adjusted$adjusted)
  frozen <- list(associations = res$pooling$frozen, impacts = res$row_se)
  sampled <- identical(plan$fit_args$backend %||% dots$backend, "sampled")
  if (sampled && is.null(plan$fit_args$seed) && is.null(dots$seed)) plan$fit_args$seed <- 1L
  n_runs <- 0L
  last_error <- NULL
  # One adjustment with some inputs changed (named by key); NULL if it fails
  # (the condition is kept in `last_error`).
  run_with <- function(values, m0 = model, fz = frozen) {
    n_runs <<- n_runs + 1L
    tryCatch(withCallingHandlers({
      m <- set_inputs(m0, values[names(values) != "risk"])
      rk <- if ("risk" %in% names(values)) values[["risk"]] else risk$value
      if (!is.null(rk) && !(rk > 0 && rk < 1)) {
        cm_abort(sprintf("The %s must lie between 0 and 1.", risk_label(risk$kind)), class = "deconflate_infeasible")
      }
      only_impacts <- length(values) > 0L && all(startsWith(names(values), "impacts:"))
      jt <- if (only_impacts && identical(m0, model) && !isTRUE(plan$shortcut)) res$joint
      check_values(run_point(m, plan, risk = rk, joint = jt, warn = FALSE, dots = dots,
                             joint_start = joint_start, start = beta_start, frozen = fz, allocate = FALSE,
                             raw = FALSE))
    }, deconflate_nonconvergence = function(w) {
      if (inherits(w, "warning")) invokeRestart("muffleWarning")
    }), deconflate_error = function(e) {
      last_error <<- e
      NULL
    })
  }

  central <- check_values(res)
  ref <- raw_reference(model, event)
  cmp0 <- compared_values(central, ref)
  c0 <- check_criterion(cmp0, ref, event)
  flag <- names(c0)[c0 < 0]
  on_log <- vapply(ref, function(rr) event && rr$on != "RD", logical(1))
  raw_nat <- vapply(names(ref), function(d) if (on_log[[d]]) exp(ref[[d]]$r) else ref[[d]]$r, numeric(1))
  side <- sign(ifelse(on_log, log(cmp0), cmp0)) != sign(vapply(ref, `[[`, numeric(1), "r"))
  flagged <- data.frame(disease = names(ref), measure = vapply(ref, `[[`, character(1), "on"),
                        raw = unname(raw_nat), adjusted = unname(cmp0),
                        criterion = unname(c0), flagged = c0 < 0,
                        reason = ifelse(c0 >= 0, "", ifelse(side, ifelse(on_log, "other side of 1", "other side of zero"),
                                                            "shrunk by more than 90%")),
                        stringsAsFactors = FALSE)
  if (!event) flagged$measure <- NULL
  rownames(flagged) <- NULL
  crit_of <- function(out) {
    if (is.null(out)) rep(NA_real_, length(flag)) else check_criterion(compared_values(out, ref), ref, event)[flag]
  }

  # Inputs with a distribution (keys in their canonical form).
  dists <- model$distributions %||% list()
  if (length(dists)) {
    dists <- dists[!startsWith(names(dists), "between:")]
    keys <- vapply(names(dists), function(k) {
      if (grepl("^(prob|associations|impacts):", k)) k else
        tryCatch(canonical_key(k, model), deconflate_error = function(e) NA_character_)
    }, character(1), USE.NAMES = FALSE)
    names(dists) <- keys
    dists <- dists[!is.na(keys)]
  }
  if (!is.null(risk$dist)) dists[["risk"]] <- risk$dist
  uncertain <- names(dists)[!vapply(dists, function(d) identical(d$type, "fixed"), logical(1))]
  effects <- NULL
  if (length(flag) && length(uncertain)) {
    rows <- lapply(uncertain, function(k) {
      d <- dists[[k]]
      x0 <- input_value(model, k, risk$value)
      xs <- tryCatch({
        u0 <- stats::qnorm(d$p(x0))
        u0 <- if (is.finite(u0)) min(max(u0, -5), 5) else if (is.na(u0)) 0 else sign(u0) * 5
        c(d$q(stats::pnorm(u0 + check_h)), d$q(stats::pnorm(u0 - check_h)))
      }, error = function(e) c(NA_real_, NA_real_))
      op <- if (is.finite(xs[1])) run_with(stats::setNames(xs[1], k))
      om <- if (is.finite(xs[2])) run_with(stats::setNames(xs[2], k))
      eff <- if (is.null(op) || is.null(om)) rep(NA_real_, length(flag)) else
        unname((compared_values(op, ref)[flag] - compared_values(om, ref)[flag]) / (2 * check_h))
      dc <- (crit_of(op) - crit_of(om)) / (2 * check_h)
      need <- ifelse(is.finite(dc) & abs(dc) > 1e-12, -c0[flag] / dc, NA_real_)
      data.frame(key = k, value = x0, disease = flag, effect_per_sd = eff, sd_needed = unname(need),
                 stringsAsFactors = FALSE)
    })
    effects <- do.call(rbind, rows)
    effects$input <- input_labels(model, effects$key, risk$kind)
    listed <- rep(FALSE, nrow(effects))
    for (d in flag) {
      i <- which(effects$disease == d & is.finite(effects$effect_per_sd))
      top <- i[order(-abs(effects$effect_per_sd[i]))][seq_len(min(check_top, length(i)))]
      listed[effects$key %in% effects$key[top]] <- TRUE
    }
    fixes <- is.finite(effects$sd_needed) & abs(effects$sd_needed) <= check_max_sd
    listed[effects$key %in% effects$key[fixes]] <- TRUE
    effects$listed <- listed
    effects <- effects[, c("key", "input", "value", "disease", "effect_per_sd", "sd_needed", "listed")]
    rownames(effects) <- NULL
  }

  # Inputs without a distribution: the value at which the flag disappears.
  thresholds <- NULL
  if (length(flag)) {
    exact <- c(paste0("prob:", model$diseases$disease),
               if (!is.null(model$associations)) paste0("associations:", model$associations$row),
               paste0("impacts:", model$impacts$row),
               if (event) "risk")
    exact <- setdiff(exact, uncertain)
    out <- list()
    add <- function(k, v, j, x, in_range) {
      out[[length(out) + 1L]] <<- data.frame(key = k, value = v, disease = flag[j], threshold = unname(x),
                                             in_range = in_range, stringsAsFactors = FALSE)
    }
    for (k in exact) {
      v <- input_value(model, k, risk$value)
      if (is.null(v) || !is.finite(v)) next
      meas <- input_measure(model, k)
      lo <- min(0.25 * v, 3 * v)
      hi <- max(0.25 * v, 3 * v)
      if (!event && startsWith(k, "impacts:")) {
        # Additive impacts are linear in an impact row: solved exactly. Shown
        # for the flagged disease's own rows, and for other rows when within
        # a quarter to three times the value.
        own <- names(ref)[vapply(ref, function(rr) k %in% paste0("impacts:", rr$rows), logical(1))]
        o1 <- run_with(stats::setNames(v + 1, k))
        if (is.null(o1)) next
        slope <- crit_of(o1) - c0[flag]
        for (j in seq_along(flag)) {
          if (!is.finite(slope[j]) || abs(slope[j]) < 1e-12) next
          x <- v - c0[flag][j] / slope[j]
          ok <- x >= lo && x <= hi
          if (ok || flag[j] %in% own) add(k, v, j, x, ok)
        }
        next
      }
      rd <- identical(meas, "RD") || identical(meas, "phi")
      xs <- if (rd) v + check_offsets else v * check_grid
      if (rd) xs <- xs[xs > -1 & xs < 1]
      if (k == "risk" || meas %in% c("prevalence", "probability")) xs <- xs[xs > 0 & xs < 1]
      if (length(xs) < 2L) next
      cs <- lapply(xs, function(x) {
        if (x == v) return(c0[flag])
        crit_of(run_with(stats::setNames(x, k)))
      })
      for (j in seq_along(flag)) {
        roots <- numeric(0)
        for (i in seq_len(length(xs) - 1L)) {
          a <- cs[[i]][j]
          b <- cs[[i + 1L]][j]
          if (!is.finite(a) || !is.finite(b) || sign(a) == sign(b)) next
          f <- function(x) {
            cc <- crit_of(run_with(stats::setNames(x, k)))[j]
            if (is.finite(cc)) cc else stop("failed run")
          }
          rt <- tryCatch(stats::uniroot(f, c(xs[i], xs[i + 1L]), f.lower = a, f.upper = b, tol = 1e-10)$root,
                         error = function(e) NA_real_)
          if (is.finite(rt)) roots <- c(roots, rt)
        }
        if (length(roots)) {
          dist <- abs(roots - v) / (if (v != 0) abs(v) else 1)
          add(k, v, j, roots[which.min(dist)], TRUE)
        }
      }
    }
    if (length(out)) {
      thresholds <- do.call(rbind, out)
      thresholds$input <- input_labels(model, thresholds$key, risk$kind)
      thresholds <- thresholds[, c("key", "input", "value", "disease", "threshold", "in_range")]
      rownames(thresholds) <- NULL
    }
  }
  fixed_by <- character(0)
  if (!is.null(effects)) {
    fixed_by <- effects$disease[is.finite(effects$sd_needed) & abs(effects$sd_needed) <= check_max_sd]
  }
  if (!is.null(thresholds)) fixed_by <- c(fixed_by, thresholds$disease[thresholds$in_range])
  no_single_fix <- setdiff(flag, fixed_by)

  # Rows that disagree with other rows for the same quantity.
  disagreeing <- NULL
  est <- res$estimates
  bad <- which(is.finite(est$std_residual) & abs(est$std_residual) >= 2)
  if (length(bad)) {
    out <- lapply(bad, function(i) {
      key <- paste0(est$table[i], ":", est$row[i])
      last_error <<- NULL
      o <- run_with(list(), m0 = drop_row(model, key), fz = NULL)
      base <- data.frame(table = est$table[i], row = est$row[i], item = est$item[i], study = est$study[i],
                         value = est$value[i], sd = est$sd[i], fitted = est$fitted[i],
                         std_residual = est$std_residual[i], stringsAsFactors = FALSE)
      common <- if (is.null(o)) character(0) else intersect(names(central$adj), names(o$adj))
      dlt <- abs(o$adj[common] - central$adj[common])
      if (is.null(o) || !any(is.finite(dlt))) {
        note <- if (inherits(last_error, "deconflate_singular")) {
          "without it, the other rows do not determine the adjusted impacts"
        } else if (!is.null(last_error)) {
          sprintf("without it, the model could not be adjusted (%s)", condition_type(last_error))
        } else {
          "without it, the adjusted impacts cannot be compared"
        }
        return(cbind(base, changes_most = NA_character_, before = NA_real_, after = NA_real_,
                     aggregate_before = central$agg, aggregate_after = NA_real_, note = note,
                     stringsAsFactors = FALSE))
      }
      t <- common[which.max(dlt)]
      cbind(base, changes_most = t, before = unname(central$adj[t]), after = unname(o$adj[t]),
            aggregate_before = central$agg, aggregate_after = o$agg, note = "", stringsAsFactors = FALSE)
    })
    disagreeing <- do.call(rbind, out)
    rownames(disagreeing) <- NULL
  }

  structure(list(flagged = flagged, inputs = effects, thresholds = thresholds,
                 no_single_fix = no_single_fix, disagreeing = disagreeing, event = event,
                 settings = list(h = check_h, max_sd = check_max_sd, top = check_top, grid = check_grid,
                                 offsets = check_offsets, sampled = sampled, runs = n_runs)),
            class = "cm_diagnosis")
}

#' @export
print.cm_diagnosis <- function(x, digits = 3, ...) {
  f <- x$flagged
  fl <- f[f$flagged, , drop = FALSE]
  fmt <- function(v) {
    out <- formatC(signif(v, digits), digits = digits, format = "fg", flag = "#")
    sub("\\.$", "", trimws(out))
  }
  if (!nrow(fl) && is.null(x$disagreeing)) {
    cat("Diagnosis: no adjusted impact is flagged and no row disagrees with other rows.\n")
    return(invisible(x))
  }
  cat("Diagnosis (diagnose = TRUE)\n")
  if (nrow(fl)) {
    if (x$event) {
      cat("\nFlagged adjusted results (raw estimate -> adjusted, on the estimate's measure; RR, OR and RD\nfor an animal with only the disease against a disease-free animal):\n")
      for (i in seq_len(nrow(fl))) {
        cat(sprintf("  %s: %s %s -> %s (%s)\n", fl$disease[i], fl$measure[i], fmt(fl$raw[i]), fmt(fl$adjusted[i]),
                    fl$reason[i]))
      }
    } else {
      cat("\nFlagged adjusted impacts (raw -> adjusted):\n")
      for (i in seq_len(nrow(fl))) {
        cat(sprintf("  %s: %s -> %s (%s)\n", fl$disease[i], fmt(fl$raw[i]), fmt(fl$adjusted[i]), fl$reason[i]))
      }
    }
    flag <- fl$disease
    max_sd <- x$settings$max_sd
    e <- x$inputs
    if (!is.null(e) && any(e$listed)) {
      cat(sprintf("\nInputs with a distribution: change in each flagged %s when the input moves up by 1 SD;\nin brackets, the shift that alone removes the flag (within %g SD).\n",
                  if (x$event) "result (on the measure above)" else "impact", max_sd))
      keys <- unique(e$key[e$listed])
      best <- vapply(keys, function(k) {
        n <- abs(e$sd_needed[e$key == k])
        n <- n[is.finite(n) & n <= max_sd]
        if (length(n)) min(n) else Inf
      }, numeric(1))
      size <- vapply(keys, function(k) {
        z <- abs(e$effect_per_sd[e$key == k])
        if (any(is.finite(z))) max(z, na.rm = TRUE) else 0
      }, numeric(1))
      keys <- keys[order(best, -size)]
      tab <- data.frame(input = vapply(keys, function(k) e$input[e$key == k][1], character(1)),
                        value = vapply(keys, function(k) fmt(e$value[e$key == k][1]), character(1)),
                        stringsAsFactors = FALSE)
      for (d in flag) {
        z <- abs(e$effect_per_sd[e$disease == d])
        big <- if (any(is.finite(z))) max(z, na.rm = TRUE) else 0
        tab[[d]] <- vapply(keys, function(k) {
          r <- e[e$key == k & e$disease == d, , drop = FALSE]
          if (!nrow(r) || !is.finite(r$effect_per_sd)) return(".")
          n <- abs(r$sd_needed)
          fixes <- is.finite(n) && n <= max_sd
          if (!fixes && abs(r$effect_per_sd) < 0.02 * big) return(".")
          s <- if (abs(r$effect_per_sd) >= 0.1) sprintf("%+.2f", r$effect_per_sd) else sprintf("%+.2g", r$effect_per_sd)
          if (fixes) {
            s <- sprintf("%s (%s SD %s)", s, if (n >= 0.095) sprintf("%.1f", n) else sprintf("%.2g", n),
                         if (r$sd_needed > 0) "up" else "down")
          }
          s
        }, character(1))
      }
      print(tab, row.names = FALSE, right = FALSE)
    }
    th <- x$thresholds
    if (!is.null(th) && nrow(th)) {
      cat("\nInputs without a distribution: the value at which the flag disappears.\n")
      keys <- unique(th$key)
      tab <- data.frame(input = vapply(keys, function(k) th$input[th$key == k][1], character(1)),
                        value = vapply(keys, function(k) fmt(th$value[th$key == k][1]), character(1)),
                        stringsAsFactors = FALSE)
      for (d in flag) {
        tab[[d]] <- vapply(keys, function(k) {
          r <- th[th$key == k & th$disease == d, , drop = FALSE]
          if (!nrow(r)) return(".")
          paste0(fmt(r$threshold), if (!r$in_range) "*" else "")
        }, character(1))
      }
      print(tab, row.names = FALSE, right = FALSE)
      if (!all(th$in_range)) {
        cat("* outside a quarter to three times the value (a disease's own impact rows are always shown).\n")
      }
    }
    if (length(x$no_single_fix)) {
      cat(sprintf("\nNo single input removes the flag for %s (within %g SD, or within a quarter to three times\nthe value of an input without a distribution): several inputs together, or a structural problem such as units.\n",
                  join_and(x$no_single_fix), max_sd))
    }
    if (isTRUE(x$settings$sampled)) {
      cat("\nThe joint distribution was fitted with the sampled backend: the effects include Monte Carlo noise.\n")
    }
  }
  dg <- x$disagreeing
  if (!is.null(dg)) {
    cat("\nRows that disagree with other rows for the same quantity (standardised residual of 2 or more):\n")
    for (i in seq_len(nrow(dg))) {
      who <- paste(c(dg$item[i], if (!is.na(dg$study[i]) && nzchar(dg$study[i])) dg$study[i]), collapse = ", ")
      cons <- if (is.na(dg$after[i])) dg$note[i] else
        sprintf("%s %s -> %s, aggregate %s -> %s", dg$changes_most[i], fmt(dg$before[i]), fmt(dg$after[i]),
                fmt(dg$aggregate_before[i]), fmt(dg$aggregate_after[i]))
      cat(sprintf("  %s row %s (%s): %s, fitted %s, residual %+.1f; without it: %s\n", dg$table[i], dg$row[i], who,
                  fmt(dg$value[i]), fmt(dg$fitted[i]), dg$std_residual[i], cons))
    }
  }
  invisible(x)
}
