#' Fit the maximum-entropy distribution of disease combinations
#'
#' The global model needs the probability of every combination of diseases.
#' Pairwise tables do not determine it (for three or more diseases there are
#' 2^n - 1 free probabilities and n(n+1)/2 constraints), so the package uses
#' the maximum-entropy distribution that matches the disease probabilities
#' and every specified pairwise table: a log-linear model with main effects
#' and the two-way terms of the constrained pairs, plus any three-way terms
#' from [cm_three_way()]. This is an assumption, not something the pairwise
#' evidence identifies. Pairs without an association (unknown pairs) are not
#' constrained: their association is whatever the maximum-entropy fit implies
#' given the others. A pair given an odds ratio of 1 is constrained to be
#' independent.
#'
#' @section Exact backend:
#' `backend = "exact"` enumerates all 2^n combinations. Starting from
#' independence (the product of marginal probabilities; slide 12 of the 2023
#' seminar), with the three-way terms put into the starting distribution,
#' iterative proportional fitting (IPF) rescales the combination
#' probabilities until every pairwise table and marginal is matched. IPF
#' keeps every log-linear term of its start that the constraints do not fix,
#' so the result is the maximum-entropy distribution with the requested
#' three-way terms. Enumeration limits this backend to about 20 diseases.
#'
#' @section Sampled backend:
#' `backend = "sampled"` fits the same log-linear model without enumerating
#' combinations, for larger numbers of diseases:
#' 1. Calibration: the model's parameters (a main effect per disease and an
#'    interaction per constrained pair; three-way terms fixed at the log of
#'    their ratios) are fitted by Monte Carlo moment matching. Gibbs sampling
#'    from many parallel chains estimates the probabilities and pairwise
#'    tables implied by the current parameters, and the parameters are
#'    updated (by the differences in logits and log odds ratios between the
#'    targets and the estimates) until the estimates match the targets
#'    within Monte Carlo error. The interaction parameters are conditional
#'    log odds ratios; they are not set equal to the marginal log odds ratios
#'    of the inputs (those are only the starting values).
#' 2. Sampling: `n_samples` combinations are drawn from the calibrated model.
#' 3. Optionally (`calibrate = TRUE`, the default), the weights of the sampled
#'    combinations are adjusted by IPF so that the disease probabilities and
#'    pairwise tables match the targets exactly on the sample (raking). The
#'    remaining Monte Carlo error then affects only higher-order structure.
#'
#' The result has the same form as the exact one (`cells` are the distinct
#' sampled combinations, `prob` their weights), so the global method,
#' interaction offsets, Shapley allocation and the snapshot hazard model of
#' event impacts work with it. Its `diagnostics` give the constraint
#' residuals of the raw sample and after raking, R-hat and effective sample
#' sizes of the fitted probabilities across chains, and the Monte Carlo
#' errors. Results are approximations with Monte Carlo error: compare runs
#' with different seeds, and with the exact backend where n allows.
#'
#' @section Convergence:
#' If the pairwise tables are jointly infeasible, the exact IPF does not
#' converge. Non-convergence after `max_iter` sweeps suggests, but does not
#' prove, infeasibility (see [check_feasibility()]). For the sampled backend,
#' a calibration that does not converge is reported as unresolved, not as
#' evidence that the inputs are infeasible.
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param tol Convergence tolerance of the IPF: maximum absolute deviation of
#'   any constrained cell probability, recomputed on the returned
#'   distribution.
#' @param max_iter Maximum number of full IPF sweeps.
#' @param max_diseases Safety limit on the number of diseases for the exact
#'   backend (the table has 2^n cells).
#' @param backend `"exact"` (default) or `"sampled"`.
#' @param n_samples Sampled backend: number of combinations drawn after
#'   calibration.
#' @param n_chains Sampled backend: number of parallel Gibbs chains.
#' @param burn_in Sampled backend: Gibbs sweeps discarded before calibration
#'   and before the final draws.
#' @param fit_iter Sampled backend: maximum number of moment-matching
#'   iterations.
#' @param calibrate Sampled backend: adjust the sample weights to match the
#'   targets exactly?
#' @param seed Optional random seed (sampled backend).
#' @param start Optional earlier exact fit to start the IPF from (e.g. the fit
#'   at the central inputs, when refitting for a draw of the inputs). It is
#'   used only when it has the same diseases, constrained pairs and three-way
#'   terms; the result is then the same maximum-entropy distribution, found
#'   in fewer sweeps.
#'
#' @return A `cm_joint` object: `cells` (0/1 matrix of combinations), `prob`
#'   (probability of each combination), `converged`, `iterations`,
#'   `max_residual`, `backend`, the `targets` it was fitted to (used by
#'   [deconflate()] to check that a supplied joint matches the model) and,
#'   for the sampled backend, `diagnostics`.
#' @export
#' @examples
#' j <- fit_joint(example_supplement())
#' j
#' \donttest{
#' js <- fit_joint(example_supplement(), backend = "sampled", n_samples = 20000, seed = 1)
#' js$diagnostics$summary
#' }
fit_joint <- function(model, tol = 1e-10, max_iter = 10000L, max_diseases = 20L,
                      backend = c("exact", "sampled"), n_samples = 50000L, n_chains = 1000L,
                      burn_in = 50L, fit_iter = 300L, calibrate = TRUE, seed = NULL,
                      start = NULL) {
  check_population(model)
  backend <- match.arg(backend)
  ids <- model$diseases$id
  p <- model$diseases$prob
  n <- length(ids)
  tw <- model$three_way
  pt <- pair_tables(model)
  cons <- pt[!is.na(pt$p11), , drop = FALSE]
  targets <- list(prob = p, pairs = cons[, c("disease1", "disease2", "p11")],
                  three_way = if (is.null(tw)) NULL else as.data.frame(tw))
  if (backend == "sampled") {
    return(fit_joint_sampled(model, ids, p, tw, cons, targets, tol = tol, max_iter = max_iter,
                             n_samples = n_samples, n_chains = n_chains, burn_in = burn_in,
                             fit_iter = fit_iter, calibrate = calibrate, seed = seed))
  }
  if (n > max_diseases) {
    cm_abort(sprintf("%d diseases give 2^%d cells; use backend = \"sampled\", or increase `max_diseases` if intended.", n, n))
  }
  cells <- disease_cells(n)
  colnames(cells) <- ids
  logp <- as.vector(cells %*% log(p) + (1 - cells) %*% log(1 - p))
  if (!is.null(tw) && nrow(tw)) {
    for (r in seq_len(nrow(tw))) {
      k <- match(c(tw$disease1[r], tw$disease2[r], tw$disease3[r]), ids)
      logp <- logp + log(tw$ratio[r]) * cells[, k[1]] * cells[, k[2]] * cells[, k[3]]
    }
  }
  prob <- exp(logp - max(logp))
  prob <- prob / sum(prob)
  if (warm_start_ok(start, ids, cons, tw)) prob <- start$prob
  fit <- ipf_cells(cells, prob, p, cons, ids, tol, max_iter)
  if (!fit$converged) {
    cm_warn(sprintf(
      "IPF did not converge after %d sweeps (max residual %.2e). The pairwise tables may be jointly infeasible; see check_feasibility().",
      fit$iterations, fit$residual), class = "deconflate_nonconvergence")
  }
  structure(list(cells = cells, prob = fit$prob, diseases = ids,
                 converged = fit$converged, iterations = fit$iterations, max_residual = fit$residual,
                 constrained_pairs = cons[, c("disease1", "disease2", "measure", "status")],
                 targets = targets, backend = "exact"),
            class = "cm_joint")
}

# Arguments of fit_joint() that may be passed through `...`.
joint_arg_names <- c("tol", "max_iter", "max_diseases", "backend", "n_samples", "n_chains",
                     "burn_in", "fit_iter", "calibrate", "seed")

# Can an earlier exact fit be the start of the IPF? IPF keeps the log-linear
# terms of its start that the constraints do not fix, so the start must be a
# converged fit of the same family: the same diseases, the same constrained
# pairs and the same three-way terms (only the targets may differ).
warm_start_ok <- function(start, ids, cons, tw) {
  if (is.null(start)) return(FALSE)
  if (!inherits(start, "cm_joint") || !identical(start$backend, "exact") || !isTRUE(start$converged) ||
      !identical(start$diseases, ids) || length(start$prob) != 2^length(ids) || any(start$prob <= 0)) {
    return(FALSE)
  }
  tg <- start$targets$pairs
  if (is.null(tg) || !setequal(pair_key(tg$disease1, tg$disease2), pair_key(cons$disease1, cons$disease2))) {
    return(FALSE)
  }
  sig <- function(t) {
    if (is.null(t) || !nrow(t)) return(character(0))
    sort(paste(apply(t[, c("disease1", "disease2", "disease3")], 1,
                     function(x) paste(sort(x), collapse = "|")), signif(t$ratio, 12)))
  }
  identical(sig(tw), sig(start$targets$three_way))
}

# Iterative proportional fitting of the weights `prob` of the combinations
# `cells` to the marginals `p` and the constrained pairwise tables `cons`.
ipf_cells <- function(cells, prob, p, cons, ids, tol, max_iter) {
  n <- length(ids)
  idx <- match(c(cons$disease1, cons$disease2), ids)
  i1 <- idx[seq_len(nrow(cons))]
  i2 <- idx[nrow(cons) + seq_len(nrow(cons))]
  # Group index 1..4 = (0,0), (0,1), (1,0), (1,1) for (d1, d2).
  groups <- lapply(seq_len(nrow(cons)), function(r) as.integer(2 * cells[, i1[r]] + cells[, i2[r]] + 1))
  targets <- lapply(seq_len(nrow(cons)), function(r) {
    p1 <- p[i1[r]]
    p2 <- p[i2[r]]
    p11 <- cons$p11[r]
    c(1 - p1 - p2 + p11, p2 - p11, p1 - p11, p11)
  })
  # Probability of each of the four cells of a pair (rowsum() works in C).
  cell_sums <- function(pr, g) {
    out <- numeric(4)
    rs <- rowsum(pr, g)
    out[as.integer(rownames(rs))] <- rs[, 1]
    out
  }
  margin_resid <- function(pr) max(abs(colSums(cells * pr) - p))
  pair_resid <- function(pr) {
    if (!length(groups)) return(0)
    max(vapply(seq_along(groups), function(r) max(abs(cell_sums(pr, groups[[r]]) - targets[[r]])),
               numeric(1)))
  }
  # A positive target with no combination to carry it cannot be matched
  # (possible on a sample, never with all 2^n combinations): do not iterate.
  has <- function(sel) any(prob[sel] > 0)
  empty <- nrow(cells) < 2^n && (
    any(vapply(seq_len(n), function(i) {
      (p[i] > tol && !has(cells[, i] == 1)) || (1 - p[i] > tol && !has(cells[, i] == 0))
    }, logical(1))) ||
      any(vapply(seq_along(groups), function(r) any(targets[[r]] > tol & cell_sums(prob, groups[[r]]) <= 0),
                 logical(1))))
  if (empty) max_iter <- 0L
  iter <- 0L
  for (it in seq_len(max_iter)) {
    iter <- it
    for (r in seq_along(groups)) {
      g <- groups[[r]]
      cur <- cell_sums(prob, g)
      fac <- ifelse(cur > 0, targets[[r]] / cur, 0)
      prob <- prob * fac[g]
    }
    # Marginals of every disease (needed for diseases without constrained
    # pairs, e.g. when all their pairs are unknown).
    for (i in seq_len(n)) {
      cur <- sum(prob[cells[, i] == 1])
      if (cur > 0 && cur < 1) {
        prob <- prob * ifelse(cells[, i] == 1, p[i] / cur, (1 - p[i]) / (1 - cur))
      }
    }
    if (max(pair_resid(prob), margin_resid(prob)) < tol) break
  }
  resid <- max(pair_resid(prob), margin_resid(prob))
  list(prob = prob, iterations = iter, residual = resid, converged = is.finite(resid) && resid < tol)
}

# One Gibbs sweep over all diseases for every chain (rows of X), for the
# log-linear model with main effects h, pairwise terms J (symmetric, zero
# diagonal) and fixed three-way terms (rows of `tri`: three column indices and
# a log ratio).
gibbs_sweep <- function(X, h, J, tri) {
  n <- ncol(X)
  for (i in seq_len(n)) {
    eta <- h[i] + as.vector(X %*% J[, i])
    if (nrow(tri)) {
      for (r in which(tri[, 1] == i | tri[, 2] == i | tri[, 3] == i)) {
        others <- setdiff(tri[r, 1:3], i)
        eta <- eta + tri[r, 4] * X[, others[1]] * X[, others[2]]
      }
    }
    X[, i] <- as.numeric(stats::runif(nrow(X)) < stats::plogis(eta))
  }
  X
}

fit_joint_sampled <- function(model, ids, p, tw, cons, targets, tol, max_iter, n_samples,
                              n_chains, burn_in, fit_iter, calibrate, seed) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(ids)
  n_chains <- as.integer(n_chains)
  if (n_chains < 10L) cm_abort("`n_chains` must be at least 10.")
  i1 <- match(cons$disease1, ids)
  i2 <- match(cons$disease2, ids)
  p11 <- cons$p11
  logit <- function(x) stats::qlogis(pmin(pmax(x, 1e-12), 1 - 1e-12))
  log_or <- function(a, b, ab) {
    t11 <- pmax(ab, 1e-12); t10 <- pmax(a - ab, 1e-12); t01 <- pmax(b - ab, 1e-12)
    t00 <- pmax(1 - a - b + ab, 1e-12)
    log(t11) + log(t00) - log(t10) - log(t01)
  }
  tri <- matrix(numeric(0), 0, 4)
  if (!is.null(tw) && nrow(tw)) {
    tri <- cbind(match(tw$disease1, ids), match(tw$disease2, ids), match(tw$disease3, ids),
                 log(tw$ratio))
  }
  # Start: independence plus the marginal log odds ratios of the targets
  # (starting values only; they are calibrated below).
  h <- logit(p)
  J <- matrix(0, n, n)
  if (length(i1)) {
    lor <- pmin(pmax(log_or(p[i1], p[i2], p11), -30), 30)
    J[cbind(i1, i2)] <- lor
    J[cbind(i2, i1)] <- lor
    for (r in seq_along(i1)) {
      h[i1[r]] <- h[i1[r]] - lor[r] * p[i2[r]]
      h[i2[r]] <- h[i2[r]] - lor[r] * p[i1[r]]
    }
  }
  X <- matrix(as.numeric(stats::runif(n_chains * n) < rep(p, each = n_chains)), n_chains, n)
  for (b in seq_len(burn_in)) X <- gibbs_sweep(X, h, J, tri)

  # Moment matching: damped updates of the logits and conditional log odds
  # ratios by the differences between targets and estimates. Each step runs
  # `sweeps` Gibbs sweeps on all chains; the Monte Carlo error of each
  # estimate comes from the spread of the chain means (the chains are
  # independent).
  tgt <- c(p, p11)
  n_mom <- length(tgt)
  z_k <- max(3, stats::qnorm(1 - 0.025 / n_mom))
  moments_of <- function(X) {
    if (length(i1)) cbind(X, X[, i1, drop = FALSE] * X[, i2, drop = FALSE]) else X
  }
  run_sweeps <- function(sweeps) {
    acc <- matrix(0, n_chains, n_mom)
    for (s in seq_len(sweeps)) {
      X <<- gibbs_sweep(X, h, J, tri)
      acc <- acc + moments_of(X)
    }
    cm <- acc / sweeps
    list(est = colMeans(cm), se = pmax(apply(cm, 2, stats::sd) / sqrt(n_chains), 1e-6))
  }
  update <- function(est, step) {
    ph <- est[seq_len(n)]
    if (length(i1)) {
      p11h <- est[n + seq_along(i1)]
      dJ <- pmin(pmax(log_or(p[i1], p[i2], p11) - log_or(ph[i1], ph[i2], p11h), -2), 2)
      J[cbind(i1, i2)] <<- pmin(pmax(J[cbind(i1, i2)] + step * dJ, -30), 30)
      J[cbind(i2, i1)] <<- J[cbind(i1, i2)]
    }
    h <<- h + step * pmin(pmax(logit(p) - logit(ph), -2), 2)
  }
  history <- numeric(0)
  fit_converged <- FALSE
  calm <- 0L
  it <- 0L
  for (it in seq_len(fit_iter)) {
    sweeps <- if (it <= 10L) 5L else if (it <= 20L) 10L else 20L
    e <- run_sweeps(sweeps)
    history <- c(history, max(abs(e$est - tgt)))
    # Converged when every estimate is within a few Monte Carlo errors of
    # its target (allowing for the number of moments and for the noise in
    # the parameters), three times in a row.
    zres <- max(abs(e$est - tgt) / e$se)
    calm <- if (it >= 15L && zres < 2 * z_k) calm + 1L else 0L
    if (calm >= 3L) {
      fit_converged <- TRUE
      break
    }
    update(e$est, if (it <= 20L) 0.8 else 0.5)
  }
  # Polyak averaging of the parameters over a final set of small steps
  # reduces their Monte Carlo noise.
  hs <- matrix(0, 10L, n)
  Js <- array(0, c(10L, n, n))
  for (k in seq_len(10L)) {
    e <- run_sweeps(20L)
    update(e$est, 0.3)
    hs[k, ] <- h
    Js[k, , ] <- J
  }
  h <- colMeans(hs)
  J <- colMeans(Js, dims = 1L)

  # Final draws: parallel chains after a further burn-in; each chain's mean
  # gives the between-chain error, R-hat and effective sample size.
  for (b in seq_len(burn_in)) X <- gibbs_sweep(X, h, J, tri)
  n_sweeps <- max(2L, ceiling(n_samples / n_chains))
  draws <- vector("list", n_sweeps)
  chain_sum <- matrix(0, n_chains, n_mom)
  for (s in seq_len(n_sweeps)) {
    X <- gibbs_sweep(X, h, J, tri)
    draws[[s]] <- X
    chain_sum <- chain_sum + moments_of(X)
  }
  Xall <- do.call(rbind, draws)
  colnames(Xall) <- ids
  target_all <- c(p, p11)
  # Chain-level diagnostics for each moment. Every moment is a 0/1
  # indicator, so its mean square equals its mean, and every chain has the
  # same number of draws, so the sample estimate is the mean of chain means.
  chain_means <- chain_sum / n_sweeps                               # chains x moments
  est_all <- colMeans(chain_means)
  within <- (chain_means - chain_means^2) * n_sweeps / (n_sweeps - 1)
  W <- colMeans(within)
  B <- n_sweeps * apply(chain_means, 2, stats::var)
  var_plus <- (n_sweeps - 1) / n_sweeps * W + B / n_sweeps
  # No variation within chains but differences between them: stuck chains.
  rhat <- ifelse(W > 0, sqrt(var_plus / W), ifelse(B > 0, Inf, 1))
  mcse <- apply(chain_means, 2, stats::sd) / sqrt(n_chains)
  total_var <- est_all * (1 - est_all)
  ess <- ifelse(mcse > 0, total_var / mcse^2, nrow(Xall))
  mom_names <- c(paste0("prob:", ids), if (length(i1)) paste0("pair:", ids[i1], ":", ids[i2]))
  moments <- data.frame(moment = mom_names, target = target_all, sampled = est_all,
                        residual = est_all - target_all, mcse = mcse, rhat = rhat, ess = ess,
                        stringsAsFactors = FALSE)

  # Distinct combinations and their frequencies.
  key <- if (n <= 50L) as.vector(Xall %*% 2^(seq_len(n) - 1)) else apply(Xall, 1, paste, collapse = "")
  ukey <- unique(key)
  pos <- match(key, ukey)
  cells <- Xall[match(ukey, key), , drop = FALSE]
  rownames(cells) <- NULL
  prob <- tabulate(pos, nbins = length(ukey)) / nrow(Xall)
  resid_raw <- max(abs(est_all - target_all))
  cal <- NULL
  if (calibrate) {
    # An unresolved calibration stays unresolved whatever the raking does, so
    # its raking (slow for many cells and pairs) is only a short diagnostic.
    cal <- ipf_cells(cells, prob, p, cons, ids, tol = max(tol, 1e-10),
                     max_iter = if (fit_converged) max_iter else min(max_iter, 100L))
    prob <- cal$prob
  }
  resid_final <- if (calibrate) cal$residual else resid_raw
  # Converged: the calibration matched the targets within Monte Carlo error
  # and, if requested, the raking matched them exactly.
  mc_ok <- max(abs(est_all - target_all) / pmax(mcse, 1e-6)) <= 2 * z_k || resid_raw <= 2e-3
  converged <- fit_converged && mc_ok && (!calibrate || isTRUE(cal$converged))
  diag <- list(
    summary = data.frame(n_samples = nrow(Xall), n_chains = n_chains, n_unique = nrow(cells),
                         fit_iterations = it, fit_converged = fit_converged,
                         residual_sample = resid_raw, residual_calibrated = if (calibrate) cal$residual else NA_real_,
                         max_rhat = max(rhat), min_ess = min(ess), max_mcse = max(mcse),
                         stringsAsFactors = FALSE),
    moments = moments, fit_history = history,
    parameters = list(h = stats::setNames(h, ids), J = J))
  if (!converged) {
    cm_warn(sprintf(
      "The sampled joint distribution is unresolved (calibration %s; sample residual %.2e; %s). This does not show that the inputs are infeasible; try more samples, chains or iterations, or check_feasibility().",
      if (fit_converged) "converged" else "did not converge", resid_raw,
      if (calibrate) sprintf("raking residual %.2e", cal$residual) else "no raking"),
      class = "deconflate_nonconvergence")
  }
  structure(list(cells = cells, prob = prob, diseases = ids, converged = converged,
                 iterations = it, max_residual = resid_final,
                 constrained_pairs = cons[, c("disease1", "disease2", "measure", "status")],
                 targets = targets, backend = "sampled", diagnostics = diag),
            class = "cm_joint")
}

# Check that a joint distribution belongs to a population: same diseases in
# the same order, matching marginals, pairwise tables and three-way terms,
# and converged. Impact-only changes do not invalidate a joint.
validate_joint <- function(joint, model, tol = 1e-7) {
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().")
  ids <- model$diseases$id
  stale <- function(why) {
    cm_abort(sprintf("The supplied joint distribution does not match the model (%s); refit it with fit_joint().", why))
  }
  if (!identical(joint$diseases, ids)) stale("different diseases or order")
  if (!isTRUE(joint$converged)) {
    cm_abort("The supplied joint distribution did not converge; its probabilities cannot be used.",
             class = "deconflate_nonconvergence")
  }
  tp <- joint$targets$prob
  if (is.null(tp) || length(tp) != length(ids) || max(abs(tp - model$diseases$prob)) > tol) {
    stale("disease probabilities")
  }
  pt <- pair_tables(model)
  cons <- pt[!is.na(pt$p11), , drop = FALSE]
  tg <- joint$targets$pairs
  if (is.null(tg)) stale("no recorded targets")
  key_m <- pair_key(cons$disease1, cons$disease2)
  key_j <- pair_key(tg$disease1, tg$disease2)
  if (!setequal(key_m, key_j)) stale("different constrained pairs")
  if (max(abs(cons$p11 - tg$p11[match(key_m, key_j)]), 0) > tol) stale("pairwise associations")
  tw_m <- model$three_way
  tw_j <- joint$targets$three_way
  sig <- function(tw) {
    if (is.null(tw) || !nrow(tw)) return(character(0))
    sort(paste(apply(tw[, c("disease1", "disease2", "disease3")], 1,
                     function(x) paste(sort(x), collapse = "|")),
               signif(tw$ratio, 10)))
  }
  if (!identical(sig(tw_m), sig(tw_j))) stale("three-way terms")
  invisible(joint)
}

#' Probabilities of disease combinations
#'
#' @param joint A [fit_joint()] result.
#' @param min_prob Drop combinations with probability below this value.
#' @return A data frame with one row per combination: a 0/1 column per
#'   disease, the number of diseases present and the probability, sorted by
#'   decreasing probability.
#' @export
combination_probs <- function(joint, min_prob = 0) {
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().")
  out <- as.data.frame(joint$cells)
  out$n_diseases <- rowSums(joint$cells)
  out$prob <- joint$prob
  out <- out[out$prob >= min_prob, , drop = FALSE]
  out <- out[order(-out$prob), , drop = FALSE]
  rownames(out) <- NULL
  out
}
