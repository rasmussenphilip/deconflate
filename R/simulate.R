#' Simulate the raw impacts that studies would report
#'
#' For validation: given a population whose disease combinations follow the
#' maximum-entropy distribution of `model` (see [fit_joint()]), and known true
#' additive impacts (and optional interactions), computes what a study would
#' estimate for each impact row: the crude difference between the two groups
#' the row compares, or the regression coefficient of an adjusted row (see
#' [cm_impacts()] for what a row means). With `n = NULL` the exact
#' population values are returned; with `n` set, `n` animals are sampled,
#' which adds sampling error.
#'
#' Adjusting the returned raw impacts with [deconflate()] recovers
#' `true_impacts` and `interactions` exactly (for `n = NULL`), whenever the
#' rows determine them.
#'
#' @param model A [cm_population()] or [cm_model()] (its impacts are ignored).
#' @param true_impacts Named numeric vector of true impacts, one per disease
#'   (any units).
#' @param interactions Optional named numeric vector of true interactions,
#'   named by combination (e.g. `c("d1;d2" = 0.5)`), in the same units.
#' @param rows Optional data frame of the rows to simulate, with column
#'   `disease` and optional `among` and `adjusted_for`. Default: one row per
#'   disease (with `adjusted_for`) and one crude row per interaction.
#' @param adjusted_for For the default rows: adjustment sets (ids separated
#'   by `";"`, or `"all"`), one per disease or recycled.
#' @param n Optional number of animals to sample.
#' @param joint Optional pre-computed [fit_joint()] result. It must match the
#'   model and have converged; a joint fitted here that does not converge
#'   stops with class `deconflate_nonconvergence`.
#' @param units Optional units label for the returned impacts.
#' @param ... Not used (arguments of earlier versions give a message).
#'
#' @return A [cm_impacts()] object with the simulated raw impacts.
#' @export
#' @examples
#' raw <- simulate_raw_impacts(example_supplement(), c(d1 = 2, d2 = 4, d3 = 6),
#'                             interactions = c("d1;d2" = 1))
#' raw
#' deconflate(cm_model(example_supplement(), raw))$adjusted
simulate_raw_impacts <- function(model, true_impacts, interactions = NULL, rows = NULL,
                                 adjusted_for = NA_character_, n = NULL, joint = NULL, units = NULL, ...) {
  old_args(list(...), c(estimand = estimand_message))
  check_population(model)
  ids <- model$diseases$disease
  if (is.null(names(true_impacts)) || !setequal(names(true_impacts), ids)) {
    cm_abort("`true_impacts` must be named, with one value per disease.", class = "deconflate_invalid_input")
  }
  if (!is.null(interactions)) {
    if (!is.numeric(interactions) || is.null(names(interactions))) {
      cm_abort("`interactions` must be a named numeric vector, named by combination (e.g. c(\"d1;d2\" = 0.5)).", class = "deconflate_invalid_input")
    }
  }
  combos <- lapply(names(interactions), function(t) parse_members(t, ids, "interactions", "name"))
  if (any(lengths(combos) < 2L)) cm_abort("Each interaction must name a combination of two or more diseases.", class = "deconflate_invalid_input")
  members <- c(as.list(ids), combos)
  names(members) <- c(ids, vapply(combos, term_label, character(1)))
  if (is.null(rows)) {
    rows <- data.frame(disease = c(ids, names(members)[-seq_along(ids)]),
                       adjusted_for = c(recycle_arg(adjusted_for, length(ids), "adjusted_for"),
                                        rep(NA_character_, length(combos))),
                       stringsAsFactors = FALSE)
  }
  if (!is.data.frame(rows) || is.null(rows$disease)) cm_abort("`rows` must be a data frame with a column `disease`.", class = "deconflate_invalid_input")
  rows$row <- seq_len(nrow(rows))
  if (is.null(rows$among)) rows$among <- NA_character_
  if (is.null(rows$adjusted_for)) rows$adjusted_for <- NA_character_
  rows$among <- clean_na(rows$among)
  rows$adjusted_for <- clean_na(rows$adjusted_for)
  rows$value <- 0
  specs <- lapply(seq_len(nrow(rows)), function(r) impact_spec(rows[r, , drop = FALSE], ids))
  # A study's "all" covers every term of the true model.
  extra <- setdiff(unlist(lapply(specs, function(s) if (length(s$g) > 1L) s$term)), names(members))
  if (length(extra)) {
    members <- c(members, stats::setNames(lapply(extra, split_ids), extra))
  }
  # A fitted joint that did not converge does not satisfy the inputs, so it
  # cannot produce known-truth data: stop in both cases.
  if (is.null(joint)) {
    joint <- withCallingHandlers(fit_joint(model),
                                 deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
  }
  validate_joint(joint, model)
  cells <- joint$cells
  y <- as.vector(cells %*% true_impacts[ids])
  for (j in seq_along(interactions)) {
    y <- y + interactions[[j]] * all_present(cells, combos[[j]])
  }
  if (is.null(n)) {
    X0 <- cells
    w <- joint$prob
    yy <- y
  } else {
    draw <- sample.int(nrow(cells), n, replace = TRUE, prob = joint$prob)
    X0 <- cells[draw, , drop = FALSE]
    w <- rep(1, n)
    yy <- y[draw]
  }
  raw <- vapply(specs, function(s) {
    A <- all_present(X0, s$apos) & none_present(X0, s$aneg) & w > 0
    XG <- all_present(X0, s$g)
    if (!s$adjusted) {
      w1 <- w * (A & XG)
      w0 <- w * (A & none_present(X0, s$g))
      return(sum(w1 * yy) / sum(w1) - sum(w0 * yy) / sum(w0))
    }
    idx <- which(A)
    cols <- list(as.numeric(XG[idx]))
    nm <- s$term
    add <- function(name, v) {
      if (!(name %in% nm)) {
        cols[[length(cols) + 1L]] <<- v
        nm <<- c(nm, name)
      }
    }
    if (length(s$g) > 1L) {
      for (d in s$g) add(d, X0[idx, d])
      for (sub in smaller_combinations(s$g)) add(term_label(sub), as.numeric(all_present(X0[idx, , drop = FALSE], sub)))
    }
    for (d in s$adj) add(d, X0[idx, d])
    if (s$adj_all) {
      for (t in names(members)) add(t, as.numeric(all_present(X0[idx, , drop = FALSE], members[[t]])))
    }
    Z <- cbind(1, do.call(cbind, cols))
    sw <- sqrt(w[idx])
    q <- qr(Z * sw, tol = 1e-9)
    keep <- sort(q$pivot[seq_len(q$rank)])
    unname(qr.coef(qr(Z[, keep, drop = FALSE] * sw), yy[idx] * sw)[which(keep == 2L)])
  }, numeric(1))
  cm_impacts(rows$disease, raw, among = rows$among, adjusted_for = rows$adjusted_for,
             source = if (is.null(n)) "simulated (expected)" else sprintf("simulated (n = %d)", n),
             units = units)
}
