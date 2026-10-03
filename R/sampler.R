#' Build a Monte Carlo sampler from input distributions
#'
#' Creates a function that, for each draw, copies `model` and replaces
#' selected inputs with random draws from [distributions]. Inputs without a
#' distribution keep their values. The result is passed to
#' [cm_monte_carlo()].
#'
#' Draws are on the scale the inputs were entered on:
#' * diseases: the `value` given to [cm_diseases()] (e.g. an incidence rate,
#'   converted to a probability as specified by its `type`);
#' * associations: the association measure (e.g. odds ratio);
#' * impacts and interactions: the input scale (e.g. percent, if entered as
#'   percent).
#'
#' A draw that produces an impossible input (a probability outside (0, 1), a
#' non-positive odds ratio) is rejected by [cm_monte_carlo()] and counted.
#'
#' @param model A [cm_model()].
#' @param diseases Named list of `cm_dist`, names = disease ids.
#' @param associations Named list of `cm_dist`, names `"d1:d2"` (either
#'   order) matching rows of the model's associations.
#' @param impacts Named list of `cm_dist`, names `"outcome:disease"`.
#' @param interactions Named list of `cm_dist`, names `"outcome:d1:d2"`.
#' @param outcome_correlation Optional correlation matrix with row and column
#'   names equal to outcome labels. Within each disease, the impacts on these
#'   outcomes are drawn with this correlation (Gaussian copula); marginal
#'   distributions are unchanged. Impacts of different diseases stay
#'   independent.
#'
#' @return A function of the draw index returning a [cm_model()], with
#'   attributes `specs` (the distributions, keyed as in `params` of
#'   [cm_monte_carlo()]) and `correlated` (keys drawn through the copula).
#'   The function also accepts `u` (named uniforms, for stratified sampling)
#'   and `values` (named input values that replace draws, for importance
#'   sampling); [cm_monte_carlo()] uses these.
#' @export
#' @examples
#' m <- example_supplement()
#' s <- cm_sampler(m, associations = list("d1:d2" = dist_lognormal_ci(2, 1.4, 2.9)))
#' s(1)$associations
cm_sampler <- function(model, diseases = list(), associations = list(),
                       impacts = list(), interactions = list(),
                       outcome_correlation = NULL) {
  check_model(model)
  specs <- list()
  rows <- list()

  check_dists <- function(x, what) {
    if (length(x) && (is.null(names(x)) || !all(vapply(x, inherits, logical(1), "cm_dist")))) {
      cm_abort(sprintf("`%s` must be a named list of cm_dist objects.", what))
    }
  }
  check_dists(diseases, "diseases")
  check_dists(associations, "associations")
  check_dists(impacts, "impacts")
  check_dists(interactions, "interactions")

  for (id in names(diseases)) {
    r <- match(id, model$diseases$id)
    if (is.na(r)) cm_abort(sprintf("Unknown disease '%s' in `diseases`.", id))
    key <- paste0("prob:", id)
    specs[[key]] <- diseases[[id]]
    rows[[key]] <- list(table = "diseases", row = r)
  }
  if (length(associations)) {
    a <- model$associations
    if (is.null(a)) cm_abort("The model has no associations to vary.")
    akey <- pair_key(a$disease1, a$disease2)
    for (nm in names(associations)) {
      ab <- strsplit(nm, ":", fixed = TRUE)[[1]]
      r <- match(pair_key(ab[1], ab[2]), akey)
      if (length(ab) != 2L || is.na(r)) {
        cm_abort(sprintf("'%s' does not match an association row (add it to cm_associations() first).", nm))
      }
      if (a$measure[r] %in% c("unknown", "independent", "table")) {
        cm_abort(sprintf("Association '%s' has measure '%s'; give it a numeric measure to vary it.", nm, a$measure[r]))
      }
      key <- paste0("assoc:", a$disease1[r], ":", a$disease2[r])
      specs[[key]] <- associations[[nm]]
      rows[[key]] <- list(table = "associations", row = r)
    }
  }
  if (length(impacts)) {
    im <- model$impacts
    if (is.null(im)) cm_abort("The model has no impacts to vary.")
    ikey <- paste(im$outcome, im$disease, sep = ":")
    for (nm in names(impacts)) {
      r <- match(nm, ikey)
      if (is.na(r)) cm_abort(sprintf("'%s' does not match an impact (use 'outcome:disease').", nm))
      key <- paste0("impact:", nm)
      specs[[key]] <- impacts[[nm]]
      rows[[key]] <- list(table = "impacts", row = r,
                          divisor = if (identical(im$input_scale[r], "percent")) 100 else 1)
    }
  }
  if (length(interactions)) {
    it <- model$interactions
    if (is.null(it)) cm_abort("The model has no interactions to vary.")
    for (nm in names(interactions)) {
      parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
      r <- if (length(parts) == 3L) {
        which(it$outcome == parts[1] & pair_key(it$disease1, it$disease2) == pair_key(parts[2], parts[3]))
      } else integer(0)
      if (length(r) != 1L) cm_abort(sprintf("'%s' does not match an interaction (use 'outcome:d1:d2').", nm))
      key <- paste0("inter:", it$outcome[r], ":", it$disease1[r], ":", it$disease2[r])
      specs[[key]] <- interactions[[nm]]
      rows[[key]] <- list(table = "interactions", row = r,
                          divisor = if (identical(it$input_scale[r], "percent")) 100 else 1)
    }
  }

  # Copula groups: for each disease, the impact keys on the correlated outcomes.
  groups <- list()
  chol_R <- NULL
  if (!is.null(outcome_correlation)) {
    R <- as.matrix(outcome_correlation)
    outs <- rownames(R)
    if (is.null(outs) || !identical(outs, colnames(R))) {
      cm_abort("`outcome_correlation` needs identical row and column names (outcomes).")
    }
    if (!isSymmetric(R) || any(abs(diag(R) - 1) > 1e-12)) {
      cm_abort("`outcome_correlation` must be a symmetric correlation matrix.")
    }
    chol_R <- tryCatch(chol(R), error = function(e) cm_abort("`outcome_correlation` is not positive definite."))
    for (id in model$diseases$id) {
      keys <- paste0("impact:", outs, ":", id)
      if (any(keys %in% names(specs))) groups[[id]] <- keys
    }
  }
  correlated <- unlist(groups, use.names = FALSE)
  correlated <- correlated[correlated %in% names(specs)]

  draw_values <- function(u_in = NULL, fixed = NULL) {
    u <- list()
    for (g in groups) {
      z <- as.vector(stats::rnorm(length(g)) %*% chol_R)
      uu <- stats::pnorm(z)
      for (j in seq_along(g)) if (g[j] %in% names(specs)) u[[g[j]]] <- uu[j]
    }
    vapply(names(specs), function(k) {
      if (!is.null(fixed) && k %in% names(fixed)) return(as.numeric(fixed[[k]]))
      uk <- if (!is.null(u[[k]])) {
        u[[k]]
      } else if (!is.null(u_in) && k %in% names(u_in)) {
        u_in[[k]]
      } else {
        stats::runif(1)
      }
      specs[[k]]$q(uk)
    }, numeric(1))
  }

  # `u`: optional named uniforms (e.g. Latin hypercube) for uncorrelated
  # inputs; `values`: optional named input values that replace draws (used
  # for importance sampling by cm_monte_carlo()).
  fn <- function(i, u = NULL, values = NULL) {
    v <- draw_values(u, values)
    m <- model
    for (k in names(v)) {
      info <- rows[[k]]
      x <- v[[k]]
      switch(info$table,
        diseases = {
          type <- m$diseases$type[info$row]
          prob <- if (type == "incidence_rate") 1 - exp(-x) else x
          if (!is.finite(prob) || prob <= 0 || prob >= 1 || x < 0) {
            cm_abort(sprintf("Draw %s = %g gives an invalid probability.", k, x),
                     class = "deconflate_infeasible")
          }
          m$diseases$value[info$row] <- x
          m$diseases$prob[info$row] <- prob
        },
        associations = {
          meas <- m$associations$measure[info$row]
          if (meas %in% c("OR", "RR") && !(x > 0)) {
            cm_abort(sprintf("Draw %s = %g is not a valid %s.", k, x, meas),
                     class = "deconflate_infeasible")
          }
          m$associations$value[info$row] <- x
        },
        impacts = {
          if (identical(m$impacts$scale[info$row], "hazard_ratio") && !(x > 0)) {
            cm_abort(sprintf("Draw %s = %g is not a valid hazard ratio.", k, x),
                     class = "deconflate_infeasible")
          }
          m$impacts$value[info$row] <- x / info$divisor
        },
        interactions = {
          m$interactions$value[info$row] <- x / info$divisor
        }
      )
    }
    m
  }
  attr(fn, "specs") <- specs
  attr(fn, "correlated") <- correlated
  class(fn) <- c("cm_sampler", "function")
  fn
}

#' @export
print.cm_sampler <- function(x, ...) {
  specs <- attr(x, "specs")
  cat(sprintf("<cm_sampler> %d uncertain inputs", length(specs)))
  corr <- attr(x, "correlated")
  if (length(corr)) cat(sprintf(" (%d drawn with outcome correlation)", length(corr)))
  cat("\n")
  for (k in utils::head(names(specs), 10)) cat("  ", k, ": ", specs[[k]]$type, "\n", sep = "")
  if (length(specs) > 10) cat(sprintf("  ... and %d more\n", length(specs) - 10))
  invisible(x)
}
