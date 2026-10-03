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
#' * three-way terms: the ratio of conditional odds ratios (see
#'   [cm_three_way()]);
#' * impacts and interactions: their own units.
#'
#' A draw that produces an impossible input (a probability outside (0, 1), a
#' non-positive odds ratio, a conditional probability below 0 or above 1) is
#' rejected by [cm_monte_carlo()] and counted.
#'
#' Keys (as in `params` of [cm_monte_carlo()]): `prob:<disease>`,
#' `assoc:<d1>:<d2>`, `three:<d1>:<d2>:<d3>`, `impact:<disease>` and
#' `inter:<d1>:<d2>`.
#'
#' @param model A [cm_model()].
#' @param diseases Named list of `cm_dist`, names = disease ids.
#' @param associations Named list of `cm_dist`, names `"d1:d2"` (either
#'   order) matching rows of the model's associations.
#' @param impacts Named list of `cm_dist`, names = disease ids.
#' @param interactions Named list of `cm_dist`, names `"d1:d2"` (either order)
#'   matching rows of the model's interactions.
#' @param three_way Named list of `cm_dist`, names `"d1:d2:d3"` (any order)
#'   matching rows of the model's three-way terms. Three-way terms affect the
#'   global method only.
#'
#' @return A function of the draw index returning a [cm_model()], with
#'   attribute `specs` (the distributions, keyed as in `params`). It also
#'   accepts `u` (named uniforms, for stratified sampling) and `values` (named
#'   input values that replace draws, for importance sampling and for batch
#'   runs); [cm_monte_carlo()] uses these.
#' @export
#' @examples
#' s <- cm_sampler(example_supplement(),
#'                 associations = list("d1:d2" = dist_lognormal_ci(2, 1.4, 2.9)),
#'                 impacts = list(d1 = dist_normal(2.5, 0.5)))
#' s(1)$associations
cm_sampler <- function(model, diseases = list(), associations = list(),
                       impacts = list(), interactions = list(), three_way = list()) {
  check_model(model)
  check_dists <- function(x, what) {
    if (length(x) && (is.null(names(x)) || !all(vapply(x, inherits, logical(1), "cm_dist")))) {
      cm_abort(sprintf("`%s` must be a named list of cm_dist objects.", what))
    }
    if (anyDuplicated(names(x))) {
      cm_abort(sprintf("`%s` names an input more than once: %s.", what,
                       paste(unique(names(x)[duplicated(names(x))]), collapse = ", ")))
    }
  }
  check_dists(diseases, "diseases")
  check_dists(associations, "associations")
  check_dists(impacts, "impacts")
  check_dists(interactions, "interactions")
  check_dists(three_way, "three_way")
  specs <- list()
  rows <- list()

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
      r <- if (length(ab) == 2L) match(pair_key(ab[1], ab[2]), akey) else NA_integer_
      if (is.na(r)) {
        cm_abort(sprintf("'%s' does not match an association row (add it to cm_associations() first).", nm))
      }
      if (a$measure[r] %in% c("unknown", "independent", "table")) {
        cm_abort(sprintf("Association '%s' has measure '%s'; give it a numeric measure to vary it.", nm, a$measure[r]))
      }
      key <- paste0("assoc:", a$disease1[r], ":", a$disease2[r])
      if (!is.null(specs[[key]])) cm_abort(sprintf("Association '%s' is given twice.", nm))
      specs[[key]] <- associations[[nm]]
      rows[[key]] <- list(table = "associations", row = r)
    }
  }
  if (length(impacts)) {
    im <- model$impacts
    if (is.null(im)) cm_abort("The model has no impacts to vary.")
    for (nm in names(impacts)) {
      r <- match(nm, im$disease)
      if (is.na(r)) cm_abort(sprintf("'%s' does not match an impact (use the disease id).", nm))
      key <- paste0("impact:", nm)
      specs[[key]] <- impacts[[nm]]
      rows[[key]] <- list(table = "impacts", row = r)
    }
  }
  if (length(interactions)) {
    it <- model$interactions
    if (is.null(it)) cm_abort("The model has no interactions to vary.")
    ikey <- pair_key(it$disease1, it$disease2)
    for (nm in names(interactions)) {
      ab <- strsplit(nm, ":", fixed = TRUE)[[1]]
      r <- if (length(ab) == 2L) match(pair_key(ab[1], ab[2]), ikey) else NA_integer_
      if (is.na(r)) cm_abort(sprintf("'%s' does not match an interaction (use 'd1:d2').", nm))
      key <- paste0("inter:", it$disease1[r], ":", it$disease2[r])
      if (!is.null(specs[[key]])) cm_abort(sprintf("Interaction '%s' is given twice.", nm))
      specs[[key]] <- interactions[[nm]]
      rows[[key]] <- list(table = "interactions", row = r)
    }
  }

  if (length(three_way)) {
    tw <- model$three_way
    if (is.null(tw) || !nrow(tw)) cm_abort("The model has no three-way terms to vary.")
    tkey <- apply(tw[, c("disease1", "disease2", "disease3")], 1,
                  function(x) paste(sort(x), collapse = "|"))
    for (nm in names(three_way)) {
      abc <- strsplit(nm, ":", fixed = TRUE)[[1]]
      r <- if (length(abc) == 3L) match(paste(sort(abc), collapse = "|"), tkey) else NA_integer_
      if (is.na(r)) cm_abort(sprintf("'%s' does not match a three-way term (use 'd1:d2:d3').", nm))
      key <- paste0("three:", tw$disease1[r], ":", tw$disease2[r], ":", tw$disease3[r])
      if (!is.null(specs[[key]])) cm_abort(sprintf("Three-way term '%s' is given twice.", nm))
      specs[[key]] <- three_way[[nm]]
      rows[[key]] <- list(table = "three_way", row = r)
    }
  }

  draw_values <- function(u_in = NULL, fixed = NULL) {
    vapply(names(specs), function(k) {
      if (!is.null(fixed) && k %in% names(fixed)) return(as.numeric(fixed[[k]]))
      uk <- if (!is.null(u_in) && k %in% names(u_in)) u_in[[k]] else stats::runif(1)
      specs[[k]]$q(uk)
    }, numeric(1))
  }

  fn <- function(i, u = NULL, values = NULL) {
    v <- draw_values(u, values)
    m <- model
    for (k in names(v)) {
      info <- rows[[k]]
      x <- v[[k]]
      if (!is.finite(x)) cm_abort(sprintf("Draw %s is not finite.", k), class = "deconflate_infeasible")
      switch(info$table,
        diseases = {
          type <- m$diseases$type[info$row]
          prob <- if (type == "incidence_rate") 1 - exp(-x) else x
          if (prob <= 0 || prob >= 1 || x < 0) {
            cm_abort(sprintf("Draw %s = %g gives an invalid probability.", k, x),
                     class = "deconflate_infeasible")
          }
          m$diseases$value[info$row] <- x
          m$diseases$prob[info$row] <- prob
        },
        associations = {
          meas <- m$associations$measure[info$row]
          ok <- switch(meas, OR = , RR = x > 0, cond_prob = x >= 0 && x <= 1,
                       RD = , phi = x >= -1 && x <= 1, TRUE)
          if (!ok) {
            cm_abort(sprintf("Draw %s = %g is not a valid %s.", k, x, meas),
                     class = "deconflate_infeasible")
          }
          m$associations$value[info$row] <- x
        },
        impacts = {
          m$impacts$value[info$row] <- x
        },
        interactions = {
          m$interactions$value[info$row] <- x
        },
        three_way = {
          if (x <= 0) {
            cm_abort(sprintf("Draw %s = %g is not a valid ratio (it must be positive).", k, x),
                     class = "deconflate_infeasible")
          }
          m$three_way$ratio[info$row] <- x
        }
      )
    }
    m
  }
  attr(fn, "specs") <- specs
  class(fn) <- c("cm_sampler", "function")
  fn
}

#' Batch sampler: several analyses on shared population draws
#'
#' Builds one [cm_sampler()] per analysis of a [cm_analyses()] object. In
#' [cm_monte_carlo()], each draw of the disease probabilities and associations
#' is shared by all analyses, and each analysis draws its own impacts and
#' interactions. Draw identifiers are the same across analyses.
#'
#' @param analyses A [cm_analyses()] object.
#' @param diseases,associations,three_way Named lists of `cm_dist` for the
#'   shared population inputs (as in [cm_sampler()]).
#' @param impacts,interactions Named lists, one element per analysis, each a
#'   named list of `cm_dist` as in [cm_sampler()].
#' @return A `cm_batch_sampler` object.
#' @export
cm_batch_sampler <- function(analyses, diseases = list(), associations = list(),
                             impacts = list(), interactions = list(), three_way = list()) {
  if (!inherits(analyses, "cm_analyses")) cm_abort("`analyses` must come from cm_analyses().")
  nms <- names(analyses$models)
  unk <- setdiff(c(names(impacts), names(interactions)), nms)
  if (length(unk)) cm_abort(sprintf("Unknown analyses: %s.", paste(unk, collapse = ", ")))
  samplers <- lapply(nms, function(nm) {
    cm_sampler(analyses$models[[nm]], diseases = diseases, associations = associations,
               impacts = impacts[[nm]] %||% list(), interactions = interactions[[nm]] %||% list(),
               three_way = three_way)
  })
  names(samplers) <- nms
  pop_keys <- names(attr(samplers[[1]], "specs"))
  pop_keys <- pop_keys[grepl("^(prob|assoc|three):", pop_keys)]
  structure(list(samplers = samplers, population_keys = pop_keys,
                 population_specs = attr(samplers[[1]], "specs")[pop_keys]),
            class = "cm_batch_sampler")
}

#' @export
print.cm_sampler <- function(x, ...) {
  specs <- attr(x, "specs")
  cat(sprintf("<cm_sampler> %d uncertain inputs\n", length(specs)))
  for (k in utils::head(names(specs), 10)) cat("  ", k, ": ", specs[[k]]$type, "\n", sep = "")
  if (length(specs) > 10) cat(sprintf("  ... and %d more\n", length(specs) - 10))
  invisible(x)
}

#' @export
print.cm_batch_sampler <- function(x, ...) {
  cat(sprintf("<cm_batch_sampler> %d analyses (%s); %d shared population inputs\n",
              length(x$samplers), paste(names(x$samplers), collapse = ", "),
              length(x$population_keys)))
  invisible(x)
}
