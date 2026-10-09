# Batch sampler for the reproduction of Rasmussen et al. (2024) ----------------
#
# Internal. Several analyses (named list of cm_model objects) share the
# population inputs: each draw of the disease and association inputs is used
# by every analysis, and each analysis draws its own impacts. Distributions
# are named lists of cm_dist objects, keyed as in cm_model(distributions = ):
# `pop` for the shared inputs ("prob:", "assoc:", "three:") and `impacts`
# (one list per analysis, keyed by disease id).
batch_sampler <- function(models, pop = list(), impacts = list()) {
  nms <- names(models)
  unk <- setdiff(names(impacts), nms)
  if (length(unk)) cm_abort(sprintf("Unknown analyses: %s.", paste(unk, collapse = ", ")))
  pop <- check_distributions(pop, models[[1]]) %||% list()
  samplers <- lapply(nms, function(nm) {
    m <- models[[nm]]
    imp <- impacts[[nm]] %||% list()
    ispecs <- if (length(imp)) stats::setNames(imp, paste0("impact:", names(imp))) else list()
    ispecs <- check_distributions(ispecs, m) %||% list()
    fn <- function(i, values) {
      x <- vapply(names(ispecs), function(k) ispecs[[k]]$r(1), numeric(1))
      set_inputs(m, c(values, x))
    }
    attr(fn, "specs") <- c(pop, ispecs)
    fn
  })
  names(samplers) <- nms
  structure(list(samplers = samplers, population_keys = names(pop), population_specs = pop),
            class = "cm_batch_sampler")
}
