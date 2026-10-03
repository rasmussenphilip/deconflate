# Shared test models.

# Supplementary File example with culling hazard ratios 1.5, 2.0 and 1.3.
supp_hr <- function() {
  m <- example_supplement()
  m$impacts <- combine_impacts(
    m$impacts,
    cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), outcome = "culling",
               scale = "hazard_ratio")
  )
  m
}
