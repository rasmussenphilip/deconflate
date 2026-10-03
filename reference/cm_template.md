# Write a template set of input files

Writes an example set of input files to a folder, to edit and read back
with
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md).
The values are illustrative. Some values have a distribution in their
row (columns `dist` and `p1`-`p4`); the others are point values.

## Usage

``` r
cm_template(dir, type = c("single", "analyses"), overwrite = FALSE)
```

## Arguments

- dir:

  Folder to write to (created if needed).

- type:

  `"single"` or `"analyses"`.

- overwrite:

  Overwrite existing files?

## Value

The file paths, invisibly.

## Details

- `type = "single"` (default): one analysis, the standard workflow.
  Three dairy diseases (`diseases.csv`) and their associations
  (`associations.csv`), one impact vector (`impacts.csv`, in percent of
  yield; any units can be used), an empty `interactions.csv` (pairwise
  interactions in the same units) and an empty `three_way.csv`. Read
  back, it gives a
  [`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md)
  and a
  [`cm_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_sampler.md),
  which supports Latin hypercube and importance sampling.

- `type = "analyses"`: two analyses on one population
  (`impacts_yield.csv` in percent of yield and
  `impacts_calving_interval.csv` in days), with
  `interactions_yield.csv`. Read back, it gives a
  [`cm_analyses()`](https://rasmussenphilip.github.io/deconflate/reference/cm_analyses.md)
  object and a
  [`cm_batch_sampler()`](https://rasmussenphilip.github.io/deconflate/reference/cm_batch_sampler.md)
  (shared population draws; simple random sampling only).

Hazard ratios are a separate model
([`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md))
and are not part of either template; see
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md)
for the `hazard_ratios.csv` columns. The folder
`system.file("extdata", "five_diseases", package = "deconflate")` has a
fuller example that uses every feature.

## Examples

``` r
dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 5 files to /tmp/RtmpNGVr3j/my-inputs
inp <- cm_read_inputs(dir = dir)
deconflate(inp$model)
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted   change
#>      LAM 4.81    2.871 -0.40320
#>      SCK 8.40    7.819 -0.06915
#>      MET 5.61    3.170 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0
# \donttest{
# Importance sampling of the lameness yield impact, with a defensive
# mixture that covers its whole support
specs <- attr(inp$sampler, "specs")
prop <- list("impact:LAM" = dist_mixture(specs[["impact:LAM"]], dist_normal(6, 1.5),
                                         weights = c(0.5, 0.5)))
mc <- cm_monte_carlo(inp$sampler, 200, proposal = prop, seed = 1)
summary(mc, diagnose = FALSE)
#>   disease       method     mean        sd       mcse    q0.025     q0.5
#> 1     LAM simultaneous 2.780003 1.0144557 0.07196501 0.7373599 2.796856
#> 2     SCK simultaneous 7.966686 1.2181885 0.09387597 5.7170650 7.940101
#> 3     MET simultaneous 3.180425 0.3882904 0.02769793 2.4727201 3.184312
#>      q0.975 trimmed_mean    rel_mcse tail_share stability
#> 1  4.709681     2.772107 0.025886662 0.05108681        ok
#> 2 10.075637     7.954862 0.011783566 0.09490480        ok
#> 3  3.881797     3.182747 0.008708877 0.05977454        ok
# }
```
