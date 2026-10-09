# Write a template set of input files

Writes an example set of input files to a folder, to edit and read back
with
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md).
The values are illustrative. Some values have a distribution in their
row (columns `dist` and `p1`-`p4`); the others are point values. The
files are:

## Usage

``` r
cm_template(dir, overwrite = FALSE)
```

## Arguments

- dir:

  Folder to write to (created if needed).

- overwrite:

  Overwrite existing files?

## Value

The file paths, invisibly.

## Details

- `diseases.csv`: three dairy diseases;

- `associations.csv`: their associations (odds ratios);

- `yield.csv`: additive impacts on milk yield, in percent of yield (any
  units can be used);

- `yield_interactions.csv`: an empty interactions table for the yield
  impacts (pairwise interactions in the same units);

- `culling.csv`: event impacts on culling (hazard ratios and a risk
  ratio), for `deconflate(..., event_model = TRUE, overall_risk = ...)`;

- `three_way.csv`: an empty three-way table.

The folder
`system.file("extdata", "five_diseases", package = "deconflate")` has a
fuller example that uses every feature.

## Examples

``` r
dir <- file.path(tempdir(), "my-inputs")
cm_template(dir, overwrite = TRUE)
#> Wrote 6 files to /tmp/RtmpNxOcWp/my-inputs
yield <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                        associations = file.path(dir, "associations.csv"),
                        impacts = file.path(dir, "yield.csv"))
deconflate(yield, n_draws = 200, seed = 1)
#> <cm_result> method: simultaneous; milk yield loss [% of yield]
#> 
#>  disease  raw adjusted lower  upper   change
#>      LAM 4.81    2.871 0.600  4.697 -0.40320
#>      SCK 8.40    7.819 5.475 10.065 -0.06915
#>      MET 5.61    3.170 2.302  4.010 -0.43496
#> 
#> Raw sum: 4.966; adjusted total: 4.015 (95% interval 3.111 to 4.932)
#> Diagnostics: residual 1.78e-15, condition number 2, sign changes 0
#> Uncertainty: 95% intervals from 200 draws (0 rejected; random sampling; seed 1).
culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = file.path(dir, "associations.csv"),
                          impacts = file.path(dir, "culling.csv"))
deconflate(culling, event_model = TRUE, overall_risk = 0.25, n_draws = 0)
#> <cm_event_result> culling; method: snapshot
#> 
#>  disease measure  raw adjusted_hr
#>      LAM      HR 1.74       1.514
#>      SCK      HR 1.92       1.782
#>      MET      RR 1.45       1.201
#> 
#> Overall risk 0.25; disease-free risk 0.1767; attributable to disease 0.07328 (29.3% of the overall risk)
#> 
#> Attributable risk by disease (Shapley allocation):
#>  disease attributable   share
#>      LAM     0.022295 0.30424
#>      SCK     0.047051 0.64205
#>      MET     0.003936 0.05371
#> 
#> Diagnostics: residual 2.22e-16, condition number 2.2, sign changes 0
```
