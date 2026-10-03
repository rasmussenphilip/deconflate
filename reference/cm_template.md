# Write a template set of input files

Writes an example set of input files to a folder, to edit and read back
with
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md):
three dairy diseases with their associations, two analyses
(`impacts_yield.csv` in percent of yield and
`impacts_calving_interval.csv` in days), culling hazard ratios, an empty
interactions file for the yield analysis, an empty three-way file and an
uncertainty file. The values are illustrative.

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

## Examples

``` r
cm_template(file.path(tempdir(), "my-inputs"), overwrite = TRUE)
#> Wrote 8 files to /tmp/Rtmp8dzg05/my-inputs
```
