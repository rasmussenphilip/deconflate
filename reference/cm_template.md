# Write a template set of input files

Writes example `diseases.csv`, `associations.csv`, `impacts.csv`,
`interactions.csv` (header only) and `uncertainty.csv` files to a
folder, to edit and read back with
[`cm_read_inputs()`](https://rasmussenphilip.github.io/deconflate/reference/cm_read_inputs.md).
The example has three dairy diseases, with yield impacts and culling
hazard ratios; the values are illustrative.

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
#> Wrote 5 files to /tmp/Rtmpbx0thx/my-inputs
```
