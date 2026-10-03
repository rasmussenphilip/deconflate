# Worked example from the Supplementary File of Rasmussen et al. (2022)

Three hypothetical diseases with prevalences 0.10, 0.15 and 0.20, yield
impacts of 2.5%, 5% and 7.5%, and odds ratios of 2 (d1:d2), 1 (d1:d3)
and 3 (d2:d3). The observed mean yield in the example is 10,000 units
per animal per year.

## Usage

``` r
example_supplement()
```

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Details

Note: the published productivity gap (10,225; 20, 61 and 143 units)
rounds the adjusted impacts to 0.02, 0.04 and 0.07 before computing the
gap. Without rounding, the published method gives 10,215.7 and 21.1,
56.7 and 137.9 units.
