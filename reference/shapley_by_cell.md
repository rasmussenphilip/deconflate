# Shapley attribution of a general loss function over disease combinations

Computes each disease's expected Shapley value for an arbitrary loss
function of the diseases an animal has. The calculation runs cell by
cell over the joint distribution
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)):
within a combination, only the diseases present are enumerated, so the
cost depends on how many diseases co-occur rather than on the total
number of diseases. Shapley values are linear in the game, so the
expected Shapley value equals the probability-weighted sum over cells.

## Usage

``` r
shapley_by_cell(joint, loss, max_present = Inf)
```

## Arguments

- joint:

  A
  [`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
  result.

- loss:

  A function of a named 0/1 vector (one element per disease) returning
  the loss for an animal with those diseases; a loss of zero for no
  disease is assumed. See
  [`loss_additive()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md)
  and
  [`loss_multiplicative()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md).

- max_present:

  Cells with more diseases present than this are skipped (the cost grows
  as 2^k for a cell with k diseases). The default skips nothing. Skipped
  cells are reported: their probability (`skipped_mass`) and their
  expected loss (`skipped_loss`), which is the part of the total that is
  not allocated.

## Value

A data frame with each disease's expected Shapley value and share (of
the allocated amount, which is the total only when no cell is skipped),
with attributes `total` (expected loss over all cells), `allocated` (the
sum of the Shapley values), `skipped_mass` and `skipped_loss`.

## Details

For additive loss functions the result has a closed form: each term is
split equally among the diseases it involves
([`loss_additive()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md)).
The cell-wise calculation is needed for non-additive losses such as
[`loss_multiplicative()`](https://rasmussenphilip.github.io/deconflate/reference/loss_functions.md).

## Examples

``` r
j <- fit_joint(example_supplement())
m <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)
shapley_by_cell(j, loss_multiplicative(m))
#>   disease     shapley      share
#> 1      d1 0.001978346 0.09981218
#> 2      d2 0.005922273 0.29879258
#> 3      d3 0.011920065 0.60139525
```
