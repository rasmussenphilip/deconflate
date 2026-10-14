# Disease combinations, subgroups and interactions

``` r

library(deconflate)
```

## The distribution of disease combinations

[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
fits the distribution of all 2^n combinations of n diseases, by
iterative proportional fitting from independence. The result is the
maximum-entropy distribution consistent with the disease probabilities
and every association row.

``` r

m <- example_supplement()
j <- fit_joint(m)
j
#> <cm_joint> 3 diseases, 8 combinations (exact backend)
#>   Converged: TRUE after 28 sweeps (max residual 8.28e-12)
#>   Association rows: 3
head(combination_probs(j))
#>   d1 d2 d3 n_diseases       prob
#> 1  0  0  0          0 0.64262256
#> 2  0  0  1          1 0.13185683
#> 3  0  1  0          1 0.07737744
#> 4  1  0  0          1 0.06410444
#> 5  0  1  1          2 0.04814317
#> 6  1  1  0          2 0.01589556
```

Enumerating 2^n combinations limits this exact fit to about 20 diseases.
For more, `fit_joint(backend = "sampled")` fits the same model by Monte
Carlo (Gibbs sampling with calibrated parameters, followed by raking to
the pairwise tables);
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
switches to it by itself above 20 diseases. The sampled backend handles
crude pair rows only. See
[`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md)
for how it works, its diagnostics and a 24-disease example.

## What maximum entropy assumes

With three or more diseases, the probabilities and pairwise associations
do not determine a unique distribution: there are 2^n - 1 free
probabilities and n(n + 1) / 2 constraints. Maximum entropy picks the
one that assumes nothing beyond the rows. It is an assumption, not
something the pairwise evidence identifies. With only rows about pairs:

- the fitted distribution is a log-linear model with all two-way terms
  and no three-way or higher terms. The odds ratio of a pair is the same
  whether or not a third disease is present (given the other diseases);
- three diseases still occur together more (or less) often than under
  independence, through the pairwise terms;
- a pair whose association is unknown has no term of its own. Its
  association is implied by the associations both diseases have with
  other diseases.

Rows about combinations or subgroups (below) add terms of their own, and
the fit matches them too.

Additive results without interactions need only the pairwise tables when
every pair has a crude row: then
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
uses those tables directly (`$diagnostics$method` is
`"pairwise shortcut"`) and gives the same answer as the full
distribution. The distribution matters for unknown pairs, for rows about
combinations or subgroups, for interactions, and for event impacts
(`deconflate(..., event_model = TRUE)`, which uses the whole
distribution; see
[`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md)).

If the rows cannot hold together, no distribution exists.
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
then stops, and
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)
identifies the conflicting pairs:

``` r

bad <- cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                     cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
check_feasibility(bad, method = "triples")
#> <cm_feasibility> NOT jointly feasible [triples]
#> 
#> Conflicting triples:
#>  disease1 disease2 disease3   gap
#>         a        b        c 0.226
```

## Unknown associations are not independence

A pair without an association row is unknown. Giving it an odds ratio of
1 states that the two diseases are unrelated, which is an assumption of
its own. Leaving it unknown lets the maximum-entropy fit imply an
association through the other diseases. Here `a` and `c` are unrelated
in one model and unknown in the other:

``` r

d <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
imp <- cm_impacts(c("a", "b", "c"), c(3, 4, 5), units = "%")
ind <- cm_model(d, imp, associations = cm_associations(c("a", "b", "a"), c("b", "c", "c"), c(3, 2, 1)))
unk <- cm_model(d, imp, associations = cm_associations(c("a", "b"), c("b", "c"), c(3, 2)))
deconflate(ind, n_draws = 0)$adjusted$adjusted
#> [1] 2.263264 2.946945 4.548702
res_unk <- deconflate(unk, n_draws = 0)
res_unk
#> <cm_result> full model [%]
#> 
#>  disease raw adjusted  change
#>        a   3    2.100 -0.3000
#>        b   4    2.987 -0.2532
#>        c   5    4.481 -0.1037
#> 
#> Raw sum: 3.05; adjusted total: 2.436
#> Diagnostics: residual 0.00e+00, condition number 1.7, sign changes 0
#> Unknown pairs: 1 without an association; the fit gave it an odds ratio of 1.19 (see $unknown_pairs).
```

`$unknown_pairs` lists the pairs that no row mentions and the odds
ratios the fit gave them. Both `a` and `c` are positively associated
with `b`, so the fit implies a positive association between them too:

``` r

res_unk$unknown_pairs
#>   disease1 disease2 fitted_or
#> 1        a        c  1.194529
```

The published approximation needs every pair, so
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
reports why it cannot run:

``` r

compare_methods(unk)$failed
#>                                                                           published 
#> "The published approximation needs a crude pair row for every pair (missing: a-c)."
```

## Association rows about combinations and subgroups

An association row relates `disease` (G) to `associated_with` (H).
Either can be a combination of diseases separated by `";"`; a
combination as the comparison (`associated_with`) compares animals with
all of its diseases to animals with none of them, and a combination as
the outcome (`disease`) is “all of them”. A row can hold within a
subgroup (`among`: `"d3"`, `"no d3"` or `"no others"`), and an odds
ratio between two diseases can be adjusted for other diseases
(`adjusted_for`). The measures:

| measure | meaning |
|----|----|
| `OR`, `RR`, `RD`, `phi` | the 2x2 table of G by H, within `among` |
| `prevalence` | P(G \| H, `among`); H may be empty: P(G \| `among`), or the prevalence of a combination |
| `OR_ratio` | the odds ratio of G and H among `among` (one disease, S) divided by that among animals without S |

The five-disease example has a second LAM-MAS odds ratio (pooled with
the first; below), and an `OR_ratio` row: the LAM-MAS odds ratio is 1.5
times higher among cows with subclinical ketosis (SCK) than among cows
without:

``` r

ex <- system.file("extdata", "five_diseases", package = "deconflate")
a <- read.csv(file.path(ex, "associations.csv"))
a[, c("disease", "associated_with", "among", "measure", "value", "study")]
#>    disease associated_with among    measure    value   study
#> 1      LAM             MAS               OR 1.800000 Study A
#> 2      LAM             MAS               OR 2.200000 Study B
#> 3      LAM             SCK               OR 2.000000 Study C
#> 4      MET              RP               OR 6.000000 Study D
#> 5      MET             SCK               RR 1.500000 Study E
#> 6      MAS             SCK       prevalence 0.350000 Study F
#> 7      MAS             MET               RD 0.080000 Study G
#> 8      LAM             MET              phi 0.050000 Study H
#> 9       RP             SCK               OR 1.429688 Study I
#> 10     LAM              RP               OR 1.000000        
#> 11     MAS              RP               OR 1.500000 Study J
#> 12     LAM             MAS   SCK   OR_ratio 1.500000 Study K
```

Pairwise rows cannot identify how often three diseases occur together;
an `OR_ratio` row can. A ratio of 1 is the maximum-entropy assumption.
[`set_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/set_three_way.md)
adds or replaces such a row (the ratio is the same whichever disease is
the subgroup), and every other row is still matched:

``` r

m3 <- set_three_way(m, "d1", "d2", "d3", ratio = 2)
m3$associations
#>   row disease associated_with among adjusted_for  measure value study   source
#> 1   1      d1              d2  <NA>         <NA>       OR     2  <NA>     <NA>
#> 2   2      d1              d3  <NA>         <NA>       OR     1  <NA>     <NA>
#> 3   3      d2              d3  <NA>         <NA>       OR     3  <NA>     <NA>
#> 4   4      d1              d2    d3         <NA> OR_ratio     2  <NA> scenario
#>   note
#> 1 <NA>
#> 2 <NA>
#> 3 <NA>
#> 4 <NA>
j3 <- fit_joint(m3)
c(maximum_entropy = sum(j$prob[rowSums(j$cells) == 3]),
  ratio_2 = sum(j3$prob[rowSums(j3$cells) == 3]))
#> maximum_entropy         ratio_2 
#>     0.008583836     0.010703834
```

When every pair of the triple has a crude row, the `OR_ratio` row keeps
the pairwise tables fixed, so it changes additive results only through
interactions or rows about combinations or subgroups; without them the
result is unchanged:

``` r

res3 <- deconflate(m3, n_draws = 0)
res3$diagnostics$method
#> [1] "full model"
all.equal(res3$adjusted$adjusted, deconflate(m, n_draws = 0)$adjusted$adjusted)
#> [1] TRUE
```

It does change results that use the whole distribution, such as event
impacts (here hazard ratios of culling, with an overall culling risk of
0.25; see
[`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md)):

``` r

hr <- cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), measure = "HR", label = "culling")
e0 <- deconflate(cm_model(m, hr), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
e3 <- deconflate(cm_model(m3, hr), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
rbind(maximum_entropy = setNames(e0$adjusted$adjusted, e0$adjusted$disease),
      ratio_2 = setNames(e3$adjusted$adjusted, e3$adjusted$disease))
#>                       d1       d2       d3
#> maximum_entropy 1.383153 1.890932 1.143987
#> ratio_2         1.379893 1.891328 1.139917
```

An odds ratio adjusted for other diseases (from a logistic regression)
is not the crude odds ratio of the pair.
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
matches it as what it is: the coefficient of H in the logistic
regression of G on H and the adjustment diseases, within the fitted
distribution. Ids in `adjusted_for` that are not diseases of the model
(parity, herd, …) cannot be matched; the row is then used as crude, with
a note.

## Several rows for one item

Several studies can estimate the same thing. Association rows that fix
the same table (the same G, H, `among` and adjustment) are pooled on the
log scale, by random effects (DerSimonian-Laird) by default, or by fixed
effect (`deconflate(..., pool = "fixed")`). A row without a distribution
is exact. `$heterogeneity` gives Cochran’s Q, I^2 and tau^2, and
`$estimates` lists each row with its fitted value:

``` r

five <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
                       associations = file.path(ex, "associations.csv"),
                       impacts = file.path(ex, "impacts_calving_interval.csv"))
res_ci <- deconflate(five, n_draws = 0)
res_ci$heterogeneity
#>          table            item rows k          scale   pooled
#> 1 associations       LAM ~ MAS 1, 2 2 log odds ratio 1.960993
#> 2      impacts all impact rows  all 6          units       NA
#>                   method         Q df         p I2 tau2
#> 1                 random 0.6177079  1 0.4319003  0    0
#> 2 weighted least squares 0.3123068  1 0.5762681  0   NA
```

With random effects, the draws add the between-study spread (tau^2) to
the pooled value. The calving-interval table has two MAS rows (one
adjusted for SCK, one crude): impact rows are fitted together by
weighted least squares, with weights 1 / SD^2 from their distributions
(rows without a distribution are exact):

``` r

e <- res_ci$estimates
e[e$table == "impacts", c("row", "item", "adjusted_for", "among", "value", "sd", "fitted", "std_residual")]
#>   row item adjusted_for     among value       sd    fitted std_residual
#> 1   1  LAM         <NA>      <NA>    12 3.000000 12.027295 -0.009098366
#> 2   2  MAS          SCK      <NA>     6 2.000000  6.618529 -0.309264437
#> 3   3  MAS         <NA>      <NA>     9 3.000000  7.603910  0.465363382
#> 4   4  MET         <NA>      <NA>    18 4.463809 18.017756 -0.003977672
#> 5   5  SCK          all      <NA>     4 0.000000  4.000000           NA
#> 6   6   RP         <NA> no others    10 2.526123 10.002042 -0.000808178
res_ci$adjusted
#>   disease    kind raw  adjusted        change rows
#> 1     LAM disease  12 10.206306 -1.494745e-01    1
#> 2     MAS disease  NA  4.768205            NA    2
#> 3     MET disease  18 14.774125 -1.792153e-01    1
#> 4     SCK disease   4  4.000000  2.220446e-16    1
#> 5      RP disease  10 10.002042  2.041557e-04    1
```

`raw` is `NA` for a disease with several rows. The RP row compares cows
with RP and no other disease against healthy cows
(`among = "no others"`), which is RP’s own impact.

## Impact interactions

An interaction is the extra impact when diseases occur together, in the
units of the impact table: positive values are synergistic and negative
values antagonistic. Interactions cannot be inferred from associations;
they need evidence or explicit scenarios. An impact row about a
combination (`"d1;d2"`) adds an interaction term:

- crude (the default): animals with all of the diseases compared with
  animals with none of them;
- adjusted for `"all"`: the interaction itself (the coefficient of the
  product term in a regression with every disease and term of the
  model).

[`set_interaction()`](https://rasmussenphilip.github.io/deconflate/reference/set_interaction.md)
adds or replaces a row adjusted for `"all"`:

``` r

mi <- set_interaction(m, "d1", "d2", 1)   # 1 percentage point
mi$impacts
#>   row disease among adjusted_for measure value study
#> 1   1      d1  <NA>         <NA>    <NA>   2.5  <NA>
#> 2   2      d2  <NA>         <NA>    <NA>   5.0  <NA>
#> 3   3      d3  <NA>         <NA>    <NA>   7.5  <NA>
#> 4   4   d1;d2  <NA>          all    <NA>   1.0  <NA>
#>                            source note
#> 1 Rasmussen et al. 2022, Table S1 <NA>
#> 2 Rasmussen et al. 2022, Table S1 <NA>
#> 3 Rasmussen et al. 2022, Table S1 <NA>
#> 4                        scenario <NA>
res <- deconflate(mi, n_draws = 0)
res
#> <cm_result> full model; yield [%]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    1.914 -0.23445
#>       d2 5.0    3.241 -0.35187
#>       d3 7.5    6.936 -0.07525
#>    d1;d2 1.0    1.000       NA
#> 
#> Raw sum: 2.5; adjusted total: 2.089; interactions: 0.02448
#> Diagnostics: residual 8.88e-16, condition number 1.6, sign changes 0
```

The adjusted values reproduce every impact row: the raw values are
`A %*% adjusted`, where each row of `A` says how much of each term (main
impact or interaction) the row’s comparison contains under the fitted
distribution:

``` r

res$conflation$A
#>                     d1           d2           d3      d1;d2
#> impacts:1 1.000000e+00 1.053266e-01 2.058839e-11 0.24479393
#> impacts:2 7.434818e-02 1.000000e+00 2.096236e-01 0.16319595
#> impacts:3 1.158099e-11 1.670438e-01 1.000000e+00 0.02304974
#> impacts:4 0.000000e+00 2.038632e-16 0.000000e+00 1.00000000
as.vector(res$conflation$A %*% res$adjusted$adjusted)
#> [1] 2.5 5.0 7.5 1.0
res$estimates$value[res$estimates$table == "impacts"]
#> [1] 2.5 5.0 7.5 1.0
```

A simulated herd with known impacts (2, 4 and 6 percent) and
interactions (1 point for d1 with d2, 1.5 for d2 with d3) checks the
method. By default
[`simulate_raw_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/simulate_raw_impacts.md)
simulates a crude row for each disease and each combination:

``` r

raw <- simulate_raw_impacts(m, c(d1 = 2, d2 = 4, d3 = 6),
                            interactions = c("d1;d2" = 1, "d2;d3" = 1.5))
raw[, c("disease", "adjusted_for", "value")]
#>   disease adjusted_for     value
#> 1      d1         <NA>  2.714619
#> 2      d2         <NA>  6.136904
#> 3      d3         <NA>  7.116677
#> 4   d1;d2         <NA>  8.608404
#> 5   d2;d3         <NA> 11.772543
deconflate(cm_model(m, raw), n_draws = 0)$adjusted
#>   disease        kind       raw adjusted     change rows
#> 1      d1     disease  2.714619      2.0 -0.2632484    1
#> 2      d2     disease  6.136904      4.0 -0.3482055    1
#> 3      d3     disease  7.116677      6.0 -0.1569099    1
#> 4   d1;d2 interaction  8.608404      1.0         NA    1
#> 5   d2;d3 interaction 11.772543      1.5         NA    1
```

The interaction terms of the model are the combinations that appear in
the `disease` column. The rows must tell the terms apart: there must be
at least as many independent rows as terms (one per disease and per
combination); otherwise
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
stops and names the terms the rows do not determine. Any row that
depends on the interaction helps to determine it: besides a row for the
combination itself, e.g. the impact of d1 among animals with d2
(`among = "d2"`), set against that of d1 overall.

The welfare table of the five-disease example has two interactions
adjusted for all, and a crude LAM;SCK row whose value implies an
antagonistic LAM:SCK interaction:

``` r

welfare <- cm_read_inputs(diseases = file.path(ex, "diseases.csv"),
                          associations = file.path(ex, "associations.csv"),
                          impacts = file.path(ex, "impacts_welfare.csv"))
res_w <- deconflate(welfare, n_draws = 0)
res_w$adjusted
#>   disease        kind  raw   adjusted     change rows
#> 1     LAM     disease 10.0  8.7442977 -0.1255702    1
#> 2     MAS     disease  6.0  3.8319918 -0.3613347    1
#> 3     MET     disease  5.0  3.1977586 -0.3604483    1
#> 4     SCK     disease  3.0  1.0855585 -0.6381472    1
#> 5      RP     disease  4.0  2.1515488 -0.4621128    1
#> 6 LAM;MAS interaction  2.0  2.0000000         NA    1
#> 7  MET;RP interaction  1.5  1.5000000         NA    1
#> 8 LAM;SCK interaction 11.5 -0.6200834         NA    1
res_w$totals
#>    raw_sum adjusted_total interaction_total
#> 1 5.845091       4.147957         0.1398854
```

For a term with one row, `raw` is that row’s value. For a crude
combination row it is the difference between animals with both diseases
and animals with neither, which includes both main impacts: 11.5 points
for LAM;SCK, which the model splits into the main impacts of LAM and SCK
and an interaction of about -0.6.

For event impacts, a combination row adds an interaction on the
log-hazard scale (the culling table of the five-disease example has
one); see
[`vignette("event-impacts")`](https://rasmussenphilip.github.io/deconflate/articles/event-impacts.md).

## Attribution

The expected aggregate is `sum_i p_i b_i + sum_T P(all of T) delta_T`
over the interaction terms T. Removing any disease of a term removes
that term, so each disease’s Shapley value is its own term plus an equal
share of every interaction term it is involved in. This is the
`contributions` table:

``` r

attribute_burden(res)
#>   disease      main interaction     total      share
#> 1      d1 0.1913880   0.0122397 0.2036277 0.09747208
#> 2      d2 0.4860961   0.0122397 0.4983358 0.23854231
#> 3      d3 1.3871243   0.0000000 1.3871243 0.66398561
```

For general loss functions,
[`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
computes Shapley values cell by cell over the joint distribution. For
example, with proportional impacts that compound (here as proportions):

``` r

shapley_by_cell(j, loss_multiplicative(c(d1 = 0.02, d2 = 0.04, d3 = 0.06)))
#>   disease     shapley      share
#> 1      d1 0.001978346 0.09981218
#> 2      d2 0.005922273 0.29879258
#> 3      d3 0.011920065 0.60139525
```

## Which interactions or three-way associations would matter?

[`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
sets an interaction of each size for each pair in turn and reports the
change in the aggregate and in the ranking of diseases;
[`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md)
does the same with an `OR_ratio` row for each triple. Without
interactions, three-way rows change additive results only through
unknown pairs:

``` r

screen_interactions(m, values = c(-1, 1))
#> <cm_screen> 6 scenarios; baseline total 2.109
#>   pair      status   scenario total  change rel_change max_rank_shift rank_corr
#>  d2:d3 interaction delta = -1  2.15  0.0383    0.01815              0         1
#>  d2:d3 interaction  delta = 1  2.07 -0.0383   -0.01815              0         1
#>  d1:d3 interaction  delta = 1  2.09 -0.0209   -0.00990              0         1
#>  d1:d3 interaction delta = -1  2.13  0.0209    0.00990              0         1
#>  d1:d2 interaction  delta = 1  2.09 -0.0201   -0.00955              0         1
#>  d1:d2 interaction delta = -1  2.13  0.0201    0.00955              0         1
#>  failed
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
screen_three_way(m, ratios = c(0.5, 2))
#> <cm_screen> 2 scenarios; baseline total 2.109
#>      pair    status    scenario total    change rel_change max_rank_shift
#>  d1:d2:d3 three-way   ratio = 2  2.11 -6.91e-10  -3.28e-10              0
#>  d1:d2:d3 three-way ratio = 0.5  2.11  1.35e-10   6.40e-11              0
#>  rank_corr failed
#>          1   <NA>
#>          1   <NA>
screen_three_way(mi, ratios = c(0.5, 2))
#> <cm_screen> 2 scenarios; baseline total 2.089
#>      pair    status    scenario total   change rel_change max_rank_shift
#>  d1:d2:d3 three-way   ratio = 2  2.09 -0.00234   -0.00112              0
#>  d1:d2:d3 three-way ratio = 0.5  2.09  0.00231    0.00111              0
#>  rank_corr failed
#>          1   <NA>
#>          1   <NA>
```
