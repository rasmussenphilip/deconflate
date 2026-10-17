# Pig example: PRRS and the respiratory disease complex

An illustrative population of growing pigs with four respiratory diseases
and three impact tables. **The values are made up to illustrate the
package, not taken from studies.** `vignette("other-settings")` works
through it.

```r
dir <- system.file("extdata", "pigs", package = "deconflate")
# One model per impact table; the disease and association tables are shared
adg <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                      associations = file.path(dir, "associations.csv"),
                      impacts = file.path(dir, "impacts_adg.csv"),
                      outcome = "average daily gain loss", units = "g/day")
deconflate(adg)
```

Diseases: PRRS (porcine reproductive and respiratory syndrome), PCV2
(porcine circovirus type 2), MHP (*Mycoplasma hyopneumoniae*, enzootic
pneumonia), SIV (swine influenza). All values refer to the grow-finish
period.

## Files

| File | What it shows |
|---|---|
| `diseases.csv` | Prevalences, three of them uncertain (`beta`, `uniform`, `pert`). |
| `associations.csv` | Odds ratios for five of the six pairs; PCV2 and SIV have no row, so that pair is unknown and the fit fills it in (`$unknown_pairs`). |
| `impacts_adg.csv` | Loss of average daily gain (g/day), one crude row per disease, and the PRRS;MHP interaction (adjusted for `all`). |
| `impacts_environment.csv` | Increase in emission intensity (% of kg CO2-eq per kg live weight), one crude row per disease. |
| `impacts_welfare.csv` | Points lost on a 0-100 welfare index, one crude row per disease. |

The outcome and units are given as arguments when a table is read
(optional); the first row's note of each impact table says what they are.
Every table has a `note` column describing its rows; the association and
impact tables also have a `study` column. Uncertain values have a
distribution in `dist` and `p1`-`p4`; rows with an empty `dist` are point
values.
