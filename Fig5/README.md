# Figure 5C, 5F, 5G, 5J, 5K and 5L

The public pipeline starts from two patient-level tab-separated matrices:

- `Figure5_KMASTER_modeling_matrix.tsv`
- `Figure5_GENIE_validation_matrix.tsv`

No MAF, RDS model, event-level clinical table, or previously summarized association
table is required by the public analysis script.

## Identifiers and endpoints

The K-MASTER identifier is `KM_ID`. The external identifier is
`GENIE_PATIENT_ID`.

`TTNT_status` uses:

- `0`: censored
- `1`: second-line treatment observed
- `2`: death before second-line treatment

The plotted treatment-durability endpoint treats status 1 or 2 as an event. `OS_event`
uses 0 for living/censored and 1 for deceased.

## Run

```bash
Rscript figure5_cfgjkl.R \
  --kmaster Figure5_KMASTER_modeling_matrix.tsv \
  --genie Figure5_GENIE_validation_matrix.tsv \
  --outdir results \
  --seed 42 \
  --stability_iterations 100
```

Required packages: `data.table`, `dplyr`, `tidyr`, `purrr`, `tibble`, `survival`,
`glmnet`, `ggplot2`, `ggridges`, `ggrepel`, `survminer`, and `patchwork`.

## Panels

- Figure 5C: cancer-specific TTNT ridgeline distributions, capped at 70 months only
  for visualization.
- Figure 5F: cancer-feature TTNT comparisons using the Wilcoxon test and median TTNT
  difference, matching the uploaded legacy association script. The manuscript legend
  calls this association plot Figure 5G; the filename retains the requested panel name.
- Figure 5G: per-cancer elastic-net stability across 100 bootstrap resamples. Every
  bootstrap iteration compares alpha 0.3-0.9 using up to five-fold CV and uses
  `lambda.min`. The manuscript legend calls this coefficient plot Figure 5H.
- Figure 5J: internal 30% test-set treatment-durability Kaplan-Meier plot.
- Figure 5K: external GENIE-BPC NSCLC treatment-durability validation.
- Figure 5L: OS plots for the internal test and external validation cohorts.

## Integrated model

The 70/30 split is stratified by cancer type and uses seed 42. The integrated pan-cancer
model compares alpha 0.3-0.9 by **10-fold CV**, because the uploaded original script
uses 10 folds when the training set has at least 100 patients (and five folds only below
100). It uses `lambda.min`. The risk-group cutoff is the training median and is applied
unchanged to the internal test and external validation cohorts.

For GENIE-BPC, `NSCLC` is mapped to the corresponding K-MASTER `LUCA` design-matrix
level during prediction; the displayed cohort label remains `NSCLC`.

## OS availability

K-MASTER OS was matched from the survival clinical table for 1,190 of 1,305 patients.
GENIE-BPC OS was matched for 640 of 696 patients using the exact first-line platinum
regimen and its corresponding regimen-specific OS fields. Patients without an exact
regimen-specific OS column remain `NA` and are excluded only from Figure 5L.

## Important methodological note

The uploaded legacy file contained several successively overwritten exploratory KM
definitions, including hard-coded cutoffs and a test-outcome-based score flip. Those lines
cannot all represent one final analysis. This public script preserves the documented
split, alpha grid, fold rule and `lambda.min`, while locking score direction and cutoff
without consulting outcomes in the test or validation cohorts. Consequently, group
assignments may differ from a plot generated from one of the exploratory hard-coded
branches. The supplied legacy score table can be retained as an audit reference, but it
is not required as an input.

The simple matrices contain clinical variables, TMB, hotspot and gene-level features.
They do not contain pathway or MSI columns. Therefore the public association/modeling
code analyzes the features actually present in the matrices; pathway/MSI points require
those binary columns to be appended with identical names to both matrices.
