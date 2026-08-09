# Figure 5C, 5F, 5G, 5J, 5K and 5L

The public pipeline starts from two patient-level tab-separated matrices:

- `Figure5_KMASTER_modeling_matrix.tsv`
- `Figure5_GENIE_validation_matrix.tsv`

No MAF, RDS model, event-level clinical table, or previously summarized association
table is required by the public analysis script.

## Identifiers and endpoints

The K-MASTER identifier is `KM_SAMPLE_ID`. The external identifier is
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
`broom`, `glmnet`, `ggplot2`, `ggridges`, `ggrepel`, `survminer`, and `patchwork`.

## Panels

- Figure 5C: cancer-specific TTNT ridgeline distributions, capped at 70 months only
  for visualization.
- Figure 5F: cancer-feature univariable Cox associations for composite treatment
  durability.
- Figure 5G: per-cancer elastic-net stability across 100 bootstrap resamples.
- Figure 5J: internal 30% test-set treatment-durability Kaplan-Meier plot.
- Figure 5K: external GENIE-BPC NSCLC treatment-durability validation.
- Figure 5L: OS plots for the internal test and external validation cohorts.

## Integrated model

The 70/30 split is stratified by cancer type and uses seed 42. Alpha and lambda are
selected using only cross-validation within the training cohort. The test cohort is not
used for model or hyperparameter selection. The risk-group cutoff is the median training
linear predictor and is applied unchanged to the internal test and external validation
cohorts.

For GENIE-BPC, `NSCLC` is mapped to the corresponding K-MASTER `LUCA` design-matrix
level during prediction; the displayed cohort label remains `NSCLC`.

## OS availability

K-MASTER OS was matched from the survival clinical table for 1,190 of 1,305 patients.
GENIE-BPC OS was matched for 640 of 696 patients using the exact first-line platinum
regimen and its corresponding regimen-specific OS fields. Patients without an exact
regimen-specific OS column remain `NA` and are excluded only from Figure 5L.

## Important methodological note

This implementation removes test-set leakage present in one exploratory legacy script:
hyperparameters and score orientation are determined from training data only. Therefore,
exact risk scores or group assignments may differ from an earlier figure produced by
test-set-guided model selection. The public code represents the leakage-free analysis.
