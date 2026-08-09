# Figure 4A, 4E and 4F

This pipeline starts from the public patient-level binary feature matrix. It does not
read the original MAF, `feature_defs.rds`, clinical source table, or any previously
summarized Cox/MPI result.

## Input

`KM00_Figure4_binary_feature_matrix.tsv` contains one row per patient and uses
`KM_SAMPLE_ID` as the only sample-ID column. Required clinical columns are:

- `KM_SAMPLE_ID`
- `OS_event` (`0` = censored, `1` = death)
- `OS_month`
- `cancer_type`
- `onco_level`

All remaining columns are treated as binary features and must contain only `0` or `1`.

## Run

```bash
Rscript figure4_aef.R \
  --input KM00_Figure4_binary_feature_matrix.tsv \
  --outdir results \
  --iterations 100 \
  --seed 123 \
  --min_cancer_n 30
```

The repeated Cox analysis can take several hours because it fits hundreds of thousands
of models. The script writes the complete iteration-level result so it does not need to
be rerun merely to inspect or validate the summary.

## Analysis

### Figure 4A

- Kaplan-Meier overall-survival curves are estimated by cancer type.
- The black curve is the pan-cancer cohort.
- The adjacent forest plot compares each cancer type with all remaining cancers using
  a univariable Cox model.

### Figures 4E/F

- Cancer types with at least 30 patients are analyzed.
- Each patient-level feature is compared as present (`1`) versus absent (`0`).
- Feature prevalence is the percentage of feature-positive patients.
- The Molecular Prognostic Index follows the legacy calculation and is

  `-log2(Cox HR) * -log10(log-rank P) * prevalence_percent`.

  The original workflow used two different P values: log-rank P for the MPI magnitude
  and Cox Wald P for counting significant resampling iterations. Both are retained in
  the output as `p_value` and `cox_pval`, respectively.

- For each cancer type, 100 random subsamples are drawn without replacement.
- A 90% subsample is used when the cancer cohort has fewer than 100 patients; otherwise
  an 80% subsample is used. This is repeated subsampling, not bootstrap resampling.
- MPI values are z-score normalized across successfully fitted features within each
  cancer type and iteration using the sample standard deviation.
- Figure 4E compares the median repeated-subsampling MPI z-score with the full-cohort
  mutant-minus-WT median OS difference. Point size is the number of iterations with
  Cox Wald `P < 0.05`.
- Figure 4F ranks features by their total number of significant iterations and displays
  the leading 42 features. The upper bar is the number of cancer types in which the
  feature was significant in at least one iteration.

## Main outputs

- `Figure4A_overall_survival.pdf/png`
- `Figure4E_MPI_OS_correlation.pdf/png`
- `Figure4F_recurrent_MPI_features.pdf/png`
- `full_cohort_feature_cox_mpi.tsv`
- `score_z_all_features_100iter.tsv`
- `MPI_100iter_summary.tsv`
- Panel-specific plot data and validation tables

## Reproducibility note

The seed is set once immediately before the nested cancer/iteration sampling loop.
Changing cancer order, feature order, seed, or input rows can therefore change the exact
subsamples. Preserve the supplied matrix row and column order for exact reruns.
