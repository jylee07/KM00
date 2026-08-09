# K-MASTER pan-cancer clinicogenomic analysis

This repository contains the analysis code accompanying the manuscript:

> **Molecular determinants of prognosis and treatment durability across cancers:  
> a nationwide clinicogenomic study of the K-MASTER program**

The repository provides code for the principal genomic, survival, platinum-treatment,
and immune-checkpoint-blockade analyses presented in Figures 2–6.

## Repository structure

| Folder | Panels or analysis | Main script(s) |
|---|---|---|
| `Fig2/` | Mutation-frequency comparison across cancer functional categories | `cancer_category_mutation_frequency.py` |
| `Fig3/` | Figure 3A–B: K-MASTER–GENIE mutation-frequency correlation and frequency differences | `figure3_ab.R` |
| `Fig4/` | Figure 4A, E and F: overall survival and Molecular Prognostic Index analyses | `figure4_aef.R` |
| `Fig5/` | Figure 5C, F, G, J, K and L: platinum TTNT, prognostic associations, elastic-net stability, and validation | `figure5_cfgjkl.R` |
| `Fig6/` | ICB pan-cancer and cancer-specific TTNT analyses | `icb_pancancer_hr.py`, `icb_percancer_hr.py` |

Panels consisting exclusively of study schematics, cohort diagrams, or manually
assembled graphical summaries are not generated programmatically. Each figure folder
contains a separate README describing its inputs, statistical procedures, and outputs.

## Data availability

The de-identified input data required by these scripts are available through Figshare:

> **Figshare record:** [DOI and public URL to be added]

Download the Figshare archive and place the files under the repository-level `Data/`
directory while preserving the directory structure described below. The data files are
not duplicated in this GitHub repository.

All identifiers in the released analysis data are study-specific identifiers, including
`KM_ID` and `KM_SAMPLE_ID`; directly identifying information is not included.

### Expected data files

| Figure | Figshare input | Description |
|---|---|---|
| Figure 2 | `KM00_Supple1_Cohortinfo.xlsx` | Subject-level cohort and cancer-type information |
| Figure 2 | `KM00_Supple2_NGSinfo.xlsx` | NGS sample metadata and selected-sample indicator |
| Figure 2 | `KM00_Panel_GeneList.csv` | Gene coverage for each sequencing panel |
| Figure 2 | `KM00_AllSamples_Combined_Mutations.maf` | Combined K-MASTER mutation calls |
| Figure 3 | `KMASTER_filtered_common_gene.maf` | Filtered K-MASTER MAF used for cohort comparison |
| Figure 3 | `GENIE_filtered_common_gene.maf` | Filtered AACR Project GENIE MAF |
| Figure 4 | `KM00_Figure4_binary_feature_matrix.tsv` | Patient-level survival and binary molecular-feature matrix |
| Figure 5 | `Figure5_KMASTER_modeling_matrix.tsv` | K-MASTER platinum TTNT, OS, clinical, and molecular features |
| Figure 5 | `Figure5_GENIE_validation_matrix.tsv` | GENIE-BPC NSCLC external-validation matrix |
| Figure 6 | `KM00_Supple7-1_ICB.xlsx` | Subject-level ICB treatment and TTNT information |
| Figure 6 | `KM00_Supple7-2_ICB_GenomicFeatures.xlsx` | Subject-level ICB genomic features |
| Figure 6 | `KM00_Panel_GeneList.csv` | Panel-specific gene coverage |
| Figure 6 | `KM00_AllSamples_Combined_Mutations.maf` | Combined mutation calls |

The exact filenames in the Figshare release should be used in the figure-specific
configuration or command-line arguments. Previously summarized figure-result tables are
not used as primary inputs for Figures 3–5: statistical testing, Cox modeling, repeated
resampling, MPI calculation, and figure generation are performed by the supplied code.

## Data organization

After downloading the Figshare archive, the recommended layout is:

```text
Data/
├── Fig2/
│   ├── KM00_Supple1_Cohortinfo.xlsx
│   ├── KM00_Supple2_NGSinfo.xlsx
│   ├── KM00_Panel_GeneList.csv
│   └── KM00_AllSamples_Combined_Mutations.maf
├── Fig3/
│   ├── KMASTER_filtered_common_gene.maf
│   └── GENIE_filtered_common_gene.maf
├── Fig4/
│   └── KM00_Figure4_binary_feature_matrix.tsv
├── Fig5/
│   ├── Figure5_KMASTER_modeling_matrix.tsv
│   └── Figure5_GENIE_validation_matrix.tsv
└── Fig6/
    ├── KM00_Supple1_Cohortinfo.xlsx
    ├── KM00_Supple2_NGSinfo.xlsx
    ├── KM00_Supple7-1_ICB.xlsx
    ├── KM00_Supple7-2_ICB_GenomicFeatures.xlsx
    ├── KM00_Panel_GeneList.csv
    └── KM00_AllSamples_Combined_Mutations.maf
```

Files shared by Figures 2 and 6 may instead be stored once in `Data/shared/` if the
corresponding paths in the Python scripts are updated.

## Software environment

Figure 2 and Figure 6 analyses are implemented in Python. Figure 3–5 analyses are
implemented in R.

### Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### R

The principal R packages are:

- `data.table`, `dplyr`, `tidyr`, `purrr`, and `tibble`
- `ggplot2`, `ggrepel`, `ggridges`, `patchwork`, and `cowplot`
- `survival`, `survminer`, `broom`, and `glmnet`
- `qvalue`

Install any missing packages before running the scripts. The R and package versions used
for the final analysis should be recorded in `sessionInfo.txt`.

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

## Reproducing the analyses

Run all commands from the repository root.

### Figure 2

```bash
python Fig2/cancer_category_mutation_frequency.py
```

The script compares mutation frequencies across cancer functional categories using
Fisher's exact tests and Benjamini–Hochberg false-discovery-rate correction.

### Figure 3A–B

```bash
Rscript Fig3/figure3_ab.R \
  --km_maf Data/Fig3/KMASTER_filtered_common_gene.maf \
  --genie_maf Data/Fig3/GENIE_filtered_common_gene.maf \
  --outdir results/Fig3
```

Mutation frequency is calculated after collapsing multiple variants in the same patient
and gene. The two cohorts are compared across matched cancer-type–gene pairs.

### Figure 4A, E and F

```bash
Rscript Fig4/figure4_aef.R \
  --input Data/Fig4/KM00_Figure4_binary_feature_matrix.tsv \
  --outdir results/Fig4 \
  --iterations 100 \
  --seed 123 \
  --min_cancer_n 30
```

The Molecular Prognostic Index (MPI) is calculated as:

```text
MPI = -log2(Cox HR) × -log10(log-rank P) × feature prevalence (%)
```

Figure 4 stability analysis uses repeated random subsampling without replacement. A 90%
subsample is used for cancer cohorts with fewer than 100 patients and an 80% subsample
for larger cohorts, repeated 100 times. Cox Wald P is retained separately and is used to
count significant iterations.

### Figure 5C, F, G, J, K and L

```bash
Rscript Fig5/figure5_cfgjkl.R \
  --kmaster Data/Fig5/Figure5_KMASTER_modeling_matrix.tsv \
  --genie Data/Fig5/Figure5_GENIE_validation_matrix.tsv \
  --outdir results/Fig5 \
  --seed 42 \
  --stability_iterations 100
```

`TTNT_status` is coded as `0` for censoring, `1` for observed second-line treatment, and
`2` for death before second-line treatment. Figures 5J and 5K use the composite event of
second-line treatment or death before second-line treatment.

For the integrated Cox elastic-net model, the 70:30 split is stratified by cancer type.
Alpha values 0.3-0.9 are compared by 10-fold cross-validation (the original conditional
rule uses five folds only for training sets smaller than 100), and `lambda.min` is used.
Feature selection, hyperparameter selection, score direction, and the risk-group cutoff
are determined using the training cohort only. The internal test and external GENIE-BPC
cohorts are not used during model selection.

The feature-association panel uses a Wilcoxon comparison and the median TTNT difference,
as in the legacy analysis. Figure 5G uses 100 bootstrap resamples with replacement; alpha
0.3-0.9 and `lambda.min` are selected using up to five-fold CV inside each iteration. This
differs from the repeated subsampling without replacement used in Figure 4.

### Figure 6

```bash
python Fig6/icb_pancancer_hr.py
python Fig6/icb_percancer_hr.py
```

The first script performs pan-cancer Cox proportional-hazards analyses of TMB, MSI,
MATH, pathway/gene mutation status, mutational-signature presence, cancer type, and ICB
treatment line against TTNT. The second script fits the same features separately within
each cancer type with at least 20 ICB-treated patients. Benjamini–Hochberg correction is
applied to the corresponding analysis family.

## Outputs

Each script writes statistical source data and figure files to the specified output
directory. Generated outputs are not required as inputs and do not need to be committed
to the repository. A typical output structure is:

```text
results/
├── Fig2/
├── Fig3/
├── Fig4/
├── Fig5/
└── Fig6/
```

The repeated Cox and elastic-net analyses can require several hours, depending on the
hardware and number of resampling iterations. Simpler plotting and frequency analyses
generally complete within several minutes.

## External datasets

The analyses use the following external resources where indicated:

- AACR Project GENIE for mutation-frequency comparison;
- GENIE-BPC NSCLC for external validation of the platinum-treatment model.

The exact data releases and access dates should match those reported in the manuscript
Methods and Data Availability sections.

## Citation

If you use this code or the accompanying Figshare data release, please cite the
corresponding manuscript and the Figshare record.

## Contact

For questions regarding the analysis code or data access, please open a GitHub issue or
contact the corresponding authors listed in the manuscript.
