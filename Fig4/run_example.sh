#!/usr/bin/env bash
set -euo pipefail

Rscript figure4_aef.R \
  --input KM00_Figure4_binary_feature_matrix.tsv \
  --outdir results \
  --iterations 100 \
  --seed 123 \
  --min_cancer_n 30
