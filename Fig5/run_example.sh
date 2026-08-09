#!/usr/bin/env bash
set -euo pipefail

Rscript figure5_cfgjkl.R \
  --kmaster Figure5_KMASTER_modeling_matrix.tsv \
  --genie Figure5_GENIE_validation_matrix.tsv \
  --outdir results \
  --seed 42 \
  --stability_iterations 100
