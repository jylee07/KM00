#!/usr/bin/env bash
set -euo pipefail

Rscript figure3_ab.R \
  --km_maf "data/Fig3_KM00.maf" \
  --genie_maf "data/Fig3_GENIE.maf" \
  --outdir results
