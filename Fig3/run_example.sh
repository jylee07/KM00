#!/usr/bin/env bash
set -euo pipefail

Rscript figure3_ab.R \
  --km_maf "data/km00_final_hyper_commongene(final_analysis).maf" \
  --genie_maf "data/GENIE_all_white_filter_hypermutated_sub_commongene(final).maf" \
  --outdir results
