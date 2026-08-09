cran_packages <- c(
  "dplyr", "tidyr", "readr", "purrr", "tibble", "stringr",
  "ggplot2", "survival", "survminer", "broom", "glmnet",
  "ggridges", "forcats", "patchwork", "scales", "yaml"
)
missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)
if (!requireNamespace("qvalue", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("qvalue", ask = FALSE, update = FALSE)
}
