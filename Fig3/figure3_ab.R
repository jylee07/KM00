#!/usr/bin/env Rscript

# Figure 3A-B: comparison of gene-level mutation frequencies between
# filtered K-MASTER and GENIE MAF files.
#
# Usage:
# Rscript figure3_ab.R \
#   --km_maf path/to/km_filtered.maf \
#   --genie_maf path/to/genie_filtered.maf \
#   --outdir results

required_packages <- c("data.table", "dplyr", "tidyr", "ggplot2", "qvalue", "ggrepel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(qvalue)
  library(ggrepel)
})

parse_args <- function(x) {
  defaults <- list(
    km_maf = "data/Fig3_KM00.maf",
    genie_maf = "data/Fig3_GENIE.maf",
    outdir = "results"
  )
  if (length(x) == 0) return(defaults)
  if (length(x) %% 2 != 0) stop("Arguments must be supplied as --name value pairs.")
  for (i in seq(1, length(x), by = 2)) {
    key <- sub("^--", "", x[[i]])
    if (!key %in% names(defaults)) stop("Unknown argument: ", x[[i]])
    defaults[[key]] <- x[[i + 1]]
  }
  defaults
}

check_columns <- function(x, required, cohort) {
  absent <- setdiff(required, names(x))
  if (length(absent) > 0) {
    stop(cohort, " MAF is missing required columns: ", paste(absent, collapse = ", "))
  }
}

safe_qvalue <- function(p) {
  p <- as.numeric(p)
  if (length(p) < 2 || all(is.na(p))) return(rep(NA_real_, length(p)))
  tryCatch(
    qvalue::qvalue(p, lambda = 0)$qvalues,
    error = function(e) p.adjust(p, method = "BH")
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

km <- data.table::fread(args$km_maf, data.table = FALSE, check.names = FALSE)
genie <- data.table::fread(args$genie_maf, data.table = FALSE, check.names = FALSE)

check_columns(km, c("SAMPLE.ID", "Tumor_Sample_Barcode", "Hugo_Symbol", "sub"), "K-MASTER")
check_columns(genie, c("Tumor_Sample_Barcode", "Hugo_Symbol", "sub"), "GENIE")

# The supplied K-MASTER MAF has a one-to-one SAMPLE.ID/barcode relationship.
# Enforce it so that a future input cannot silently change the denominator.
km_id_map <- km %>% distinct(SAMPLE.ID, Tumor_Sample_Barcode)
if (anyDuplicated(km_id_map$SAMPLE.ID) || anyDuplicated(km_id_map$Tumor_Sample_Barcode)) {
  stop("K-MASTER SAMPLE.ID and Tumor_Sample_Barcode must have a one-to-one relationship.")
}

km_clean <- km %>%
  transmute(cohort = "K-MASTER", patient_id = as.character(SAMPLE.ID),
            sample_id = as.character(Tumor_Sample_Barcode),
            cancer_type = as.character(sub), gene = as.character(Hugo_Symbol)) %>%
  filter(!is.na(patient_id), patient_id != "", !is.na(cancer_type), cancer_type != "",
         !is.na(gene), gene != "")

genie_clean <- genie %>%
  transmute(cohort = "GENIE", patient_id = as.character(Tumor_Sample_Barcode),
            sample_id = as.character(Tumor_Sample_Barcode),
            cancer_type = as.character(sub), gene = as.character(Hugo_Symbol)) %>%
  filter(!is.na(patient_id), patient_id != "", !is.na(cancer_type), cancer_type != "",
         !is.na(gene), gene != "")

all_maf <- bind_rows(km_clean, genie_clean)

# Legacy Figure 3 denominator: unique IDs represented in each filtered MAF,
# separately for each cancer type.
cohort_sizes <- all_maf %>%
  distinct(cohort, cancer_type, patient_id) %>%
  count(cohort, cancer_type, name = "cohort_n")

# A patient with multiple variants in the same gene is counted once.
mutated_counts <- all_maf %>%
  distinct(cohort, cancer_type, patient_id, gene) %>%
  count(cohort, cancer_type, gene, name = "mutated_n")

frequency_long <- mutated_counts %>%
  left_join(cohort_sizes, by = c("cohort", "cancer_type")) %>%
  mutate(frequency_percent = 100 * mutated_n / cohort_n)

counts_wide <- mutated_counts %>%
  select(cohort, cancer_type, gene, mutated_n) %>%
  pivot_wider(
    names_from = cohort,
    values_from = mutated_n,
    values_fill = 0
  ) %>%
  rename(
    km_mutated_n = `K-MASTER`, genie_mutated_n = GENIE
  )

sizes_wide <- cohort_sizes %>%
  pivot_wider(names_from = cohort, values_from = cohort_n) %>%
  rename(km_n = `K-MASTER`, genie_n = GENIE)

comparison <- counts_wide %>%
  left_join(sizes_wide, by = "cancer_type") %>%
  filter(!is.na(km_n), !is.na(genie_n)) %>%
  mutate(
    km_frequency = 100 * km_mutated_n / km_n,
    genie_frequency = 100 * genie_mutated_n / genie_n,
    km_wildtype_n = km_n - km_mutated_n,
    genie_wildtype_n = genie_n - genie_mutated_n,
    frequency_difference = km_frequency - genie_frequency,
    p_value = mapply(
      function(a, b, c, d) suppressWarnings(chisq.test(matrix(c(a, b, c, d), nrow = 2))$p.value),
      km_mutated_n, genie_mutated_n, km_wildtype_n, genie_wildtype_n
    )
  ) %>%
  group_by(cancer_type) %>%
  mutate(q_value = safe_qvalue(p_value)) %>%
  ungroup() %>%
  mutate(
    log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
    log10_q = -log10(pmax(q_value, .Machine$double.xmin)),
    enrichment = case_when(
      p_value < 0.05 & frequency_difference > 5 ~ "K-MASTER enriched",
      p_value < 0.05 & frequency_difference < -5 ~ "GENIE enriched",
      TRUE ~ "Not significant"
    ),
    label_a = if_else(km_frequency > 25 | genie_frequency > 25, gene, NA_character_),
    label_b = if_else(enrichment != "Not significant",
                      paste0(cancer_type, "-", gene), NA_character_)
  )

data.table::fwrite(cohort_sizes, file.path(args$outdir, "cohort_sizes.tsv"), sep = "\t")
data.table::fwrite(frequency_long, file.path(args$outdir, "mutation_frequency_long.tsv"), sep = "\t")
data.table::fwrite(comparison, file.path(args$outdir, "mutation_frequency_comparison.tsv"), sep = "\t")

correlation <- cor.test(comparison$km_frequency, comparison$genie_frequency,
                        method = "pearson")
correlation_table <- data.frame(
  method = "Pearson",
  n_gene_cancer_pairs = nrow(comparison),
  estimate = unname(correlation$estimate),
  p_value = correlation$p.value,
  conf_low = correlation$conf.int[[1]],
  conf_high = correlation$conf.int[[2]]
)
data.table::fwrite(correlation_table, file.path(args$outdir, "figure3A_correlation.tsv"), sep = "\t")

set.seed(1)
plot_a_data <- comparison[sample(seq_len(nrow(comparison))), ]
plot_a <- ggplot(plot_a_data, aes(km_frequency, genie_frequency)) +
  geom_abline(slope = 1, intercept = 0, color = "grey55", linewidth = 0.5) +
  geom_point(aes(size = km_n, color = cancer_type), alpha = 0.8,
             position = position_jitter(width = 1.2, height = 1.2)) +
  ggrepel::geom_text_repel(aes(label = label_a), na.rm = TRUE,
                           size = 3, fontface = "italic", color = "black",
                           max.overlaps = 30, show.legend = FALSE) +
  annotate("text", x = 3, y = 97, hjust = 0,
           label = sprintf("Pearson R = %.2f\nP = %.2g", correlation_table$estimate,
                           correlation_table$p_value)) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
  scale_size_continuous(name = "K-MASTER cohort size") +
  labs(x = "K-MASTER mutation frequency (%)",
       y = "GENIE mutation frequency (%)", color = "Cancer type") +
  theme_bw(base_size = 11)

plot_b <- ggplot(comparison, aes(frequency_difference, log10_p)) +
  geom_vline(xintercept = c(-5, 5), color = "grey50", linetype = "longdash") +
  geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "longdash") +
  geom_point(aes(size = log10_q, color = enrichment), alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = label_b), na.rm = TRUE,
                           size = 2.8, fontface = "italic", color = "black",
                           max.overlaps = 100, show.legend = FALSE) +
  scale_color_manual(values = c("K-MASTER enriched" = "#AB2A30",
                                "GENIE enriched" = "#33B977",
                                "Not significant" = "grey75")) +
  labs(x = "Mutation-frequency difference (K-MASTER - GENIE, percentage points)",
       y = expression(-log[10](italic(P))), color = NULL,
       size = expression(-log[10](italic(q)))) +
  theme_classic(base_size = 11)

ggsave(file.path(args$outdir, "Figure3A_mutation_frequency_correlation.pdf"),
       plot_a, width = 8.2, height = 6.5, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure3A_mutation_frequency_correlation.png"),
       plot_a, width = 8.2, height = 6.5, dpi = 300)
ggsave(file.path(args$outdir, "Figure3B_mutation_frequency_difference.pdf"),
       plot_b, width = 10, height = 7, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure3B_mutation_frequency_difference.png"),
       plot_b, width = 10, height = 7, dpi = 300)

message("Completed Figure 3A-B analysis. Results: ", normalizePath(args$outdir))
