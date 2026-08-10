#!/usr/bin/env Rscript

# Figure 4A, 4E and 4F from a patient-level binary feature matrix.
# Run: Rscript figure4_aef.R --input KM00_Figure4_binary_feature_matrix.tsv --outdir results

required <- c(
  "data.table", "dplyr", "tidyr", "purrr", "tibble", "survival",
  "broom", "ggplot2", "survminer", "cowplot", "ggrepel", "scales"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(survival)
  library(broom)
  library(ggplot2)
  library(survminer)
  library(cowplot)
  library(ggrepel)
  library(scales)
})

parse_args <- function(x) {
  ans <- list(
    input = "KM00_Figure4_binary_feature_matrix.tsv",
    outdir = "results",
    iterations = 100L,
    seed = 123L,
    min_cancer_n = 30L
  )
  if (!length(x)) return(ans)
  if (length(x) %% 2L) stop("Arguments must be --name value pairs.")
  for (i in seq(1L, length(x), 2L)) {
    nm <- sub("^--", "", x[[i]])
    if (!nm %in% names(ans)) stop("Unknown argument: ", x[[i]])
    ans[[nm]] <- x[[i + 1L]]
  }
  ans$iterations <- as.integer(ans$iterations)
  ans$seed <- as.integer(ans$seed)
  ans$min_cancer_n <- as.integer(ans$min_cancer_n)
  ans
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

clinical_columns <- c("KM_SAMPLE_ID", "OS_event", "OS_month", "cancer_type", "onco_level")
plot_exclusions <- c(
  "Dual primary cancer", "Others", "Intestinal cancer",
  "Cholangiocarcinoma and Gallbladder cancer"
)

cancer_map <- c(
  "Adrenal cancer"="AC", "Ampullary carcinoma"="AMPCA", "Breast cancer"="BRCA",
  "Cervical cancer"="CC", "Dual primary cancer"="DPCA", "Esophageal cancer"="ESCA",
  "Gastric cancer"="GC", "Head and Neck cancer"="HNCA",
  "Hepatocellular carcinoma"="HCC", "Kidney Renal cancer"="KIRC",
  "Melanoma"="MEL", "Neuroendocrine tumor"="NET", "Others"="OTHR",
  "Ovarian cancer"="OVCA", "Pancreatic cancer"="PACA", "Prostate cancer"="PRCA",
  "Sarcoma"="SARC", "Thymic carcinoma"="THYM",
  "Bladder and Urethral cancer"="URTC", "Uterine cancer"="UTER",
  "Skin cancer"="SKCA", "Thyroid cancer"="THCA",
  "Peripheral nerve cancer"="PNST", "Testicular cancer"="TECA",
  "Penile cancer"="PSCC", "Vaginal cancer"="VC", "Mesothelioma"="MESO",
  "CNS cancer"="CNST", "Lung cancer"="LUCA",
  "CHOL"="CHOL", "COADREAD"="COADREAD", "GBC"="GBC", "NSCLC"="NSCLC",
  "SBC"="SBC", "SCLC"="SCLC"
)

cancer_colors <- c(
  "Pan-cancer"="#000000", "CNST"="#83b0c4", "HNCA"="#4c6af0",
  "THCA"="#00adfb", "THYM"="#7498ff", "MESO"="#cbcdec",
  "LUCA"="#f07f4b", "NSCLC"="#f07f4b", "SCLC"="#ffe3b9",
  "BRCA"="#fb8072", "ESCA"="#908ccd", "GC"="#d2bae0",
  "COADREAD"="#c76033", "SBC"="#d38736", "HCC"="#005f49",
  "CHOL"="#009462", "GBC"="#4c9798", "AMPCA"="#d1e7e5",
  "PACA"="#00c0ba", "NET"="#fe664f", "KIRC"="#f8766d",
  "AC"="#e58b8c", "URTC"="#d39490", "PRCA"="#5e662b",
  "TECA"="#5d8a00", "PSCC"="#8b9e00", "OVCA"="#a7a400",
  "UTER"="#9dc787", "CC"="#a0cf2f", "VC"="#d9e98f",
  "SARC"="#8462ba", "PNST"="#c07fff", "MEL"="#ffb5df",
  "SKCA"="#ffdedf"
)

abbr <- function(x) {
  y <- unname(cancer_map[x])
  ifelse(is.na(y), x, y)
}

write_tsv <- function(x, filename) {
  data.table::fwrite(x, file.path(args$outdir, filename), sep = "\t", na = "NA")
}

safe_median_survival <- function(time, event) {
  if (length(time) < 2L) return(NA_real_)
  fit <- tryCatch(survfit(Surv(time, event) ~ 1), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  tryCatch(as.numeric(summary(fit)$table[["median"]]), error = function(e) NA_real_)
}

safe_coxph <- function(formula, data = NULL) {
  infinite_warning <- FALSE
  fit <- withCallingHandlers(
    tryCatch(coxph(formula, data = data), error = function(e) NULL),
    warning = function(w) {
      if (grepl("coefficient may be infinite", conditionMessage(w), fixed = TRUE)) {
        infinite_warning <<- TRUE
        invokeRestart("muffleWarning")
      }
    }
  )
  if (is.null(fit) || infinite_warning || any(!is.finite(coef(fit)))) return(NULL)
  fit
}

analyse_feature <- function(dat, feature) {
  x <- as.integer(dat[[feature]])
  ok <- is.finite(dat$OS_month) & !is.na(dat$OS_event) & !is.na(x)
  x <- x[ok]
  time <- dat$OS_month[ok]
  event <- dat$OS_event[ok]
  n_mut <- sum(x == 1L)
  n_wt <- sum(x == 0L)
  if (n_mut < 1L || n_wt < 1L || length(unique(x)) < 2L || sum(event) < 2L) return(NULL)

  fit <- safe_coxph(Surv(time, event) ~ x)
  if (is.null(fit) || !is.finite(coef(fit)[[1]])) return(NULL)
  td <- tryCatch(broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
                 error = function(e) NULL)
  if (is.null(td) || !nrow(td)) return(NULL)

  # Legacy rule: MPI uses the Cox Wald P value. The original column named
  # `cox_pval` was actually a log-rank P value and was used for repeated-
  # significance counts; use an unambiguous public name here.
  lr <- tryCatch(survdiff(Surv(time, event) ~ x), error = function(e) NULL)
  logrank_p <- if (is.null(lr)) NA_real_ else
    pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
  if (!is.finite(logrank_p)) return(NULL)

  hr <- td$estimate[[1]]
  cox_p <- td$p.value[[1]]
  prevalence <- 100 * n_mut / (n_mut + n_wt)
  score_raw <- (-log2(hr)) *
    (-log10(max(cox_p, .Machine$double.xmin))) * prevalence

  tibble(
    feature = feature,
    HR = hr,
    CI_lower = td$conf.low[[1]],
    CI_upper = td$conf.high[[1]],
    p_value = cox_p,
    logrank_p = logrank_p,
    n_mut = n_mut,
    median_os_mut = safe_median_survival(time[x == 1L], event[x == 1L]),
    n_wt = n_wt,
    median_os_wt = safe_median_survival(time[x == 0L], event[x == 0L]),
    prevalence_percent = prevalence,
    score_raw = score_raw
  )
}

zscore_safe <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------
dat <- data.table::fread(args$input, data.table = FALSE, check.names = FALSE)
absent <- setdiff(clinical_columns, names(dat))
if (length(absent)) stop("Missing input columns: ", paste(absent, collapse = ", "))
if (anyDuplicated(dat$KM_SAMPLE_ID)) stop("KM_SAMPLE_ID must be unique.")

feature_columns <- setdiff(names(dat), clinical_columns)
if (!length(feature_columns)) stop("No binary feature columns were found.")
nonbinary <- feature_columns[!vapply(dat[feature_columns], function(x) all(x %in% c(0, 1)), logical(1))]
if (length(nonbinary)) stop("Non-binary feature columns: ", paste(nonbinary, collapse = ", "))

dat <- dat %>%
  mutate(
    KM_SAMPLE_ID = as.character(KM_SAMPLE_ID),
    OS_event = as.integer(OS_event),
    OS_month = as.numeric(OS_month),
    cancer_type = as.character(cancer_type),
    cancer_abbr = abbr(cancer_type)
  ) %>%
  filter(!is.na(KM_SAMPLE_ID), KM_SAMPLE_ID != "", !is.na(OS_month),
         OS_month >= 0, OS_event %in% c(0L, 1L), !is.na(cancer_type), cancer_type != "")

write_tsv(
  dat %>% count(cancer_type, cancer_abbr, name = "n") %>% arrange(desc(n)),
  "cohort_sizes.tsv"
)

# -----------------------------------------------------------------------------
# Figure 4A: cancer-type overall survival curves and cancer-versus-rest HRs
# -----------------------------------------------------------------------------
a_dat <- dat %>% filter(!cancer_type %in% plot_exclusions)

a_long <- bind_rows(
  a_dat %>% transmute(KM_SAMPLE_ID, OS_month, OS_event, group = cancer_abbr),
  a_dat %>% transmute(KM_SAMPLE_ID, OS_month, OS_event, group = "Pan-cancer")
)
a_long$group <- factor(a_long$group, levels = intersect(names(cancer_colors), unique(a_long$group)))
a_fit <- survfit(Surv(OS_month, OS_event) ~ group, data = a_long)
a_test <- survdiff(Surv(OS_month, OS_event) ~ group, data = a_dat %>% mutate(group = cancer_abbr))
a_logp <- pchisq(a_test$chisq, df = length(a_test$n) - 1L, lower.tail = FALSE, log.p = TRUE)
a_p_text <- if (is.finite(a_logp)) sprintf("P = %.2g", exp(a_logp)) else
  sprintf("P < 10^-%d", floor(-a_logp / log(10)))

curve <- ggsurvplot(
  a_fit, data = a_long, risk.table = FALSE, conf.int = FALSE,
  censor = TRUE, palette = unname(cancer_colors[levels(a_long$group)]),
  legend = "none", xlab = "Months", ylab = "Overall survival",
  ggtheme = theme_classic(base_size = 10)
)$plot +
  annotate("text", x = Inf, y = Inf, label = a_p_text, hjust = 1.1, vjust = 1.5,
           parse = FALSE, size = 3)

a_forest <- map_dfr(sort(unique(a_dat$cancer_abbr)), function(ct) {
  tmp <- a_dat %>% mutate(in_cancer = as.integer(cancer_abbr == ct))
  fit <- safe_coxph(Surv(OS_month, OS_event) ~ in_cancer, data = tmp)
  if (is.null(fit)) return(tibble())
  td <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  tibble(cancer_type = ct, HR = td$estimate[[1]], CI_lower = td$conf.low[[1]],
         CI_upper = td$conf.high[[1]], p_value = td$p.value[[1]],
         n = sum(tmp$in_cancer))
}) %>% arrange(HR)
write_tsv(a_forest, "Figure4A_cancer_vs_rest_cox.tsv")

a_forest$cancer_type <- factor(a_forest$cancer_type, levels = a_forest$cancer_type)
forest <- ggplot(a_forest, aes(HR, cancer_type, color = cancer_type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.15) +
  geom_point(size = 2) +
  geom_text(aes(label = format.pval(p_value, digits = 2, eps = 1e-99)),
            x = Inf, hjust = 1, size = 2.5, color = "black") +
  scale_color_manual(values = cancer_colors, guide = "none") +
  scale_x_log10() +
  labs(x = "Hazard ratio", y = NULL) +
  theme_classic(base_size = 9) +
  theme(plot.margin = margin(5.5, 32, 5.5, 5.5))

panel_a <- cowplot::plot_grid(
  curve, forest, nrow = 1, rel_widths = c(1.35, 1), align = "h", axis = "tb"
)
ggsave(file.path(args$outdir, "Figure4A_overall_survival.pdf"), panel_a,
       width = 12, height = 6.2, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure4A_overall_survival.png"), panel_a,
       width = 12, height = 6.2, dpi = 300)

# -----------------------------------------------------------------------------
# Full-cohort feature analysis
# -----------------------------------------------------------------------------
eligible_cancers <- dat %>% count(cancer_type) %>%
  filter(n >= args$min_cancer_n) %>% pull(cancer_type)

full_results <- map_dfr(eligible_cancers, function(ct) {
  tmp <- filter(dat, cancer_type == ct)
  map_dfr(feature_columns, ~ analyse_feature(tmp, .x)) %>%
    mutate(cancer_type = ct, cancer_abbr = abbr(ct), .before = 1)
}) %>%
  mutate(delta_median_OS = median_os_mut - median_os_wt)
write_tsv(full_results, "full_cohort_feature_cox_mpi.tsv")

# -----------------------------------------------------------------------------
# 100 repeated random subsamples: 90% if n < 100, otherwise 80%
# -----------------------------------------------------------------------------
set.seed(args$seed)
score_results <- map_dfr(feature_columns, function(feature_name) {
  map_dfr(eligible_cancers, function(ct) {
    full <- filter(dat, cancer_type == ct)
    sampled_n <- round(ifelse(nrow(full) < 100, 0.9, 0.8) * nrow(full))
    map_dfr(seq_len(args$iterations), function(iteration_number) {
      # The legacy scripts sampled independently for every feature/iteration.
      sampled_ids <- sample(full$KM_SAMPLE_ID, size = sampled_n, replace = FALSE)
      tmp <- filter(full, KM_SAMPLE_ID %in% sampled_ids)
      one <- analyse_feature(tmp, feature_name)
      if (is.null(one) || !nrow(one)) return(tibble())
      one %>% mutate(
        cancer_type = ct,
        cancer_abbr = abbr(ct),
        iteration = iteration_number,
        .before = 1
      )
    })
  })
}) %>%
  group_by(cancer_type, cancer_abbr, iteration) %>%
  mutate(score_z = zscore_safe(score_raw)) %>%
  ungroup()
write_tsv(score_results, "score_z_all_features_100iter.tsv")

mpi_summary <- score_results %>%
  group_by(cancer_type, cancer_abbr, feature) %>%
  summarise(
    median_score_z = median(score_z, na.rm = TRUE),
    n_iterations = sum(!is.na(score_z)),
    n_significant = sum(logrank_p < 0.05, na.rm = TRUE),
    median_HR = median(HR, na.rm = TRUE),
    median_prevalence = median(prevalence_percent, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(mpi_summary, "MPI_100iter_summary.tsv")

# -----------------------------------------------------------------------------
# Figure 4E: median MPI z-score versus mutant-minus-WT median OS
# -----------------------------------------------------------------------------
e_df <- full_results %>%
  select(cancer_type, cancer_abbr, feature, HR, p_value, n_mut, n_wt,
         prevalence_percent, delta_median_OS) %>%
  inner_join(mpi_summary, by = c("cancer_type", "cancer_abbr", "feature")) %>%
  filter(
    !cancer_type %in% plot_exclusions,
    !grepl("1_2", feature, fixed = TRUE),
    is.finite(delta_median_OS), is.finite(median_score_z)
  ) %>%
  mutate(label = paste(cancer_abbr, feature, sep = "\n"))

cor_e <- cor.test(e_df$median_score_z, e_df$delta_median_OS, method = "pearson")
write_tsv(e_df, "Figure4E_plot_data.tsv")
write_tsv(tibble(R = unname(cor_e$estimate), p_value = cor_e$p.value, n = nrow(e_df)),
          "Figure4E_correlation.tsv")

label_df <- e_df %>%
  filter(p_value < 0.05 | n_significant >= 50) %>%
  arrange(desc(abs(median_score_z))) %>%
  slice_head(n = 30)

p_scatter <- ggplot(e_df, aes(median_score_z, delta_median_OS, color = cancer_abbr)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey65") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey65") +
  geom_point(aes(size = n_significant), alpha = 0.65) +
  geom_smooth(method = "lm", se = FALSE, color = "#4C78A8", linewidth = 0.6) +
  ggrepel::geom_text_repel(data = label_df, aes(label = label), size = 2.4,
                           fontface = "italic", max.overlaps = Inf, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R = %.2f", unname(cor_e$estimate)),
           hjust = -0.1, vjust = 1.3, size = 3) +
  scale_color_manual(values = cancer_colors, guide = "none") +
  scale_size_continuous(name = "Significant iterations\n(out of 100)", range = c(1, 5)) +
  labs(x = "Molecular Prognostic Index (median z-score)",
       y = expression(Delta*"OS (Mut - WT), months")) +
  theme_bw(base_size = 10)

p_top <- ggplot(e_df, aes(median_score_z, fill = cancer_abbr)) +
  geom_density(aes(y = after_stat(scaled)), alpha = 0.4, linewidth = 0.15) +
  scale_fill_manual(values = cancer_colors, guide = "none") + theme_void()
p_right <- ggplot(e_df, aes(delta_median_OS, fill = cancer_abbr)) +
  geom_density(aes(y = after_stat(scaled)), alpha = 0.4, linewidth = 0.15) +
  scale_fill_manual(values = cancer_colors, guide = "none") +
  coord_flip() + theme_void()

e_main <- cowplot::plot_grid(
  cowplot::plot_grid(p_top, NULL, rel_widths = c(4, 1)),
  cowplot::plot_grid(p_scatter, p_right, rel_widths = c(4, 1), align = "hv"),
  ncol = 1, rel_heights = c(1, 4)
)
ggsave(file.path(args$outdir, "Figure4E_MPI_OS_correlation.pdf"), e_main,
       width = 9, height = 8, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure4E_MPI_OS_correlation.png"), e_main,
       width = 9, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Figure 4F: recurrent favorable/unfavorable prognostic features
# -----------------------------------------------------------------------------
f_rank <- mpi_summary %>%
  filter(
    !cancer_type %in% plot_exclusions,
    !feature %in% c("onco_level1_2", "onco_level1_2_3")
  ) %>%
  group_by(feature) %>%
  summarise(
    total_significant_iterations = sum(n_significant, na.rm = TRUE),
    significant_cancer_types = sum(n_significant > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_significant_iterations), desc(significant_cancer_types), feature)
selected_features <- head(f_rank$feature, 42)

f_df <- mpi_summary %>%
  filter(
    !cancer_type %in% plot_exclusions,
    feature %in% selected_features,
    n_significant > 0
  ) %>%
  mutate(
    direction = if_else(median_score_z >= 0, "Favorable", "Unfavorable"),
    abs_mpi = abs(median_score_z)
  )

feature_order <- f_rank$feature[f_rank$feature %in% selected_features]
cancer_order <- rev(unique(c(
  "HNCA","THCA","THYM","NSCLC","SCLC","BRCA","ESCA","GC","COADREAD","SBC",
  "HCC","CHOL","GBC","AMPCA","PACA","NET","KIRC","AC","URTC","PRCA",
  "OVCA","UTER","CC","SARC","MEL","SKCA"
)))
f_df <- f_df %>%
  mutate(feature = factor(feature, levels = feature_order),
         cancer_abbr = factor(cancer_abbr, levels = cancer_order))

f_bar <- f_df %>% distinct(feature, cancer_abbr, direction) %>%
  count(feature, direction, name = "n_tumor_types")
write_tsv(f_rank, "Figure4F_feature_rank.tsv")
write_tsv(f_df, "Figure4F_plot_data.tsv")

p_bar <- ggplot(f_bar, aes(feature, n_tumor_types, fill = direction)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = c("Favorable"="#F9B711", "Unfavorable"="#007114")) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = "No. tumor types") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none")

p_dot <- ggplot(f_df, aes(feature, cancer_abbr)) +
  geom_point(aes(size = abs_mpi, color = direction), alpha = 0.9) +
  scale_color_manual(values = c("Favorable"="#F9B711", "Unfavorable"="#007114"),
                     name = "Prognostic significance") +
  scale_size_continuous(name = "Normalized absolute\nMPI score", range = c(1, 6)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_line(linewidth = 0.2, color = "grey85"),
        axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, face = "italic"),
        legend.position = "right")

panel_f <- cowplot::plot_grid(
  p_bar, p_dot, ncol = 1, rel_heights = c(1, 4.5), align = "v", axis = "lr"
)
ggsave(file.path(args$outdir, "Figure4F_recurrent_MPI_features.pdf"), panel_f,
       width = 13, height = 8, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure4F_recurrent_MPI_features.png"), panel_f,
       width = 13, height = 8, dpi = 300)

message("Figure 4A/E/F analysis completed: ", normalizePath(args$outdir))
