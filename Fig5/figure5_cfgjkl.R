#!/usr/bin/env Rscript

# Figure 5C, 5F, 5G, 5J, 5K and 5L from two patient-level TSV matrices.

required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "tibble", "survival",
  "glmnet", "ggplot2", "ggridges", "ggrepel", "survminer", "patchwork"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(survival); library(glmnet); library(ggplot2)
  library(ggridges); library(ggrepel); library(survminer); library(patchwork)
})

parse_args <- function(x) {
  out <- list(
    kmaster = "Figure5_KMASTER_modeling_matrix.tsv",
    genie = "Figure5_GENIE_validation_matrix.tsv",
    outdir = "results",
    seed = 42L,
    stability_iterations = 100L
  )
  if (!length(x)) return(out)
  if (length(x) %% 2L) stop("Arguments must be --name value pairs.")
  for (i in seq(1L, length(x), 2L)) {
    nm <- sub("^--", "", x[[i]])
    if (!nm %in% names(out)) stop("Unknown argument: ", x[[i]])
    out[[nm]] <- x[[i + 1L]]
  }
  out$seed <- as.integer(out$seed)
  out$stability_iterations <- as.integer(out$stability_iterations)
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, filename) {
  data.table::fwrite(x, file.path(args$outdir, filename), sep = "\t", na = "NA")
}

km_meta <- c(
  "KM_SAMPLE_ID", "TTNT_month", "TTNT_status", "OS_month", "OS_event",
  "cancer_type", "female", "late"
)
genie_meta <- c(
  "GENIE_PATIENT_ID", "TTNT_month", "TTNT_status", "OS_month", "OS_event",
  "cancer_type", "female", "late", "platinum_regimen", "OS_mapping"
)

km <- data.table::fread(args$kmaster, data.table = FALSE, check.names = FALSE)
genie <- data.table::fread(args$genie, data.table = FALSE, check.names = FALSE)
if (length(setdiff(km_meta, names(km)))) stop("K-MASTER matrix is missing required columns.")
if (length(setdiff(genie_meta, names(genie)))) stop("GENIE matrix is missing required columns.")
if (anyDuplicated(km$KM_SAMPLE_ID)) stop("KM_SAMPLE_ID must be unique.")
if (anyDuplicated(genie$GENIE_PATIENT_ID)) stop("GENIE_PATIENT_ID must be unique.")

feature_columns <- setdiff(names(km), km_meta)
if (!identical(feature_columns, setdiff(names(genie), genie_meta))) {
  stop("K-MASTER and GENIE feature columns or their order differ.")
}
bad_binary <- feature_columns[
  !vapply(km[feature_columns], function(z) all(z %in% c(0, 1)), logical(1)) |
  !vapply(genie[feature_columns], function(z) all(z %in% c(0, 1)), logical(1))
]
if (length(bad_binary)) stop("Non-binary feature columns: ", paste(bad_binary, collapse = ", "))

km <- km %>% mutate(
  KM_SAMPLE_ID = as.character(KM_SAMPLE_ID), TTNT_month = as.numeric(TTNT_month),
  TTNT_status = as.integer(TTNT_status), OS_month = as.numeric(OS_month),
  OS_event = as.integer(OS_event), cancer_type = as.character(cancer_type),
  event_comp = as.integer(TTNT_status %in% c(1L, 2L))
)
genie <- genie %>% mutate(
  GENIE_PATIENT_ID = as.character(GENIE_PATIENT_ID), TTNT_month = as.numeric(TTNT_month),
  TTNT_status = as.integer(TTNT_status), OS_month = as.numeric(OS_month),
  OS_event = as.integer(OS_event), cancer_type = as.character(cancer_type),
  event_comp = as.integer(TTNT_status %in% c(1L, 2L))
)

cancer_colors <- c(
  "LUCA"="#f07f4b", "OVCA"="#a7a400", "CHOL"="#009462", "URTC"="#d39490",
  "HNCA"="#4c6af0", "GC"="#d2bae0", "GBC"="#4c9798", "BRCA"="#fb8072",
  "UTER"="#9dc787", "ESCA"="#908ccd", "CC"="#a0cf2f", "NET"="#fe664f",
  "SARC"="#8462ba", "AMPCA"="#d1e7e5", "HCC"="#005f49", "PACA"="#00c0ba",
  "PRCA"="#5e662b", "KIRC"="#f8766d", "THCA"="#00adfb", "SCLC"="#ffe3b9"
)

# -----------------------------------------------------------------------------
# Figure 5C: TTNT distributions by cancer type
# -----------------------------------------------------------------------------
c_dat <- km %>%
  filter(
    !is.na(TTNT_month), !is.na(cancer_type),
    !cancer_type %in% c("Intestinal cancer", "Cholangiocarcinoma and Gallbladder cancer")
  ) %>%
  group_by(cancer_type) %>% filter(n() >= 5L) %>% ungroup() %>%
  mutate(TTNT_cap = pmin(TTNT_month, 70))
c_order <- c_dat %>% group_by(cancer_type) %>%
  summarise(median_TTNT = median(TTNT_month), .groups = "drop") %>%
  arrange(median_TTNT) %>% pull(cancer_type)
c_dat$cancer_type <- factor(c_dat$cancer_type, levels = c_order)
c_summary <- c_dat %>% group_by(cancer_type) %>% summarise(
  n = n(), median_TTNT_month = median(TTNT_month),
  Q1 = quantile(TTNT_month, 0.25), Q3 = quantile(TTNT_month, 0.75), .groups = "drop"
)
write_tsv(c_dat, "Figure5C_plot_data.tsv")
write_tsv(c_summary, "Figure5C_summary.tsv")

p_c <- ggplot(c_dat, aes(TTNT_cap, cancer_type, fill = cancer_type)) +
  geom_density_ridges(scale = 0.8, rel_min_height = 0.01, alpha = 0.8,
                      linewidth = 0.3, bandwidth = 2.5) +
  geom_point(data = c_summary, aes(median_TTNT_month, cancer_type),
             inherit.aes = FALSE, shape = 124, size = 5, color = "red") +
  coord_cartesian(xlim = c(0, 70)) +
  labs(x = "TTNT (months; capped at 70 for visualization)", y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "none", panel.grid.minor = element_blank())
ggsave(file.path(args$outdir, "Figure5C_TTNT_distribution.pdf"), p_c,
       width = 10, height = 8, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure5C_TTNT_distribution.png"), p_c,
       width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Figure 5F: original-style feature association analysis
# The legacy script compared TTNT distributions with a Wilcoxon test and plotted
# the median TTNT difference. It did not fit a univariable Cox model for this panel.
# -----------------------------------------------------------------------------
fit_univariable <- function(dat, feature) {
  z <- as.integer(dat[[feature]])
  n_pos <- sum(z == 1L); n_neg <- sum(z == 0L)
  if (n_pos < 5L || n_neg < 5L || length(unique(z)) < 2L || sum(dat$event_comp) < 3L) {
    return(NULL)
  }
  ttnt_pos <- dat$TTNT_month[z == 1L]
  ttnt_neg <- dat$TTNT_month[z == 0L]
  p <- tryCatch(wilcox.test(ttnt_pos, ttnt_neg, exact = FALSE)$p.value,
                error = function(e) NA_real_)
  tibble(
    feature = feature, median_positive = median(ttnt_pos, na.rm = TRUE),
    median_negative = median(ttnt_neg, na.rm = TRUE), p_value = p,
    n_positive = n_pos, n_negative = n_neg,
    prevalence_percent = 100 * n_pos / (n_pos + n_neg)
  )
}

eligible_cancers <- km %>% count(cancer_type) %>% filter(n >= 30L) %>% pull(cancer_type)
f_results <- map_dfr(eligible_cancers, function(ct) {
  dat <- filter(km, cancer_type == ct)
  map_dfr(c("female", "late", feature_columns), ~ fit_univariable(dat, .x)) %>%
    mutate(cancer_type = ct, .before = 1)
}) %>% mutate(
  median_difference_days = 30.4375 * (median_positive - median_negative),
  log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
  significant = p_value < 0.05,
  label = paste(cancer_type, feature, sep = "_")
)
write_tsv(f_results, "Figure5F_TTNT_Wilcoxon.tsv")
f_labels <- f_results %>% filter(significant) %>% arrange(p_value) %>% slice_head(n = 35)
p_f <- ggplot(f_results, aes(median_difference_days, log10_p, color = cancer_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  geom_point(aes(size = prevalence_percent), alpha = 0.75) +
  ggrepel::geom_text_repel(data = f_labels, aes(label = label), color = "black",
                           size = 2.7, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = cancer_colors) +
  scale_size_continuous(name = "Feature prevalence (%)", range = c(1, 6)) +
  labs(x = "Median TTNT difference (days; feature-positive minus negative)",
       y = expression(-log[10](italic(P))), color = "Cancer type") +
  theme_bw(base_size = 10)
ggsave(file.path(args$outdir, "Figure5F_TTNT_associations.pdf"), p_f,
       width = 11, height = 8, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure5F_TTNT_associations.png"), p_f,
       width = 11, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Shared model helpers
# -----------------------------------------------------------------------------
is_binary01 <- function(x) {
  u <- sort(unique(x[!is.na(x)])); length(u) <= 2L && all(u %in% c(0, 1))
}
align_matrix <- function(x, target_columns) {
  extra <- setdiff(colnames(x), target_columns)
  if (length(extra)) x <- x[, setdiff(colnames(x), extra), drop = FALSE]
  absent <- setdiff(target_columns, colnames(x))
  if (length(absent)) {
    x <- cbind(x, matrix(0, nrow(x), length(absent), dimnames = list(NULL, absent)))
  }
  x[, target_columns, drop = FALSE]
}
make_design <- function(dat, features, cancer_levels = NULL) {
  x <- dat[, features, drop = FALSE]
  for (nm in features) {
    x[[nm]] <- as.numeric(x[[nm]])
    x[[nm]][is.na(x[[nm]])] <- 0
  }
  cancer_model <- ifelse(dat$cancer_type == "NSCLC", "LUCA", dat$cancer_type)
  x$.__cancer__ <- if (is.null(cancer_levels)) factor(cancer_model) else
    factor(cancer_model, levels = cancer_levels)
  model.matrix(~ . - 1, data = x)
}

# -----------------------------------------------------------------------------
# Figure 5G: per-cancer 100-resample elastic-net stability
# -----------------------------------------------------------------------------
set.seed(args$seed)
g_cancers <- km %>% count(cancer_type) %>% filter(n >= 40L) %>% pull(cancer_type)
g_results <- map_dfr(g_cancers, function(ct) {
  dat <- filter(km, cancer_type == ct)
  candidate <- c("female", "late", feature_columns)
  candidate <- candidate[vapply(dat[candidate], function(x) {
    length(unique(x)) > 1L && sum(x == 1, na.rm = TRUE) >= 2L
  }, logical(1))]
  if (length(candidate) < 2L) return(tibble())

  map_dfr(seq_len(args$stability_iterations), function(iteration_number) {
    idx <- sample(seq_len(nrow(dat)), nrow(dat), replace = TRUE)
    boot <- dat[idx, , drop = FALSE]
    valid <- candidate[vapply(boot[candidate], function(x) length(unique(x)) > 1L, logical(1))]
    if (length(valid) < 2L || sum(boot$event_comp) < 5L) return(tibble())
    x <- as.matrix(boot[, valid, drop = FALSE])
    y <- Surv(boot$TTNT_month, boot$event_comp)
    alpha_values <- seq(0.3, 0.9, 0.1)
    fits <- lapply(alpha_values, function(alpha_value) {
      tryCatch(cv.glmnet(x, y, family = "cox", alpha = alpha_value,
                         nfolds = min(5L, nrow(boot) - 1L),
                         standardize = TRUE), error = function(e) NULL)
    })
    names(fits) <- as.character(alpha_values)
    fits <- Filter(Negate(is.null), fits)
    if (!length(fits)) return(tibble())
    score <- vapply(fits, function(z) min(z$cvm), numeric(1))
    best_index <- which.min(score)
    best <- fits[[best_index]]
    beta <- as.matrix(coef(best, s = "lambda.min"))[, 1]
    tibble(feature = names(beta), coefficient = as.numeric(beta),
           iteration = iteration_number, selected_alpha = as.numeric(names(fits)[[best_index]]),
           selected_lambda = best$lambda.min, nfolds = min(5L, nrow(boot) - 1L))
  }) %>% mutate(cancer_type = ct, .before = 1)
})
write_tsv(g_results, "Figure5G_all_resample_coefficients.tsv")
g_summary <- g_results %>% group_by(cancer_type, feature) %>% summarise(
  selection_count = sum(coefficient != 0),
  selection_frequency = selection_count / args$stability_iterations,
  mean_coefficient = ifelse(any(coefficient != 0), mean(coefficient[coefficient != 0]), 0),
  .groups = "drop"
) %>% mutate(
  stability_group = case_when(
    selection_frequency >= 0.5 & abs(mean_coefficient) >= 0.05 ~ "Strong & frequent",
    selection_frequency >= 0.5 ~ "Frequent only",
    abs(mean_coefficient) >= 0.05 ~ "Strong only",
    TRUE ~ "Low"
  ), label = paste(cancer_type, feature, sep = "_")
)
write_tsv(g_summary, "Figure5G_elastic_net_stability.tsv")
g_labels <- g_summary %>% filter(stability_group == "Strong & frequent") %>%
  mutate(rank_score = selection_frequency * abs(mean_coefficient)) %>%
  arrange(desc(rank_score)) %>% slice_head(n = 30)
p_g <- ggplot(g_summary, aes(mean_coefficient, selection_count)) +
  geom_point(aes(fill = stability_group, size = abs(mean_coefficient)),
             shape = 21, color = "grey30", alpha = 0.8) +
  ggrepel::geom_text_repel(data = g_labels, aes(label = label), size = 2.7,
                           fontface = "italic", max.overlaps = Inf) +
  scale_fill_manual(values = c(
    "Strong & frequent"="firebrick", "Frequent only"="darkorange",
    "Strong only"="steelblue", "Low"="grey80"
  )) +
  guides(size = "none") +
  labs(x = "Mean non-zero coefficient", y = "Non-zero coefficient count",
       fill = "Feature stability") + theme_bw(base_size = 10)
ggsave(file.path(args$outdir, "Figure5G_elastic_net_stability.pdf"), p_g,
       width = 10, height = 8, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure5G_elastic_net_stability.png"), p_g,
       width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Leakage-free pan-cancer Cox elastic-net model for Figures 5J/K/L
# -----------------------------------------------------------------------------
set.seed(args$seed)
model_km <- km %>%
  filter(is.finite(TTNT_month), TTNT_month > 0, TTNT_status %in% c(0L, 1L, 2L))
model_genie <- genie %>%
  filter(is.finite(TTNT_month), TTNT_month > 0, TTNT_status %in% c(0L, 1L, 2L))

train_ids <- model_km %>% group_by(cancer_type) %>%
  group_modify(~ tibble(KM_SAMPLE_ID = sample(
    .x$KM_SAMPLE_ID, size = max(1L, floor(0.70 * nrow(.x))), replace = FALSE
  ))) %>% ungroup() %>% pull(KM_SAMPLE_ID)
model_km <- model_km %>% mutate(split = ifelse(KM_SAMPLE_ID %in% train_ids, "train", "test"))
train <- filter(model_km, split == "train")
test <- filter(model_km, split == "test")

model_features <- c("female", "late", feature_columns)
model_features <- model_features[vapply(train[model_features], function(x) {
  length(unique(x)) > 1L && (!is_binary01(x) || sum(x == 1, na.rm = TRUE) >= 2L)
}, logical(1))]

x_train <- make_design(train, model_features)
cancer_levels <- levels(factor(train$cancer_type))
y_train <- Surv(train$TTNT_month, train$event_comp)
alpha_grid <- seq(0.3, 0.9, 0.1)
cv_fits <- lapply(alpha_grid, function(alpha_value) {
  cv.glmnet(x_train, y_train, family = "cox", alpha = alpha_value,
            nfolds = 10L, standardize = TRUE)
})
best_i <- which.min(vapply(cv_fits, function(z) min(z$cvm), numeric(1)))
best_cv <- cv_fits[[best_i]]
best_alpha <- alpha_grid[[best_i]]
best_lambda <- best_cv$lambda.min

score_data <- function(dat) {
  x <- make_design(dat, model_features, cancer_levels)
  x <- align_matrix(x, colnames(x_train))
  as.numeric(predict(best_cv, newx = x, s = best_lambda, type = "link"))
}
train$risk_bad <- score_data(train)
test$risk_bad <- score_data(test)
model_genie$risk_bad <- score_data(model_genie)
risk_cutoff <- median(train$risk_bad, na.rm = TRUE)
assign_group <- function(z) factor(ifelse(z >= risk_cutoff, "High risk", "Low risk"),
                                   levels = c("Low risk", "High risk"))
train$risk_group <- assign_group(train$risk_bad)
test$risk_group <- assign_group(test$risk_bad)
model_genie$risk_group <- assign_group(model_genie$risk_bad)

coef_matrix <- as.matrix(coef(best_cv, s = best_lambda))
model_coefficients <- tibble(
  feature = rownames(coef_matrix), coefficient = as.numeric(coef_matrix[, 1])
) %>% filter(coefficient != 0)
write_tsv(model_coefficients, "integrated_model_coefficients.tsv")
write_tsv(tibble(best_alpha = best_alpha, best_lambda = best_lambda,
                 training_risk_cutoff = risk_cutoff,
                 train_n = nrow(train), test_n = nrow(test)),
          "integrated_model_parameters.tsv")
write_tsv(test, "Figure5J_internal_test_scores.tsv")
write_tsv(model_genie, "Figure5K_GENIE_validation_scores.tsv")

saveRDS(list(
  cv_fit = best_cv, alpha = best_alpha, lambda = best_lambda,
  model_features = model_features, design_columns = colnames(x_train),
  cancer_levels = cancer_levels, training_risk_cutoff = risk_cutoff,
  seed = args$seed
), file.path(args$outdir, "platinum_integrated_cox_model.rds"))

make_km <- function(dat, time, event, title, ylab, risk_table = TRUE) {
  use <- dat %>% filter(is.finite(.data[[time]]), .data[[time]] >= 0,
                        .data[[event]] %in% c(0L, 1L), !is.na(risk_group))
  fit <- survfit(Surv(use[[time]], use[[event]]) ~ risk_group, data = use)
  p <- ggsurvplot(
    fit, data = use, pval = TRUE, risk.table = risk_table, conf.int = FALSE,
    palette = c("#00BFC4", "#F8766D"), xlab = "Months", ylab = ylab,
    legend.title = NULL, title = title, ggtheme = theme_classic(base_size = 10),
    risk.table.height = 0.25
  )
  list(plot = if (risk_table) p$plot / p$table else p$plot, data = use)
}

j <- make_km(test, "TTNT_month", "event_comp", "K-MASTER test cohort",
             "Treatment durability", TRUE)
k <- make_km(model_genie, "TTNT_month", "event_comp", "GENIE-BPC NSCLC",
             "Treatment durability", TRUE)
l_internal <- make_km(test, "OS_month", "OS_event", "K-MASTER test cohort",
                      "Overall survival", FALSE)
l_external <- make_km(model_genie, "OS_month", "OS_event", "GENIE-BPC NSCLC",
                      "Overall survival", FALSE)
write_tsv(j$data, "Figure5J_plot_data.tsv")
write_tsv(k$data, "Figure5K_plot_data.tsv")
write_tsv(l_internal$data, "Figure5L_internal_OS_plot_data.tsv")
write_tsv(l_external$data, "Figure5L_external_OS_plot_data.tsv")
ggsave(file.path(args$outdir, "Figure5J_internal_TTNT.pdf"), j$plot,
       width = 6.5, height = 6.5, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure5K_external_TTNT.pdf"), k$plot,
       width = 6.5, height = 6.5, useDingbats = FALSE)
ggsave(file.path(args$outdir, "Figure5L_overall_survival.pdf"),
       l_internal$plot + l_external$plot, width = 11, height = 4.5, useDingbats = FALSE)

message("Figure 5C/F/G/J/K/L completed: ", normalizePath(args$outdir))
