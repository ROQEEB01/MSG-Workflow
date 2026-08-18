#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript filter_snps.R ",
    "snp_quality_features.rds per_pool_AB_table.tsv sample_order.tsv output_dir"
  )
}

features_file <- args[1]
ab_file <- args[2]
sample_order_file <- args[3]
output_dir <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

features <- readRDS(features_file)
ab <- read_tsv(ab_file, show_col_types = FALSE)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)

required <- c(
  "chr", "pos", "mean_depth", "cv_depth", "depth_ratio_B3_B1",
  "in_homopolymer", "pct_near_ends"
)
if (!all(required %in% names(features))) {
  stop("SNP feature table is missing required final-filter metrics.")
}

q_depth_cv <- quantile(features$cv_depth, 0.95, na.rm = TRUE)
q_expression <- quantile(features$depth_ratio_B3_B1, c(0.025, 0.975), na.rm = TRUE)
q_read_ends <- quantile(features$pct_near_ends, 0.95, na.rm = TRUE)

features <- features %>%
  select(-matches("^flag_"), -any_of(c("n_flags", "quality_tier"))) %>%
  mutate(
    flag_low_depth = mean_depth < 30,
    flag_depth_cv = cv_depth > q_depth_cv,
    flag_expression = depth_ratio_B3_B1 < q_expression[[1]] |
      depth_ratio_B3_B1 > q_expression[[2]],
    flag_homopolymer = as.character(in_homopolymer) %in% c("TRUE", "1"),
    flag_read_ends = pct_near_ends > q_read_ends
  )

flag_cols <- grep("^flag_", names(features), value = TRUE)
features <- features %>%
  mutate(
    n_flags = rowSums(across(all_of(flag_cols)), na.rm = TRUE),
    quality_tier = case_when(
      n_flags == 0 ~ "High quality",
      n_flags == 1 ~ "Minor concern",
      TRUE ~ "Potential artifact"
    ),
    key = paste(chr, pos, sep = "_")
  )

thresholds <- tibble(
  metric = c(
    "mean_depth_min", "depth_cv_95th", "depth_ratio_2.5th",
    "depth_ratio_97.5th", "read_end_95th"
  ),
  value = c(30, q_depth_cv, q_expression[[1]], q_expression[[2]], q_read_ends)
)

hq <- features %>% filter(quality_tier == "High quality")
if (nrow(hq) == 0) stop("No SNPs passed the final quality filter.")

ab <- ab %>% mutate(key = paste(chr, pos, sep = "_"))
match_idx <- match(hq$key, ab$key)
if (anyNA(match_idx)) {
  stop("One or more high-quality SNPs were not found in the per-pool table.")
}
ab_hq <- ab[match_idx, , drop = FALSE]

n_pools <- nrow(sample_order)
af_cols <- paste0("AB_", sample_order$pool_index)
depth_cols <- paste0("depth_", sample_order$pool_index)
if (!all(c(af_cols, depth_cols) %in% names(ab_hq))) {
  stop("Per-pool table does not match sample_order.tsv pool indices.")
}

af_hq <- as.matrix(ab_hq[, af_cols])
depth_hq <- as.matrix(ab_hq[, depth_cols])
storage.mode(af_hq) <- "numeric"
storage.mode(depth_hq) <- "numeric"
rownames(af_hq) <- hq$key
rownames(depth_hq) <- hq$key
colnames(af_hq) <- sample_order$site_batch
colnames(depth_hq) <- sample_order$site_batch

snp_metadata <- ab_hq %>%
  select(-all_of(c(af_cols, depth_cols)))

write_tsv(features, file.path(output_dir, "snp_quality_classification.tsv"))
write_tsv(hq, file.path(output_dir, "snps_high_quality.tsv"))
write_tsv(thresholds, file.path(output_dir, "snp_filter_thresholds.tsv"))
write_tsv(ab_hq, file.path(output_dir, "per_pool_AB_high_quality.tsv"))
write_tsv(snp_metadata, file.path(output_dir, "snp_metadata_high_quality.tsv"))
saveRDS(features, file.path(output_dir, "snp_quality_classification.rds"))
saveRDS(af_hq, file.path(output_dir, "af_high_quality.rds"))
saveRDS(depth_hq, file.path(output_dir, "depth_high_quality.rds"))
