#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: Rscript build_snp_quality_features.R ",
    "per_pool_AB_table.tsv sample_order.tsv snp_context.tsv ",
    "read_features_dir output_dir"
  )
}

ab_file <- args[1]
sample_order_file <- args[2]
context_file <- args[3]
read_features_dir <- args[4]
output_dir <- args[5]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ab <- read_tsv(ab_file, show_col_types = FALSE)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)
context <- read_tsv(context_file, show_col_types = FALSE)

required_sample_cols <- c(
  "pool_index", "sample_id", "pool_size", "ploidy",
  "chromosome_count", "batch", "site", "site_batch"
)
if (!all(required_sample_cols %in% names(sample_order))) {
  stop("sample_order.tsv is missing required columns.")
}

n_pools <- nrow(sample_order)
af_cols <- paste0("AB_", sample_order$pool_index)
depth_cols <- paste0("depth_", sample_order$pool_index)
if (!all(c("chr", "pos", af_cols, depth_cols) %in% names(ab))) {
  stop("Per-pool table does not match sample_order.tsv pool indices.")
}

if (!all(c("chr", "pos", "in_homopolymer") %in% names(context))) {
  stop("snp_context.tsv must contain chr, pos, and in_homopolymer.")
}

read_feature_files <- list.files(
  read_features_dir,
  pattern = "_read_features\\.tsv$",
  full.names = TRUE
)
if (length(read_feature_files) == 0) {
  stop("No *_read_features.tsv files found in read_features_dir.")
}

read_required <- c(
  "chr", "pos", "rel_read_pos", "base_qual", "map_qual",
  "read_length", "is_reverse", "n_mismatches"
)

read_summaries <- lapply(read_feature_files, function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  if (!all(read_required %in% names(x))) {
    stop("Missing required read-feature columns in: ", path)
  }

  x %>%
    group_by(chr, pos) %>%
    summarise(
      sum_rel_read_pos = sum(rel_read_pos, na.rm = TRUE),
      sum_sq_rel_read_pos = sum(rel_read_pos^2, na.rm = TRUE),
      n_near_ends = sum(rel_read_pos < 0.1 | rel_read_pos > 0.9, na.rm = TRUE),
      sum_base_qual = sum(base_qual, na.rm = TRUE),
      n_low_bq = sum(base_qual < 10, na.rm = TRUE),
      sum_map_qual = sum(map_qual, na.rm = TRUE),
      n_low_mq = sum(map_qual < 20, na.rm = TRUE),
      sum_read_length = sum(read_length, na.rm = TRUE),
      n_reverse = sum(is_reverse, na.rm = TRUE),
      sum_mismatches = sum(n_mismatches, na.rm = TRUE),
      n_reads = n(),
      .groups = "drop"
    )
})

read_stats <- bind_rows(read_summaries) %>%
  group_by(chr, pos) %>%
  summarise(
    sum_rel_read_pos = sum(sum_rel_read_pos),
    sum_sq_rel_read_pos = sum(sum_sq_rel_read_pos),
    n_near_ends = sum(n_near_ends),
    sum_base_qual = sum(sum_base_qual),
    n_low_bq = sum(n_low_bq),
    sum_map_qual = sum(sum_map_qual),
    n_low_mq = sum(n_low_mq),
    sum_read_length = sum(sum_read_length),
    n_reverse = sum(n_reverse),
    sum_mismatches = sum(sum_mismatches),
    n_reads = sum(n_reads),
    .groups = "drop"
  ) %>%
  mutate(
    mean_rel_read_pos = sum_rel_read_pos / n_reads,
    sd_rel_read_pos = sqrt(pmax(
      sum_sq_rel_read_pos / n_reads - mean_rel_read_pos^2,
      0
    )),
    pct_near_ends = n_near_ends / n_reads,
    mean_base_qual = sum_base_qual / n_reads,
    pct_low_bq = n_low_bq / n_reads,
    mean_map_qual = sum_map_qual / n_reads,
    pct_low_mq = n_low_mq / n_reads,
    mean_read_length = sum_read_length / n_reads,
    pct_reverse = n_reverse / n_reads,
    strand_bias = abs(pct_reverse - 0.5) * 2,
    mean_mismatches = sum_mismatches / n_reads
  ) %>%
  select(
    chr, pos, mean_rel_read_pos, sd_rel_read_pos, pct_near_ends,
    mean_base_qual, pct_low_bq, mean_map_qual, pct_low_mq,
    mean_read_length, pct_reverse, strand_bias, mean_mismatches, n_reads
  )

af_matrix <- as.matrix(ab[, af_cols])
depth_matrix <- as.matrix(ab[, depth_cols])
storage.mode(af_matrix) <- "numeric"
storage.mode(depth_matrix) <- "numeric"

collection_order <- readr::parse_number(sample_order$batch)
if (anyNA(collection_order)) {
  stop("Could not derive collection order from the batch column.")
}

weights <- sample_order$chromosome_count
collection_beta <- rep(NA_real_, nrow(af_matrix))
collection_p <- rep(NA_real_, nrow(af_matrix))

for (i in seq_len(nrow(af_matrix))) {
  model_data <- tibble(
    freq = af_matrix[i, ],
    collection_order = collection_order,
    weight = weights
  ) %>%
    filter(complete.cases(.))

  if (nrow(model_data) < 3 || length(unique(model_data$collection_order)) < 2) {
    next
  }

  fit <- try(
    lm(freq ~ collection_order, weights = weight, data = model_data),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) next

  coefs <- summary(fit)$coefficients
  if ("collection_order" %in% rownames(coefs)) {
    collection_beta[i] <- coefs["collection_order", "Estimate"]
    collection_p[i] <- coefs["collection_order", "Pr(>|t|)"]
  }
}

b1 <- which(sample_order$batch == "B1")
b3 <- which(sample_order$batch == "B3")
if (length(b1) == 0 || length(b3) == 0) {
  stop("B1 and B3 pools are required for the B3:B1 depth-ratio metric.")
}

mean_af <- rowMeans(af_matrix, na.rm = TRUE)
allele_states_col <- if ("allele_states" %in% names(ab)) {
  "allele_states"
} else if ("minor_allele" %in% names(ab)) {
  "minor_allele"
} else {
  NA_character_
}

snp_features <- tibble(
  chr = ab$chr,
  pos = ab$pos,
  ref = if ("ref" %in% names(ab)) ab$ref else NA_character_,
  allele_states = if (!is.na(allele_states_col)) ab[[allele_states_col]] else NA_character_,
  collection_beta = collection_beta,
  collection_p = collection_p,
  collection_padj = p.adjust(collection_p, method = "BH"),
  abs_collection_beta = abs(collection_beta),
  global_maf = pmin(mean_af, 1 - mean_af),
  mean_depth = rowMeans(depth_matrix, na.rm = TRUE),
  median_depth = apply(depth_matrix, 1, median, na.rm = TRUE),
  sd_depth = apply(depth_matrix, 1, sd, na.rm = TRUE),
  min_depth = apply(depth_matrix, 1, min, na.rm = TRUE),
  max_depth = apply(depth_matrix, 1, max, na.rm = TRUE),
  mean_depth_B1 = rowMeans(depth_matrix[, b1, drop = FALSE], na.rm = TRUE),
  mean_depth_B3 = rowMeans(depth_matrix[, b3, drop = FALSE], na.rm = TRUE)
) %>%
  mutate(
    cv_depth = sd_depth / mean_depth,
    depth_ratio_B3_B1 = mean_depth_B3 / pmax(mean_depth_B1, 1),
    key = paste(chr, pos, sep = "_")
  ) %>%
  left_join(context %>% select(chr, pos, in_homopolymer), by = c("chr", "pos")) %>%
  left_join(read_stats, by = c("chr", "pos"))

if (!all(is.na(snp_features$ref)) && !all(is.na(snp_features$allele_states))) {
  alt <- map2_chr(snp_features$ref, snp_features$allele_states, function(ref, states) {
    if (is.na(ref) || is.na(states)) return(NA_character_)
    alleles <- str_split(states, "/", simplify = FALSE)[[1]]
    other <- alleles[alleles != ref]
    if (length(other) == 0) NA_character_ else other[1]
  })

  snp_features <- snp_features %>%
    mutate(
      alt = alt,
      mutation = if_else(!is.na(ref) & !is.na(alt), paste0(ref, ">", alt), NA_character_),
      is_transition = mutation %in% c("A>G", "G>A", "C>T", "T>C"),
      is_AG_or_TC = mutation %in% c("A>G", "T>C"),
      is_CT_or_GA = mutation %in% c("C>T", "G>A")
    )
}

correlation_features <- intersect(
  c(
    "mean_depth", "cv_depth", "depth_ratio_B3_B1", "global_maf",
    "pct_near_ends", "mean_base_qual", "mean_map_qual",
    "strand_bias", "mean_mismatches"
  ),
  names(snp_features)
)

feature_correlations <- map_dfr(correlation_features, function(feature) {
  x <- snp_features[[feature]]
  y <- snp_features$abs_collection_beta
  keep <- complete.cases(x, y)

  if (sum(keep) < 3 || sd(x[keep]) == 0) {
    return(tibble(feature = feature, rho = NA_real_, p_value = NA_real_))
  }

  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  tibble(
    feature = feature,
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(desc(abs(rho)))

write_tsv(snp_features, file.path(output_dir, "snp_quality_features.tsv"))
write_tsv(feature_correlations, file.path(output_dir, "feature_correlations.tsv"))
saveRDS(snp_features, file.path(output_dir, "snp_quality_features.rds"))
saveRDS(af_matrix, file.path(output_dir, "af_matrix_all_biallelic.rds"))
saveRDS(depth_matrix, file.path(output_dir, "depth_matrix_all_biallelic.rds"))
