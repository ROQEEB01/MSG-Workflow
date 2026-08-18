#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript genetic_diversity.R ",
    "af_high_quality.rds depth_high_quality.rds sample_order.tsv output_dir"
  )
}

af_file <- args[1]
depth_file <- args[2]
sample_order_file <- args[3]
output_dir <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

af <- readRDS(af_file)
depth <- readRDS(depth_file)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)

required <- c("pool_index", "pool_size", "chromosome_count", "batch", "site_batch")
if (!all(required %in% names(sample_order))) {
  stop("sample_order.tsv is missing required metadata columns.")
}
if (ncol(af) != nrow(sample_order) || !all(dim(af) == dim(depth))) {
  stop("Allele-frequency/depth matrices do not match sample_order.tsv.")
}

mean_he <- map_dbl(seq_len(nrow(sample_order)), function(i) {
  p <- af[, i]
  D <- depth[, i]
  C <- sample_order$chromosome_count[i]

  he <- rep(NA_real_, length(p))
  valid <- !is.na(p) & !is.na(D) & D > 1 & C > 1
  he[valid] <- (2 * p[valid] * (1 - p[valid])) *
    (C * D[valid]) / ((C - 1) * (D[valid] - 1))

  mean(he, na.rm = TRUE)
})

pool_summary <- sample_order %>%
  mutate(
    mean_He = mean_he,
    mean_depth = colMeans(depth, na.rm = TRUE),
    log_pool_size = log10(pool_size)
  )

he_cor <- cor.test(
  pool_summary$log_pool_size,
  pool_summary$mean_He,
  method = "pearson"
)

results <- tibble(
  n_pools = nrow(pool_summary),
  r = unname(he_cor$estimate),
  p_value = he_cor$p.value
)

write_tsv(pool_summary, file.path(output_dir, "pool_he_summary.tsv"))
write_tsv(results, file.path(output_dir, "pool_size_he_correlation.tsv"))

p_pool <- ggplot(pool_summary, aes(x = reorder(site_batch, mean_He), y = mean_He)) +
  geom_col() +
  geom_text(aes(label = paste0("n=", pool_size)), hjust = -0.1, size = 3) +
  coord_flip() +
  expand_limits(y = max(pool_summary$mean_He, na.rm = TRUE) * 1.15) +
  labs(
    x = NULL,
    y = expression("Mean " * H[e]),
    title = "Expected Heterozygosity by Mosquito Pool"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "expected_heterozygosity_by_pool.png"),
  p_pool, width = 9, height = 7, dpi = 300
)

p_size <- ggplot(pool_summary, aes(x = pool_size, y = mean_He, color = batch)) +
  geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = site_batch), size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(
    x = "Mosquitoes per pool (log scale)",
    y = expression("Mean " * H[e]),
    title = "Pool Size and Expected Heterozygosity",
    subtitle = paste0("r = ", round(he_cor$estimate, 3), ", p = ", signif(he_cor$p.value, 3)),
    color = "Collection period"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "pool_size_vs_he.png"),
  p_size, width = 9, height = 6, dpi = 300
)
