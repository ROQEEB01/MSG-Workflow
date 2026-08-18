#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(
    "Usage: Rscript allelic_richness.R ",
    "af_high_quality.rds sample_order.tsv output_dir"
  )
}

af_file <- args[1]
sample_order_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

af <- readRDS(af_file)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)

required <- c("pool_index", "pool_size", "chromosome_count", "batch", "site_batch")
if (!all(required %in% names(sample_order))) {
  stop("sample_order.tsv is missing required metadata columns.")
}
if (ncol(af) != nrow(sample_order)) {
  stop("Allele-frequency matrix does not match sample_order.tsv.")
}

rarefy_ar_locus <- function(p_hat, chromosome_count, rarefy_to) {
  if (is.na(p_hat) || chromosome_count < rarefy_to) return(NA_real_)
  if (p_hat <= 0 || p_hat >= 1) return(1)

  minor_count <- round(p_hat * chromosome_count)
  major_count <- chromosome_count - minor_count
  if (minor_count == 0 || major_count == 0) return(1)

  p_miss_minor <- dhyper(0, minor_count, major_count, rarefy_to)
  p_miss_major <- dhyper(0, major_count, minor_count, rarefy_to)
  2 - p_miss_minor - p_miss_major
}

n_snps <- nrow(af)
n_pools <- ncol(af)
rarefy_to <- min(sample_order$chromosome_count)

raw_polymorphic <- map_int(seq_len(n_pools), function(i) {
  x <- af[, i]
  sum(x > 0 & x < 1, na.rm = TRUE)
})

raw_mean_ar <- map_dbl(seq_len(n_pools), function(i) {
  x <- af[, i]
  mean(ifelse(x > 0 & x < 1, 2, 1), na.rm = TRUE)
})

rarefied_ar <- map_dbl(seq_len(n_pools), function(i) {
  mean(
    map_dbl(
      af[, i],
      rarefy_ar_locus,
      chromosome_count = sample_order$chromosome_count[i],
      rarefy_to = rarefy_to
    ),
    na.rm = TRUE
  )
})

ar_summary <- sample_order %>%
  mutate(
    log_pool_size = log10(pool_size),
    raw_polymorphic = raw_polymorphic,
    raw_pct_polymorphic = 100 * raw_polymorphic / n_snps,
    raw_mean_ar = raw_mean_ar,
    rarefied_mean_ar = rarefied_ar,
    rarefied_to_chromosomes = rarefy_to
  )

ar_cor <- cor.test(
  ar_summary$log_pool_size,
  ar_summary$rarefied_mean_ar,
  method = "pearson"
)

anova_fit <- aov(rarefied_mean_ar ~ factor(batch), data = ar_summary)
anova_table <- as.data.frame(summary(anova_fit)[[1]]) %>%
  rownames_to_column("term")

write_tsv(ar_summary, file.path(output_dir, "allelic_richness_summary.tsv"))
write_tsv(
  tibble(
    n_pools = n_pools,
    rarefied_to_chromosomes = rarefy_to,
    r = unname(ar_cor$estimate),
    p_value = ar_cor$p.value
  ),
  file.path(output_dir, "pool_size_rarefied_ar_correlation.tsv")
)
write_tsv(anova_table, file.path(output_dir, "rarefied_ar_collection_period_anova.tsv"))

p_size <- ggplot(ar_summary, aes(x = pool_size, y = rarefied_mean_ar, color = batch)) +
  geom_point(size = 3) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = site_batch), size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(
    x = "Mosquitoes per pool (log scale)",
    y = paste0("Mean allelic richness (rarefied to ", rarefy_to, " chromosomes)"),
    title = "Pool Size and Rarefied Allelic Richness",
    subtitle = paste0("r = ", round(ar_cor$estimate, 3), ", p = ", signif(ar_cor$p.value, 3)),
    color = "Collection period"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "pool_size_vs_rarefied_allelic_richness.png"),
  p_size, width = 9, height = 6, dpi = 300
)

p_batch <- ggplot(ar_summary, aes(x = factor(batch), y = rarefied_mean_ar)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = batch), width = 0.12, size = 2.5, show.legend = FALSE) +
  geom_text_repel(aes(label = site_batch), size = 2.5, max.overlaps = 20) +
  labs(
    x = "Collection period",
    y = paste0("Mean allelic richness (rarefied to ", rarefy_to, " chromosomes)"),
    title = "Rarefied Allelic Richness by Collection Period"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "rarefied_allelic_richness_by_collection_period.png"),
  p_batch, width = 8, height = 6, dpi = 300
)

curve_start <- max(2, rarefy_to - 4)
rarefy_levels <- seq(curve_start, rarefy_to, by = 2)
if (tail(rarefy_levels, 1) != rarefy_to) rarefy_levels <- c(rarefy_levels, rarefy_to)

curve_data <- map_dfr(seq_len(n_pools), function(i) {
  map_dfr(rarefy_levels, function(k) {
    if (k > sample_order$chromosome_count[i]) return(tibble())
    values <- map_dbl(
      af[, i],
      rarefy_ar_locus,
      chromosome_count = sample_order$chromosome_count[i],
      rarefy_to = k
    )
    tibble(
      site_batch = sample_order$site_batch[i],
      batch = sample_order$batch[i],
      chromosomes_sampled = k,
      mean_ar = mean(values, na.rm = TRUE)
    )
  })
})

curve_summary <- curve_data %>%
  group_by(batch, chromosomes_sampled) %>%
  summarise(
    mean_ar = mean(mean_ar, na.rm = TRUE),
    se_ar = sd(mean_ar, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_tsv(curve_data, file.path(output_dir, "rarefaction_curves_by_pool.tsv"))

p_curve <- ggplot(curve_summary, aes(x = chromosomes_sampled, y = mean_ar, color = batch)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_ribbon(
    aes(ymin = mean_ar - se_ar, ymax = mean_ar + se_ar, fill = batch),
    alpha = 0.15, color = NA, show.legend = FALSE
  ) +
  labs(
    x = "Chromosome copies sampled",
    y = "Mean allelic richness",
    title = "Allelic Richness Rarefaction Curves",
    color = "Collection period"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "allelic_richness_rarefaction_curves.png"),
  p_curve, width = 8, height = 6, dpi = 300
)
