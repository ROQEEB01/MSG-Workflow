#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript pcoa.R fst_high_quality.rds sample_order.tsv output_dir")
}

fst_file <- args[1]
sample_order_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fst_matrix <- readRDS(fst_file)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)

if (!all(c("site_batch", "site", "batch") %in% names(sample_order))) {
  stop("sample_order.tsv must contain site_batch, site, and batch.")
}

pool_names <- sample_order$site_batch
if (!all(pool_names %in% rownames(fst_matrix))) {
  stop("FST matrix and sample_order.tsv contain different pool labels.")
}
fst_matrix <- fst_matrix[pool_names, pool_names, drop = FALSE]

pcoa <- cmdscale(as.dist(fst_matrix), k = 2, eig = TRUE)
positive_eig <- pcoa$eig[pcoa$eig > 0]
if (length(positive_eig) < 2) {
  stop("PCoA requires at least two positive eigenvalues.")
}

variance <- positive_eig / sum(positive_eig) * 100
coords <- as.data.frame(pcoa$points[, 1:2, drop = FALSE]) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  bind_cols(sample_order %>% select(pool_index, site_batch, site, batch))

write_tsv(coords, file.path(output_dir, "pcoa_coordinates.tsv"))
write_tsv(
  tibble(axis = seq_along(positive_eig), eigenvalue = positive_eig, percent = variance),
  file.path(output_dir, "pcoa_positive_eigenvalues.tsv")
)

p <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = batch, label = site)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    x = sprintf("PCoA1 (%.1f%%)", variance[1]),
    y = sprintf("PCoA2 (%.1f%%)", variance[2]),
    title = "PCoA of Pairwise FST",
    color = "Collection period"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "pcoa_pairwise_fst.png"),
  p, width = 8, height = 6, dpi = 300
)
