#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript fst_analysis.R ",
    "snps_high_quality.tsv tier1_sync.fst sample_order.tsv output_dir"
  )
}

hq_file <- args[1]
fst_file <- args[2]
sample_order_file <- args[3]
output_dir <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hq <- read_tsv(hq_file, show_col_types = FALSE)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)
fst_raw <- read.table(
  fst_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

if (!all(c("chr", "pos") %in% names(hq))) {
  stop("High-quality SNP table must contain chr and pos.")
}
if (!all(c("pool_index", "site_batch") %in% names(sample_order))) {
  stop("sample_order.tsv must contain pool_index and site_batch.")
}

hq_key <- paste(hq$chr, hq$pos, sep = "_")
fst_key <- paste(fst_raw[[1]], fst_raw[[2]], sep = "_")
idx <- match(hq_key, fst_key)
idx <- idx[!is.na(idx)]
if (length(idx) == 0) stop("No high-quality SNPs matched the FST file.")

fst_hq <- fst_raw[idx, , drop = FALSE]
if (ncol(fst_hq) < 6) stop("Unexpected PoPoolation2 FST format.")

pair_columns <- 6:ncol(fst_hq)
pair_ids <- map_chr(fst_hq[1, pair_columns, drop = FALSE], ~ sub("=.*", "", as.character(.x)))

pairwise <- map2_dfr(pair_columns, pair_ids, function(col_idx, pair_id) {
  values <- suppressWarnings(as.numeric(sub(".*=", "", fst_hq[[col_idx]])))
  parts <- str_split_fixed(pair_id, ":", 2)

  tibble(
    pair_id = pair_id,
    pool1_index = as.integer(parts[1]),
    pool2_index = as.integer(parts[2]),
    average_FST = mean(values, na.rm = TRUE),
    n_sites = sum(!is.na(values))
  )
}) %>%
  mutate(
    pool1 = sample_order$site_batch[pool1_index],
    pool2 = sample_order$site_batch[pool2_index]
  )

n_pools <- nrow(sample_order)
fst_matrix <- matrix(0, nrow = n_pools, ncol = n_pools)
rownames(fst_matrix) <- sample_order$site_batch
colnames(fst_matrix) <- sample_order$site_batch

for (i in seq_len(nrow(pairwise))) {
  a <- pairwise$pool1_index[i]
  b <- pairwise$pool2_index[i]
  fst_matrix[a, b] <- pairwise$average_FST[i]
  fst_matrix[b, a] <- pairwise$average_FST[i]
}

write_tsv(pairwise, file.path(output_dir, "pairwise_fst_high_quality.tsv"))
write_csv(
  as.data.frame(fst_matrix) %>% rownames_to_column("pool"),
  file.path(output_dir, "pairwise_fst_matrix_high_quality.csv")
)
saveRDS(fst_matrix, file.path(output_dir, "fst_high_quality.rds"))

fst_long <- as.data.frame(as.table(fst_matrix)) %>%
  as_tibble() %>%
  rename(pool1 = Var1, pool2 = Var2, FST = Freq) %>%
  filter(pool1 != pool2)

p_heatmap <- ggplot(fst_long, aes(x = pool1, y = pool2, fill = FST)) +
  geom_tile() +
  labs(
    x = NULL,
    y = NULL,
    title = "Pairwise FST Among Mosquito Pools",
    fill = "FST"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(
  file.path(output_dir, "pairwise_fst_heatmap.png"),
  p_heatmap, width = 10, height = 9, dpi = 300
)

if (anyNA(fst_matrix)) {
  warning("Dendrogram skipped because the FST matrix contains missing values.")
} else {
  hc <- hclust(as.dist(fst_matrix), method = "average")
  pdf(file.path(output_dir, "fst_dendrogram.pdf"), width = 11, height = 7)
  plot(hc, main = "Hierarchical Clustering of Pairwise FST", xlab = "", sub = "")
  dev.off()
}
