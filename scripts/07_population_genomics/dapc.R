#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(adegenet))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript dapc.R af_high_quality.rds sample_order.tsv output_dir")
}

af_file <- args[1]
sample_order_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

af <- readRDS(af_file)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)

if (!all(c("site_batch", "site", "batch") %in% names(sample_order))) {
  stop("sample_order.tsv must contain site_batch, site, and batch.")
}
if (ncol(af) != nrow(sample_order)) {
  stop("Allele-frequency matrix does not match sample_order.tsv.")
}

x <- t(af)
rownames(x) <- sample_order$site_batch

all_na <- apply(x, 2, function(v) all(is.na(v)))
x <- x[, !all_na, drop = FALSE]
if (ncol(x) < 2) stop("Too few SNPs available for DAPC.")

for (j in seq_len(ncol(x))) {
  if (anyNA(x[, j])) {
    x[is.na(x[, j]), j] <- mean(x[, j], na.rm = TRUE)
  }
}

groups <- factor(sample_order$batch)
group_sizes <- table(groups)
if (length(group_sizes) < 2) stop("DAPC requires at least two collection groups.")

max_pcs <- min(nrow(x) - nlevels(groups), ncol(x), 10)
if (max_pcs < 2) stop("Too few degrees of freedom for DAPC cross-validation.")
training_prop <- max(0.7, 1 - (2 / min(group_sizes)))

set.seed(123)
xval <- xvalDapc(
  x,
  groups,
  n.pca.max = max_pcs,
  training.set = training_prop,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.pca = NULL,
  n.rep = 30,
  xval.plot = FALSE
)

lowest_name <- grep("Number of PCs Achieving Lowest MSE", names(xval), value = TRUE)
mse_name <- grep("Root Mean Squared Error", names(xval), value = TRUE)
success_name <- grep("Mean Successful Assignment", names(xval), value = TRUE)

if (length(lowest_name) > 0 && length(xval[[lowest_name[1]]]) > 0) {
  optimal_pcs <- as.integer(round(xval[[lowest_name[1]]][1]))
} else if (length(mse_name) > 0) {
  mse <- as.numeric(xval[[mse_name[1]]])
  pcs_tested <- seq_along(mse) + 1
  optimal_pcs <- pcs_tested[which.min(mse)]
} else {
  stop("Could not determine the number of PCs from xvalDapc output.")
}

optimal_pcs <- max(2, min(optimal_pcs, max_pcs))
n_da <- min(nlevels(groups) - 1, 2)

dapc_fit <- dapc(
  x,
  grp = groups,
  n.pca = optimal_pcs,
  n.da = n_da
)

if (length(mse_name) > 0 || length(success_name) > 0) {
  mse <- if (length(mse_name) > 0) as.numeric(xval[[mse_name[1]]]) else NULL
  success <- if (length(success_name) > 0) as.numeric(xval[[success_name[1]]]) else NULL
  n_vals <- max(length(mse), length(success))
  xval_table <- tibble(
    n_pca = seq_len(n_vals) + 1,
    mean_success = if (is.null(success)) rep(NA_real_, n_vals) else c(success, rep(NA_real_, n_vals - length(success))),
    rmse = if (is.null(mse)) rep(NA_real_, n_vals) else c(mse, rep(NA_real_, n_vals - length(mse)))
  )
  write_tsv(xval_table, file.path(output_dir, "dapc_cross_validation.tsv"))
}

variance <- dapc_fit$eig / sum(dapc_fit$eig) * 100
coords <- tibble(
  site_batch = sample_order$site_batch,
  site = sample_order$site,
  batch = sample_order$batch,
  LD1 = dapc_fit$ind.coord[, 1],
  LD2 = if (ncol(dapc_fit$ind.coord) > 1) dapc_fit$ind.coord[, 2] else 0
)
write_tsv(coords, file.path(output_dir, "dapc_coordinates.tsv"))
saveRDS(dapc_fit, file.path(output_dir, "dapc_result.rds"))

p_dapc <- ggplot(coords, aes(x = LD1, y = LD2, color = batch, label = site)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    x = sprintf("LD1 (%.1f%%)", variance[1]),
    y = sprintf("LD2 (%.1f%%)", ifelse(length(variance) > 1, variance[2], 0)),
    title = "DAPC of High-Quality SNPs",
    subtitle = paste0(optimal_pcs, " PCs retained"),
    color = "Collection period"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "dapc_collection_period.png"),
  p_dapc, width = 8, height = 6, dpi = 300
)

membership <- as.data.frame(dapc_fit$posterior) %>%
  rownames_to_column("site_batch") %>%
  as_tibble() %>%
  left_join(sample_order %>% select(site_batch, batch), by = "site_batch")

membership_long <- membership %>%
  pivot_longer(
    cols = all_of(levels(groups)),
    names_to = "predicted_group",
    values_to = "probability"
  )
write_tsv(membership_long, file.path(output_dir, "dapc_membership_probabilities.tsv"))

p_membership <- ggplot(
  membership_long,
  aes(x = site_batch, y = probability, fill = predicted_group)
) +
  geom_col() +
  facet_wrap(~ batch, scales = "free_x") +
  labs(
    x = "Mosquito pool",
    y = "Membership probability",
    title = "DAPC Membership Probabilities",
    fill = "Predicted group"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(
  file.path(output_dir, "dapc_membership_probabilities.png"),
  p_membership, width = 11, height = 6, dpi = 300
)
