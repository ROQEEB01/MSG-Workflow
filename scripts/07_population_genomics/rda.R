#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vegan))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript rda.R ",
    "af_high_quality.rds sample_order.tsv site_coordinates.csv output_dir"
  )
}

af_file <- args[1]
sample_order_file <- args[2]
coords_file <- args[3]
output_dir <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load high-quality allele frequencies and metadata
# ------------------------------------------------------------

af <- readRDS(af_file)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)
coords <- read_csv(coords_file, show_col_types = FALSE)

required_sample <- c("pool_index", "site_batch", "site", "batch")
required_coords <- c("SiteID", "Lat", "Long")

if (!all(required_sample %in% names(sample_order))) {
  stop("sample_order.tsv is missing required metadata columns.")
}
if (!all(required_coords %in% names(coords))) {
  stop("site_coordinates.csv must contain SiteID, Lat, and Long.")
}
if (ncol(af) != nrow(sample_order)) {
  stop("Allele-frequency matrix does not match sample_order.tsv.")
}
if (is.null(rownames(af))) {
  stop("Allele-frequency matrix must have SNP IDs as row names.")
}

# Keep the same pool order throughout the analysis.
if (!is.null(colnames(af))) {
  if (!all(sample_order$site_batch %in% colnames(af))) {
    stop("Allele-frequency matrix and sample_order.tsv contain different pool labels.")
  }
  af <- af[, sample_order$site_batch, drop = FALSE]
}

# RDA requires pools as rows and SNPs as columns.
Y <- t(af)
if (anyNA(Y)) {
  stop("RDA input contains missing allele frequencies.")
}

# Reconstruct the Hellinger-transformed SNP matrix used in the original analysis.
Y_hel <- decostand(Y, method = "hellinger")
GEA_SNPs <- as.data.frame(Y_hel)

meta <- sample_order %>%
  select(pool_index, site_batch, site, batch) %>%
  left_join(coords, by = c("site" = "SiteID")) %>%
  mutate(time_num = readr::parse_number(batch))

if (anyNA(meta$Lat) || anyNA(meta$Long) || anyNA(meta$time_num)) {
  stop("Missing time, latitude, or longitude values for one or more pools.")
}

rownames(GEA_SNPs) <- meta$site_batch

# ------------------------------------------------------------
# Standardize predictors
# ------------------------------------------------------------

env_scaled <- scale(
  meta %>% select(time_num, Lat, Long),
  center = TRUE,
  scale = TRUE
)
Variables <- as.data.frame(env_scaled)
rownames(Variables) <- meta$site_batch

# ------------------------------------------------------------
# PCA-derived background structure
# ------------------------------------------------------------

pca <- rda(Y_hel, scale = TRUE)
pca_importance <- summary(pca)$cont$importance

write_tsv(
  as.data.frame(pca_importance) %>%
    rownames_to_column("metric") %>%
    as_tibble(),
  file.path(output_dir, "pca_importance.tsv")
)

PC1 <- scores(pca, choices = 1, display = "sites", scaling = 0)[, 1]
Variables$PC1 <- PC1[rownames(Variables)]

# ------------------------------------------------------------
# Forward selection
# ------------------------------------------------------------

RDA0 <- rda(GEA_SNPs ~ 1, data = Variables)
RDAfull <- rda(GEA_SNPs ~ PC1 + time_num + Lat + Long, data = Variables)

set.seed(123)
forward_model <- ordiR2step(
  RDA0,
  scope = formula(RDAfull),
  Pin = 0.01,
  R2permutations = 1000,
  R2scope = TRUE,
  trace = FALSE
)

forward_table <- as.data.frame(forward_model$anova) %>%
  rownames_to_column("step")
write_tsv(forward_table, file.path(output_dir, "rda_forward_selection.tsv"))

# ------------------------------------------------------------
# Variance partitioning / partial RDA
# ------------------------------------------------------------

RDAfull <- rda(GEA_SNPs ~ PC1 + time_num + Lat + Long, data = Variables)
pRDA_time <- rda(
  GEA_SNPs ~ time_num + Condition(Lat + Long + PC1),
  data = Variables
)
pRDA_geo <- rda(
  GEA_SNPs ~ Lat + Long + Condition(time_num + PC1),
  data = Variables
)
pRDA_struct <- rda(
  GEA_SNPs ~ PC1 + Condition(time_num + Lat + Long),
  data = Variables
)

model_summary <- function(model, label, conditioned_on) {
  set.seed(123)
  test <- anova(model, permutations = 999)
  tibble(
    model = label,
    conditioned_on = conditioned_on,
    adjusted_R2 = unname(RsquareAdj(model)$adj.r.squared),
    p_value = test$`Pr(>F)`[1],
    permutations = 999
  )
}

variance_results <- bind_rows(
  model_summary(RDAfull, "Full", "None"),
  model_summary(pRDA_time, "Collection period", "Latitude + Longitude + PC1"),
  model_summary(pRDA_geo, "Geography", "Collection period + PC1"),
  model_summary(pRDA_struct, "PC1", "Collection period + Latitude + Longitude")
)

write_tsv(
  variance_results,
  file.path(output_dir, "rda_variance_partitioning.tsv")
)

# ------------------------------------------------------------
# Extreme RDA1 loadings
# ------------------------------------------------------------

extract_extreme_loadings <- function(model, model_name, z = 2.5) {
  loadings <- scores(model, choices = 1, display = "species")
  values <- as.numeric(loadings[, 1])
  snp_ids <- rownames(loadings)

  lower <- mean(values) - z * sd(values)
  upper <- mean(values) + z * sd(values)

  tibble(
    snp = snp_ids,
    loading = values,
    model = model_name,
    lower_threshold = lower,
    upper_threshold = upper,
    extreme = values < lower | values > upper
  )
}

# Standard time-only RDA.
RDA_time <- rda(GEA_SNPs ~ time_num, data = Variables)
standard_loadings <- extract_extreme_loadings(
  RDA_time,
  "Time only"
)

# Time effect after conditioning on PC1.
RDA_time_pc1 <- rda(
  GEA_SNPs ~ time_num + Condition(PC1),
  data = Variables
)
corrected_loadings <- extract_extreme_loadings(
  RDA_time_pc1,
  "Time | PC1"
)

standard_extreme <- standard_loadings %>% filter(extreme)
corrected_extreme <- corrected_loadings %>% filter(extreme)

consensus_ids <- intersect(
  standard_extreme$snp,
  corrected_extreme$snp
)

consensus <- corrected_extreme %>%
  filter(snp %in% consensus_ids) %>%
  arrange(desc(abs(loading)))

write_tsv(
  standard_loadings,
  file.path(output_dir, "rda_time_loadings.tsv")
)
write_tsv(
  corrected_loadings,
  file.path(output_dir, "rda_time_pc1_loadings.tsv")
)
write_tsv(
  standard_extreme,
  file.path(output_dir, "rda_time_extreme_loadings.tsv")
)
write_tsv(
  corrected_extreme,
  file.path(output_dir, "rda_time_pc1_extreme_loadings.tsv")
)
write_tsv(
  consensus,
  file.path(output_dir, "rda_consensus_extreme_loadings.tsv")
)

# ------------------------------------------------------------
# Loading distribution
# ------------------------------------------------------------

plot_data <- corrected_loadings %>%
  mutate(consensus = snp %in% consensus_ids)

thresholds <- distinct(
  plot_data,
  lower_threshold,
  upper_threshold
)

p <- ggplot(
  plot_data,
  aes(x = seq_along(loading), y = loading)
) +
  geom_point(aes(shape = consensus), alpha = 0.7) +
  geom_hline(
    yintercept = thresholds$lower_threshold[1],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = thresholds$upper_threshold[1],
    linetype = "dashed"
  ) +
  labs(
    x = "SNP index",
    y = "RDA1 loading",
    title = "Structure-Corrected RDA1 SNP Loadings",
    shape = "Extreme in both models"
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "rda_loading_distribution.png"),
  p,
  width = 9,
  height = 6,
  dpi = 300
)

saveRDS(RDAfull, file.path(output_dir, "rda_full_model.rds"))
saveRDS(RDA_time, file.path(output_dir, "rda_time_model.rds"))
saveRDS(RDA_time_pc1, file.path(output_dir, "rda_time_pc1_model.rds"))
