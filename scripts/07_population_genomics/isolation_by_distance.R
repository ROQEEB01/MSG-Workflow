#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(geosphere))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript isolation_by_distance.R ",
    "fst_high_quality.rds sample_order.tsv site_coordinates.csv output_dir"
  )
}

fst_file <- args[1]
sample_order_file <- args[2]
coords_file <- args[3]
output_dir <- args[4]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fst_matrix <- readRDS(fst_file)
sample_order <- read_tsv(sample_order_file, show_col_types = FALSE) %>%
  arrange(pool_index)
coords <- read_csv(coords_file, show_col_types = FALSE)

if (!all(c("site_batch", "site") %in% names(sample_order))) {
  stop("sample_order.tsv must contain site_batch and site.")
}
if (!all(c("SiteID", "Lat", "Long") %in% names(coords))) {
  stop("site_coordinates.csv must contain SiteID, Lat, and Long.")
}

pool_names <- sample_order$site_batch
if (!all(pool_names %in% rownames(fst_matrix))) {
  stop("FST matrix and sample_order.tsv contain different pool labels.")
}
fst_matrix <- fst_matrix[pool_names, pool_names, drop = FALSE]

pool_coords <- sample_order %>%
  select(site_batch, site) %>%
  left_join(coords, by = c("site" = "SiteID"))
if (anyNA(pool_coords$Lat) || anyNA(pool_coords$Long)) {
  stop("Coordinates are missing for one or more sampling sites.")
}

coord_matrix <- as.matrix(pool_coords %>% select(Long, Lat))
rownames(coord_matrix) <- pool_coords$site_batch
geo_matrix <- distm(coord_matrix, fun = distHaversine) / 1000
rownames(geo_matrix) <- pool_names
colnames(geo_matrix) <- pool_names

set.seed(123)
mantel_result <- mantel(
  as.dist(fst_matrix),
  as.dist(geo_matrix),
  method = "pearson",
  permutations = 999
)

write_tsv(
  tibble(
    n_pools = length(pool_names),
    mantel_r = unname(mantel_result$statistic),
    p_value = mantel_result$signif,
    permutations = 999
  ),
  file.path(output_dir, "mantel_results.tsv")
)

upper <- upper.tri(fst_matrix)
pairs <- tibble(
  pool1 = rownames(fst_matrix)[row(fst_matrix)[upper]],
  pool2 = colnames(fst_matrix)[col(fst_matrix)[upper]],
  geographic_distance_km = geo_matrix[upper],
  FST = fst_matrix[upper]
)
write_tsv(pairs, file.path(output_dir, "ibd_pairwise_values.tsv"))

p_ibd <- ggplot(pairs, aes(x = geographic_distance_km, y = FST)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed") +
  labs(
    x = "Geographic distance (km)",
    y = "Pairwise FST",
    title = "Isolation by Distance",
    subtitle = paste0(
      "Mantel r = ", round(mantel_result$statistic, 3),
      ", p = ", signif(mantel_result$signif, 3)
    )
  ) +
  theme_bw()

ggsave(
  file.path(output_dir, "isolation_by_distance.png"),
  p_ibd, width = 8, height = 6, dpi = 300
)
