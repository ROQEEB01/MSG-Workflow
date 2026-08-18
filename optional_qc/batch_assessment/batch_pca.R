library(DESeq2)
library(tidyverse)
library(ggrepel)

# Inputs
count_file <- "primary_counts.txt"
metadata_file <- "batch_metadata.tsv"

# Output
output_file <- "pca_mosquitoB2_labeled.pdf"
output_table <- "pca_mosquitoB2_coordinates.tsv"

# Read featureCounts output
fc <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE
)

required_annotation <- c(
  "Geneid", "Chr", "Start", "End", "Strand", "Length"
)

if (!all(required_annotation %in% names(fc))) {
  stop("featureCounts file is missing expected annotation columns.")
}

sample_cols <- setdiff(names(fc), required_annotation)

counts <- as.matrix(fc[, sample_cols, drop = FALSE])
rownames(counts) <- fc$Geneid

# featureCounts may store full BAM paths as column names
colnames(counts) <- basename(colnames(counts))
colnames(counts) <- sub("_primary\\.bam$", "", colnames(counts))

# Metadata should contain:
# sample, biological_sample, batch, collection_period
metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE
) %>%
  filter(collection_period == "B2")

required_metadata <- c(
  "sample",
  "biological_sample",
  "batch",
  "collection_period"
)

if (!all(required_metadata %in% names(metadata))) {
  stop("batch_metadata.tsv is missing required columns.")
}

missing_samples <- setdiff(metadata$sample, colnames(counts))

if (length(missing_samples) > 0) {
  stop(
    "Samples in metadata not found in count matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

# Keep B2 samples in metadata order
counts_b2 <- counts[, metadata$sample, drop = FALSE]

metadata <- metadata %>%
  mutate(
    batch = factor(batch),
    biological_sample = factor(biological_sample)
  )

metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df$sample

# DESeq2 normalization and variance-stabilizing transformation
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_b2),
  colData = metadata_df,
  design = ~ batch + biological_sample
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

vsd <- varianceStabilizingTransformation(
  dds,
  blind = FALSE
)

# PCA
pca_data <- plotPCA(
  vsd,
  intgroup = c("batch", "biological_sample"),
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar"),
  1
)

pca_data <- pca_data %>%
  mutate(
    sample_label = sub(
      ".*barcode",
      "",
      as.character(biological_sample)
    )
  )

write_tsv(
  as_tibble(pca_data),
  output_table
)

p <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = batch
  )
) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = sample_label),
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Sequencing run"
  ) +
  theme_bw()

ggsave(
  output_file,
  p,
  width = 12,
  height = 10
)
