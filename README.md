# MSG workflow

**MSG workflow** is a reproducible bioinformatics workflow for repurposing long-read RNA sequencing data generated for arbovirus surveillance for mosquito population-genomic analysis.

The workflow was developed using pooled Oxford Nanopore direct-RNA sequencing from the *Culex pipiens/restuans* surveillance group as a case study. It combines long-read RNA preprocessing, splice-aware alignment, pooled allele-frequency estimation, SNP-level technical quality assessment, and population-genomic analyses.

> **Scope:** MSG is intended for exploratory population-genomic analysis of host-derived reads recovered from surveillance RNA-seq. RNA-derived allele frequencies can be influenced by transcript abundance, allele-specific expression, mapping bias, RNA editing, sequencing error, and unequal pool size. High-quality SNP filtering reduces identifiable technical risk but does not prove that every retained SNP is a genomic polymorphism.

## Workflow

```text
Raw ONT FASTQ
    │
    ▼
01_quality_control
    │
    ▼
02_read_processing
    │
    ▼
03_alignment
    │
    ├──────────── optional_qc/batch_assessment
    │
    ▼
04_sample_processing
    │
    ▼
05_popoolation2
    │
    ▼
06_snp_processing
    │
    ▼
07_population_genomics
```

The core workflow uses all biological pools that pass the pool-level filters. Collection period is retained as metadata and used where required by the downstream analyses; the workflow does not subset the final analysis to a single collection period.

## Repository structure

```text
MSG-workflow/
├── README.md
├── .gitignore
├── config/
│   ├── barcodes.txt
│   ├── sample_metadata_template.tsv
│   └── site_coordinates_template.csv
├── scripts/
│   ├── 01_quality_control/
│   ├── 02_read_processing/
│   ├── 03_alignment/
│   ├── 04_sample_processing/
│   ├── 05_popoolation2/
│   ├── 06_snp_processing/
│   └── 07_population_genomics/
├── optional_qc/
│   └── batch_assessment/
├── metadata/
├── docs/
└── example_data/
```

## Core inputs

### `config/barcodes.txt`

One barcode or demultiplexing directory name per line. Include `unclassified` if it is present in the raw sequencing output and should be processed during read QC.

### `sample_metadata.tsv`

The sample metadata table links sequencing libraries to biological mosquito pools. Required columns are:

```text
sample_id
run_id
barcode
pool_size
ploidy
site
batch
site_batch
bam_path
```

For the *Culex* case study, `ploidy = 2`. `pool_size` is the number of mosquitoes in the biological pool. The workflow derives:

```text
chromosome_count = pool_size × ploidy
```

The chromosome count is used for pooled finite-sample corrections and the PoPoolation2 FST pool-size vector.

### `site_coordinates.csv`

Required for isolation-by-distance and RDA:

```text
SiteID,Lat,Long
```

## Running the workflow

The examples below illustrate the intended order. Replace paths with your project paths and use the scheduler/module configuration appropriate for your computing environment.

### 1. Raw-read QC

```bash
N=$(wc -l < config/barcodes.txt)

sbatch --array=1-"$N"%20 \
  scripts/01_quality_control/run_fastqc.slurm \
  /path/to/raw_runs \
  results/qc/fastqc \
  raw \
  config/barcodes.txt

sbatch scripts/01_quality_control/run_multiqc.slurm \
  results/qc/fastqc/raw \
  results/qc/multiqc/raw

sbatch --array=1-"$N"%20 \
  scripts/01_quality_control/run_nanoplot.slurm \
  /path/to/raw_runs \
  results/qc/nanoplot \
  raw \
  config/barcodes.txt
```

### 2. Adapter trimming and read filtering

```bash
sbatch --array=1-"$N"%20 \
  scripts/02_read_processing/run_porechop.slurm \
  /path/to/raw_runs \
  results/porechop \
  config/barcodes.txt

sbatch --array=1-"$N"%20 \
  scripts/02_read_processing/run_nanofilt.slurm \
  results/porechop \
  results/nanofilt \
  config/barcodes.txt

sbatch --array=1-"$N"%20 \
  scripts/02_read_processing/run_concatenate_reads.slurm \
  results/nanofilt \
  results/concatenated \
  config/barcodes.txt
```

The filtering parameters used in the case study were:

```text
mean read quality ≥ 11
read length 300–15,000 bp
trim 20 bp from the 5′ end
trim 20 bp from the 3′ end
```

Post-filtering FastQC/MultiQC/NanoPlot can be run with the same QC scripts using `filtered` mode.

### 3. Alignment

Build a minimap2 index:

```bash
sbatch scripts/03_alignment/build_minimap2_index.slurm \
  reference.fna \
  reference.mmi
```

Create a FASTQ list:

```bash
find results/concatenated \
  -maxdepth 1 \
  -name "*.fastq.gz" \
  -type f | sort > concatenated_files.txt
```

Submit the alignment array:

```bash
N=$(wc -l < concatenated_files.txt)

sbatch --array=1-"$N"%20 \
  scripts/03_alignment/run_minimap2_alignment.slurm \
  reference.mmi \
  concatenated_files.txt \
  results/alignment
```

Alignment uses the splice-aware minimap2 preset:

```text
-ax splice
```

### 4. Merge technical replicates and select biological pools

Prepare the merge manifest:

```bash
python scripts/04_sample_processing/prepare_merge_manifest.py \
  --metadata sample_metadata.tsv \
  --output-dir results/sample_processing
```

Merge technical replicates:

```bash
N=$(($(wc -l < results/sample_processing/merge_manifest.tsv) - 1))

sbatch --array=1-"$N"%20 \
  scripts/04_sample_processing/run_merge_technical_replicates.slurm \
  results/sample_processing/merge_manifest.tsv \
  results/merged_bams
```

Select biological pools:

```bash
python scripts/04_sample_processing/select_samples.py \
  --manifest results/sample_processing/merge_manifest.tsv \
  --stats-dir results/merged_bams/stats \
  --output-dir results/selected_samples
```

Default case-study pool-level filters:

```text
mapped reads ≥ 50,000
mapping rate ≥ 50%
pool size ≥ 5 mosquitoes
```

`final_bams.txt` is the authoritative BAM order for mpileup. `sample_order.tsv` preserves that same population order for every downstream analysis.

### Optional: sequencing-run batch assessment

For datasets with technical replicate libraries across sequencing runs, see:

```text
optional_qc/batch_assessment/
```

The case study used B2 libraries sequenced across multiple runs. Gene-level counts were generated with featureCounts and a DESeq2 PCA was used to assess run-associated variation before technical replicates were merged.

### 5. PoPoolation2

Generate mpileup:

```bash
sbatch scripts/05_popoolation2/run_mpileup.slurm \
  reference.fna \
  results/selected_samples/final_bams.txt \
  results/popoolation2/all_samples.mpileup
```

Convert to sync:

```bash
sbatch scripts/05_popoolation2/run_mpileup2sync.slurm \
  /path/to/popoolation2 \
  results/popoolation2/all_samples.mpileup \
  results/popoolation2/all_samples.sync
```

Assess alternative coverage tiers:

```bash
python scripts/05_popoolation2/assess_sync_coverage.py \
  --sync results/popoolation2/all_samples.sync \
  --output results/popoolation2/coverage_tiers.tsv
```

The case study selected the tier requiring at least 5× coverage in every retained pool:

```bash
python scripts/05_popoolation2/filter_sync_coverage.py \
  --sync results/popoolation2/all_samples.sync \
  --output results/popoolation2/tier1_all_pools_5x.sync \
  --min-coverage 5 \
  --required-fraction 1
```

Estimate SNP frequencies:

```bash
sbatch scripts/05_popoolation2/run_snp_frequency_diff.slurm \
  /path/to/popoolation2 \
  results/popoolation2/tier1_all_pools_5x.sync \
  results/popoolation2/snp_frequency
```

Calculate PoPoolation2 FST:

```bash
sbatch scripts/05_popoolation2/run_fst.slurm \
  /path/to/popoolation2 \
  results/popoolation2/tier1_all_pools_5x.sync \
  results/selected_samples/sample_order.tsv \
  results/popoolation2/tier1_sync.fst
```

### 6. Build and quality-filter the SNP dataset

Parse the PoPoolation2 `_rc` file:

```bash
python scripts/06_snp_processing/build_snp_summary.py \
  --rc results/popoolation2/snp_frequency_rc \
  --sample-order results/selected_samples/sample_order.tsv \
  --output results/snp_processing/snp_summary.tsv
```

Retain biallelic SNPs:

```bash
python scripts/06_snp_processing/filter_biallelic_snps.py \
  --input results/snp_processing/snp_summary.tsv \
  --output results/snp_processing/per_pool_AB_table.tsv
```

Create SNP positions for read-level auditing:

```bash
python scripts/06_snp_processing/prepare_snp_positions.py \
  --snps results/snp_processing/per_pool_AB_table.tsv \
  --output results/snp_processing/snp_positions.bed
```

Extract sequence context:

```bash
python scripts/06_snp_processing/extract_snp_context.py \
  --reference reference.fna \
  --bed results/snp_processing/snp_positions.bed \
  --output results/snp_processing/snp_context.tsv
```

Extract read-level SNP features from each final BAM:

```bash
N=$(wc -l < results/selected_samples/final_bams.txt)

sbatch --array=1-"$N"%20 \
  scripts/06_snp_processing/run_extract_read_features.slurm \
  results/selected_samples/final_bams.txt \
  results/snp_processing/snp_positions.bed \
  results/snp_processing/read_features
```

Build the SNP audit feature table:

```bash
Rscript scripts/06_snp_processing/build_snp_quality_features.R \
  results/snp_processing/per_pool_AB_table.tsv \
  results/selected_samples/sample_order.tsv \
  results/snp_processing/snp_context.tsv \
  results/snp_processing/read_features \
  results/snp_processing/audit
```

Apply the final SNP filter:

```bash
Rscript scripts/06_snp_processing/filter_snps.R \
  results/snp_processing/audit/snp_quality_features.rds \
  results/snp_processing/per_pool_AB_table.tsv \
  results/selected_samples/sample_order.tsv \
  results/snp_processing/final
```

The final filter flags:

1. mean depth < 30;
2. depth CV above the 95th percentile;
3. B3:B1 depth ratio outside the 2.5th–97.5th percentile range;
4. homopolymer context;
5. read-end enrichment above the 95th percentile.

SNPs with zero flags are retained as high quality for the primary analyses. Low MAF, strand imbalance, and substitution classes compatible with RNA editing can be assessed during the audit but are not final retention criteria.

### 7. Population-genomic analyses

Expected heterozygosity:

```bash
Rscript scripts/07_population_genomics/genetic_diversity.R \
  results/snp_processing/final/af_high_quality.rds \
  results/snp_processing/final/depth_high_quality.rds \
  results/selected_samples/sample_order.tsv \
  results/population_genomics/genetic_diversity
```

Allelic richness:

```bash
Rscript scripts/07_population_genomics/allelic_richness.R \
  results/snp_processing/final/af_high_quality.rds \
  results/selected_samples/sample_order.tsv \
  results/population_genomics/allelic_richness
```

Pairwise FST:

```bash
Rscript scripts/07_population_genomics/fst_analysis.R \
  results/snp_processing/final/snps_high_quality.tsv \
  results/popoolation2/tier1_sync.fst \
  results/selected_samples/sample_order.tsv \
  results/population_genomics/fst
```

Isolation by distance:

```bash
Rscript scripts/07_population_genomics/isolation_by_distance.R \
  results/population_genomics/fst/fst_high_quality.rds \
  results/selected_samples/sample_order.tsv \
  site_coordinates.csv \
  results/population_genomics/ibd
```

PCoA:

```bash
Rscript scripts/07_population_genomics/pcoa.R \
  results/population_genomics/fst/fst_high_quality.rds \
  results/selected_samples/sample_order.tsv \
  results/population_genomics/pcoa
```

DAPC:

```bash
Rscript scripts/07_population_genomics/dapc.R \
  results/snp_processing/final/af_high_quality.rds \
  results/selected_samples/sample_order.tsv \
  results/population_genomics/dapc
```

RDA:

```bash
Rscript scripts/07_population_genomics/rda.R \
  results/snp_processing/final/af_high_quality.rds \
  results/selected_samples/sample_order.tsv \
  site_coordinates.csv \
  results/population_genomics/rda
```

The RDA workflow performs PCA-derived background-structure estimation, forward selection, full and partial RDA models, and identification of SNPs with extreme RDA1 loadings using a ±2.5 SD threshold. The final table contains SNPs retained in both the standard collection-period RDA and the PC1-conditioned sensitivity analysis.


## Software

See [`docs/software_versions.md`](docs/software_versions.md).


## License

MIT License — Copyright © 2026 Roqeeb Akinbile
