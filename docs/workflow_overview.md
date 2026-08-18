# Workflow overview

## Core workflow

1. **Quality control**
   - FastQC
   - MultiQC
   - NanoPlot

2. **Read processing**
   - Porechop adapter trimming
   - NanoFilt quality/length/end trimming
   - concatenate filtered chunks by sequencing run and barcode

3. **Alignment**
   - minimap2 splice-aware alignment
   - retain primary alignments used for downstream pooled allele counting
   - SAMtools sorting/indexing/statistics

4. **Sample processing**
   - map technical sequencing libraries to biological mosquito pools
   - merge technical replicates
   - recalculate mapping statistics
   - apply pool-level filters
   - establish authoritative downstream sample order

5. **PoPoolation2**
   - SAMtools mpileup
   - mpileup2sync
   - evaluate coverage thresholds
   - filter sync positions
   - snp-frequency-diff
   - pairwise FST

6. **SNP processing**
   - parse PoPoolation2 read-count output
   - retain biallelic SNPs
   - extract sequence context and read-level features
   - build SNP-quality feature table
   - apply final five-flag SNP filter

7. **Population genomics**
   - expected heterozygosity
   - rarefied allelic richness
   - pairwise FST
   - isolation by distance
   - PCoA
   - DAPC
   - RDA / variance partitioning / extreme RDA loadings

## Optional case-study QC

The B2 batch-assessment branch uses featureCounts and DESeq2 PCA to evaluate sequencing-run variation before technical replicate libraries are merged.

## Not included

The public core workflow does not include:

- SNP functional annotation;
- SFS analysis;
- physical-distance pruning;
- B2-only downstream population-genomic analyses;
- automated biological interpretation of collection-period patterns.
