# Optional sequencing-run batch assessment

This branch reproduces the case-study assessment performed before technical replicate libraries were merged.

Only libraries that have technical replication across sequencing runs need to be included. In the case study, this was the B2 collection.

## 1. Generate gene-level counts

Create a BAM list containing the pre-merge primary BAMs:

```bash
sbatch run_featurecounts.slurm \
  annotation.gtf \
  b2_bams.txt \
  primary_counts.txt
```

The case-study featureCounts settings were:

```text
feature type: exon
grouping attribute: gene_id
--maxMOp 200
```

## 2. Prepare batch metadata

Copy `batch_metadata_template.tsv` and record:

```text
sample
biological_sample
batch
collection_period
```

## 3. Run PCA

```bash
Rscript batch_pca.R
```

The final case-study figure was a B2 PCA with points coloured by sequencing run and labelled by biological mosquito pool.

This is an optional QC branch and does not alter the final population-genomic analysis unless it identifies a reason not to merge technical replicate libraries.
