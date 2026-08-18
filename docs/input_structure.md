# Input structure

## Raw ONT layout

The QC and preprocessing scripts assume a run-oriented directory structure similar to:

```text
raw_runs/
├── run01/
│   └── fastq_pass/
│       ├── barcode01/
│       │   ├── chunk_001.fastq.gz
│       │   └── chunk_002.fastq.gz
│       └── unclassified/
└── run02/
    └── fastq_pass/
        └── ...
```

After Porechop/NanoFilt, the workflow uses:

```text
processed/
├── run01/
│   ├── barcode01/
│   │   └── ...
│   └── unclassified/
└── run02/
    └── ...
```

## Sample metadata

One row represents one sequencing library. Multiple rows can map to the same biological `sample_id` when technical replicates exist.

Required fields:

| Column | Meaning |
|---|---|
| `sample_id` | Biological mosquito pool identifier |
| `run_id` | Sequencing run identifier |
| `barcode` | Demultiplexed barcode/directory |
| `pool_size` | Number of mosquitoes in the biological pool |
| `ploidy` | Ploidy; 2 for the case-study mosquitoes |
| `site` | Sampling-site identifier |
| `batch` | Collection-period/group label |
| `site_batch` | Human-readable pool label used in plots |
| `bam_path` | Primary-alignment BAM for this sequencing library |

## Authoritative population order

After biological-pool filtering, `select_samples.py` creates:

- `final_bams.txt`
- `sample_order.tsv`

`final_bams.txt` is supplied to `samtools mpileup`. PoPoolation2 preserves this input order in the sync file, so `sample_order.tsv` is used to map population indices back to biological pool identities for all later analyses.
