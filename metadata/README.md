# Metadata

Project-specific sample metadata are not included in the public repository.

Copy `config/sample_metadata_template.tsv` and populate one row per sequencing library. Technical replicate libraries should share the same `sample_id`, `pool_size`, `ploidy`, `site`, `batch`, and `site_batch`.

`pool_size` is the number of mosquitoes represented in the biological pool. `chromosome_count` is derived later as:

```text
pool_size × ploidy
```

For the *Culex* case study, ploidy was 2.

The `batch` column is the collection-period/grouping variable used by the downstream case-study analyses. The public workflow keeps this column explicit rather than parsing biological meaning from filenames.
