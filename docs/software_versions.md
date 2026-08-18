# Software versions

The workflow preserves the software choices used in the thesis case study wherever they are known.

| Software | Version / note |
|---|---|
| FastQC | 0.12.0 |
| MultiQC | 1.25.2 |
| NanoPlot | 1.41.1 in final thesis methods; verify against final execution environment |
| Porechop | 0.2.4 |
| NanoFilt | 2.8.0 |
| minimap2 | 2.28 |
| SAMtools | 1.22.1 in final thesis methods; an earlier alignment script loaded 1.20 |
| Subread / featureCounts | 2.0.6 in the batch-effect script |
| PoPoolation2 | version used in thesis workflow; verify exact release string before tagging |
| R | 4.5.2 in final thesis methods |
| Python | Python 3 |
| pysam | required for SNP read/context extraction |

## R packages

Core downstream scripts use:

- tidyverse
- ggrepel
- vegan
- geosphere
- adegenet

Optional batch assessment uses:

- DESeq2
- tidyverse
- ggrepel

## Reproducibility policy

The repository preserves the actual thesis-era tools rather than silently replacing archived software with newer alternatives. For example, Porechop and NanoFilt are retained because they were used in the original workflow.

Modern alternatives can be documented separately without changing the commands used to reproduce the thesis analysis.
