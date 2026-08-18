# Provenance checks before public release

The workflow is assembled from the thesis methods and the original analysis scripts. A few implementation details should be checked against the final analysis environment before creating a versioned release.

## 1. SAMtools version

The final thesis methods record SAMtools 1.22.1, while an earlier alignment script loaded SAMtools 1.20. Confirm the version used to create the final alignment BAMs.

## 2. NanoPlot version

The final thesis methods record NanoPlot 1.41.1; an earlier workflow note recorded 1.41.0. Confirm the final executed version.

## 3. Primary-alignment filter

The development alignment code removed secondary and supplementary alignments (`-F 2304`). The thesis prose explicitly mentions removal of secondary alignments (`-F 256`). Confirm which command produced the BAMs used for the final analysis and make the prose/repository agree.

## 4. RDA input transformation

The supplied final RDA code starts from an already-created object named `Y_hel`. The current public `rda.R` reconstructs this object with a Hellinger transformation of the high-quality allele-frequency matrix. Verify that this is how `Y_hel` was generated in the final analysis. If not, replace that single input-transformation step with the original calculation before release.

These checks are provenance details, not requests to rerun the biological analyses unless the verified final commands differ from the current reconstruction.
