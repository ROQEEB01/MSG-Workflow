#!/bin/bash

#SBATCH --job-name=msg_featurecounts
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=featurecounts_%j.log

# MSG Workflow
# Generate gene-level counts for sequencing-run batch assessment.
#
# Usage:
#   sbatch run_featurecounts.slurm \
#       /path/to/annotation.gtf \
#       /path/to/bam_list.txt \
#       /path/to/primary_counts.txt

set -euo pipefail

GTF="$1"
BAM_LIST="$2"
OUTPUT="$3"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF file not found: $GTF"
    exit 1
fi

if [[ ! -f "$BAM_LIST" ]]; then
    echo "ERROR: BAM list not found: $BAM_LIST"
    exit 1
fi

mapfile -t BAM_FILES < <(grep -v '^[[:space:]]*$' "$BAM_LIST")

if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: No BAM files found in $BAM_LIST"
    exit 1
fi

for bam in "${BAM_FILES[@]}"; do
    if [[ ! -f "$bam" ]]; then
        echo "ERROR: BAM file not found: $bam"
        exit 1
    fi
done

mkdir -p "$(dirname "$OUTPUT")"

featureCounts \
    -T "$THREADS" \
    -a "$GTF" \
    -o "$OUTPUT" \
    -t exon \
    -g gene_id \
    --maxMOp 200 \
    "${BAM_FILES[@]}"

echo "featureCounts output: $OUTPUT"
echo "BAM files counted: ${#BAM_FILES[@]}"
