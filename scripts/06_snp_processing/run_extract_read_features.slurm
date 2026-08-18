#!/bin/bash

#SBATCH --job-name=msg_snp_read_features
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=snp_read_features_%A_%a.log

# Usage:
#   sbatch --array=1-N%20 run_extract_read_features.slurm \
#       final_bams.txt snp_positions.bed /path/to/output

set -euo pipefail

BAM_LIST="$1"
BED_FILE="$2"
OUTPUT_DIR="$3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "$BAM_LIST" ]]; then
    echo "ERROR: BAM list not found: $BAM_LIST"
    exit 1
fi

if [[ ! -f "$BED_FILE" ]]; then
    echo "ERROR: BED file not found: $BED_FILE"
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: Submit this script as a SLURM array."
    exit 1
fi

if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 not found."
    exit 1
fi

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAM_LIST")
if [[ -z "$BAM" || ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found for task ${SLURM_ARRAY_TASK_ID}: $BAM"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

SAMPLE_ID=$(basename "$BAM" .bam)

python3 "$SCRIPT_DIR/extract_read_features.py" \
    --bam "$BAM" \
    --bed "$BED_FILE" \
    --output "$OUTPUT_DIR/${SAMPLE_ID}_read_features.tsv"
