#!/bin/bash
#SBATCH --job-name=msg_concat
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=concat_%A_%a.log

set -euo pipefail
shopt -s nullglob

INPUT_DIR="$1"
OUTPUT_DIR="$2"
BARCODE_LIST="$3"

BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_LIST")
[[ -n "$BARCODE" ]] || { echo "ERROR: Empty barcode entry."; exit 1; }

mkdir -p "$OUTPUT_DIR"

mapfile -t RUN_DIRS < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

for RUN_DIR in "${RUN_DIRS[@]}"; do
    RUN_NAME=$(basename "$RUN_DIR")
    SAMPLE_DIR="$RUN_DIR/$BARCODE"
    [[ -d "$SAMPLE_DIR" ]] || continue

    FASTQS=("$SAMPLE_DIR"/*_filter.fastq.gz)
    [[ ${#FASTQS[@]} -gt 0 ]] || continue

    OUT="$OUTPUT_DIR/${RUN_NAME}_${BARCODE}.fastq.gz"
    cat "${FASTQS[@]}" > "$OUT"
done
