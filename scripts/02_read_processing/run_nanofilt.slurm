#!/bin/bash
#SBATCH --job-name=msg_nanofilt
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=nanofilt_%A_%a.log

set -euo pipefail
shopt -s nullglob

INPUT_DIR="$1"
OUTPUT_DIR="$2"
BARCODE_LIST="$3"

MIN_Q=11
MIN_LEN=300
MAX_LEN=15000
HEADCROP=20
TAILCROP=20

BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_LIST")
[[ -n "$BARCODE" ]] || { echo "ERROR: Empty barcode entry."; exit 1; }

command -v NanoFilt >/dev/null || { echo "ERROR: NanoFilt not found."; exit 1; }

mapfile -t RUN_DIRS < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

for RUN_DIR in "${RUN_DIRS[@]}"; do
    RUN_NAME=$(basename "$RUN_DIR")
    SAMPLE_DIR="$RUN_DIR/$BARCODE"
    [[ -d "$SAMPLE_DIR" ]] || continue

    FASTQS=("$SAMPLE_DIR"/*_trim.fastq.gz)
    [[ ${#FASTQS[@]} -gt 0 ]] || continue

    OUT_DIR="$OUTPUT_DIR/$RUN_NAME/$BARCODE"
    mkdir -p "$OUT_DIR"

    for FASTQ in "${FASTQS[@]}"; do
        BASE=$(basename "$FASTQ" _trim.fastq.gz)

        gzip -dc "$FASTQ" |
        NanoFilt \
            -q "$MIN_Q" \
            -l "$MIN_LEN" \
            --maxlength "$MAX_LEN" \
            --headcrop "$HEADCROP" \
            --tailcrop "$TAILCROP" |
        gzip > "$OUT_DIR/${BASE}_filter.fastq.gz"
    done
done
