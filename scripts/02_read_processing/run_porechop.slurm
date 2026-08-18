#!/bin/bash
#SBATCH --job-name=msg_porechop
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=porechop_%A_%a.log

set -euo pipefail
shopt -s nullglob

INPUT_DIR="$1"
OUTPUT_DIR="$2"
BARCODE_LIST="$3"

BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_LIST")
[[ -n "$BARCODE" ]] || { echo "ERROR: Empty barcode entry."; exit 1; }

command -v porechop >/dev/null || { echo "ERROR: porechop not found."; exit 1; }

mapfile -t RUN_DIRS < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

for RUN_DIR in "${RUN_DIRS[@]}"; do
    RUN_NAME=$(basename "$RUN_DIR")
    SAMPLE_DIR="$RUN_DIR/fastq_pass/$BARCODE"
    [[ -d "$SAMPLE_DIR" ]] || continue

    FASTQS=("$SAMPLE_DIR"/*.fastq.gz "$SAMPLE_DIR"/*.fq.gz)
    [[ ${#FASTQS[@]} -gt 0 ]] || continue

    OUT_DIR="$OUTPUT_DIR/$RUN_NAME/$BARCODE"
    mkdir -p "$OUT_DIR"

    for FASTQ in "${FASTQS[@]}"; do
        BASE=$(basename "$FASTQ")
        BASE=${BASE%.fastq.gz}
        BASE=${BASE%.fq.gz}

        porechop \
            -i "$FASTQ" \
            -o "$OUT_DIR/${BASE}_trim.fastq.gz" \
            --threads "${SLURM_CPUS_PER_TASK:-8}"
    done
done
