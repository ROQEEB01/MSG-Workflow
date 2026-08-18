#!/bin/bash
#SBATCH --job-name=msg_nanoplot
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=nanoplot_%A_%a.log

set -euo pipefail
shopt -s nullglob

INPUT_DIR="$1"
OUTPUT_DIR="$2"
STAGE="$3"
BARCODE_LIST="$4"

if [[ "$STAGE" != "raw" && "$STAGE" != "filtered" ]]; then
    echo "ERROR: STAGE must be raw or filtered."
    exit 1
fi

BARCODE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BARCODE_LIST")
[[ -n "$BARCODE" ]] || { echo "ERROR: Empty barcode entry."; exit 1; }

command -v NanoPlot >/dev/null || { echo "ERROR: NanoPlot not found."; exit 1; }

mapfile -t RUN_DIRS < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

for RUN_DIR in "${RUN_DIRS[@]}"; do
    RUN_NAME=$(basename "$RUN_DIR")

    if [[ "$STAGE" == "raw" ]]; then
        SAMPLE_DIR="$RUN_DIR/fastq_pass/$BARCODE"
    else
        SAMPLE_DIR="$RUN_DIR/$BARCODE"
    fi

    [[ -d "$SAMPLE_DIR" ]] || continue

    FASTQS=("$SAMPLE_DIR"/*.fastq.gz "$SAMPLE_DIR"/*.fq.gz)
    [[ ${#FASTQS[@]} -gt 0 ]] || continue

    OUT="$OUTPUT_DIR/$STAGE/$RUN_NAME/$BARCODE"
    mkdir -p "$OUT"

    NanoPlot \
        --fastq "${FASTQS[@]}" \
        --threads "${SLURM_CPUS_PER_TASK:-4}" \
        --outdir "$OUT"
done
