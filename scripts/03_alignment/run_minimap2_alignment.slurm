#!/bin/bash
#SBATCH --job-name=msg_align
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --output=alignment_%A_%a.log

set -euo pipefail

REFERENCE="$1"
SAMPLE_LIST="$2"
OUTPUT_DIR="$3"

MM2_THREADS=8
SORT_THREADS=8

[[ -f "$REFERENCE" ]] || { echo "ERROR: minimap2 index not found: $REFERENCE"; exit 1; }
[[ -f "$SAMPLE_LIST" ]] || { echo "ERROR: sample list not found: $SAMPLE_LIST"; exit 1; }

FASTQ_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
[[ -n "$FASTQ_FILE" && -f "$FASTQ_FILE" ]] || {
    echo "ERROR: FASTQ not found for task ${SLURM_ARRAY_TASK_ID}: $FASTQ_FILE"
    exit 1
}

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found."; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found."; exit 1; }

SAMPLE_ID=$(basename "$FASTQ_FILE")
SAMPLE_ID=${SAMPLE_ID%.fastq.gz}
SAMPLE_ID=${SAMPLE_ID%.fq.gz}

BAM_DIR="$OUTPUT_DIR/bam"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$BAM_DIR" "$LOG_DIR"

BAM_FILE="$BAM_DIR/${SAMPLE_ID}_primary.bam"

minimap2 \
    -ax splice \
    -t "$MM2_THREADS" \
    "$REFERENCE" \
    "$FASTQ_FILE" \
    2> "$LOG_DIR/${SAMPLE_ID}_minimap2.log" |
samtools view \
    -u \
    -F 2304 \
    - |
samtools sort \
    -@ "$SORT_THREADS" \
    -o "$BAM_FILE" \
    -

samtools index "$BAM_FILE"
samtools quickcheck "$BAM_FILE"
