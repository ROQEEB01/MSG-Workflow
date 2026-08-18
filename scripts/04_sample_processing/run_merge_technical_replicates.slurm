#!/bin/bash
#SBATCH --job-name=msg_merge_bams
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=merge_bams_%A_%a.log

set -euo pipefail

MANIFEST="$1"
OUTPUT_DIR="$2"

[[ -f "$MANIFEST" ]] || { echo "ERROR: Manifest not found."; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools not found."; exit 1; }

ROW=$(awk -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR==n {print; exit}' "$MANIFEST")
[[ -n "$ROW" ]] || { echo "ERROR: No manifest row for array task."; exit 1; }

IFS=$'\t' read -r \
    SAMPLE_ID POOL_SIZE PLOIDY CHROMOSOME_COUNT SITE BATCH SITE_BATCH N_BAMS BAM_LIST \
    <<< "$ROW"

[[ -f "$BAM_LIST" ]] || { echo "ERROR: BAM list not found: $BAM_LIST"; exit 1; }

mapfile -t BAMS < <(grep -v '^[[:space:]]*$' "$BAM_LIST")
[[ ${#BAMS[@]} -eq "$N_BAMS" ]] || { echo "ERROR: BAM count mismatch."; exit 1; }

for BAM in "${BAMS[@]}"; do
    [[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM"; exit 1; }
    samtools quickcheck "$BAM"
done

BAM_DIR="$OUTPUT_DIR/bam"
STATS_DIR="$OUTPUT_DIR/stats"
mkdir -p "$BAM_DIR" "$STATS_DIR"

OUT_BAM="$BAM_DIR/${SAMPLE_ID}.bam"

if [[ ${#BAMS[@]} -eq 1 ]]; then
    cp "${BAMS[0]}" "$OUT_BAM"
else
    samtools merge \
        -f \
        -@ "${SLURM_CPUS_PER_TASK:-8}" \
        "$OUT_BAM" \
        "${BAMS[@]}"
fi

samtools index "$OUT_BAM"
samtools quickcheck "$OUT_BAM"

TOTAL=$(samtools view -c "$OUT_BAM")
MAPPED=$(samtools view -c -F 4 "$OUT_BAM")
UNMAPPED=$(samtools view -c -f 4 "$OUT_BAM")

if [[ "$TOTAL" -gt 0 ]]; then
    MAPPING_RATE=$(awk -v m="$MAPPED" -v t="$TOTAL" 'BEGIN {printf "%.4f", 100*m/t}')
else
    MAPPING_RATE="0.0000"
fi

SUMMARY="$STATS_DIR/${SAMPLE_ID}_mapping_summary.tsv"

printf "sample_id\tpool_size\tploidy\tchromosome_count\tsite\tbatch\tsite_batch\ttotal_reads\tmapped_reads\tunmapped_reads\tmapping_rate\tbam_path\n" > "$SUMMARY"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE_ID" "$POOL_SIZE" "$PLOIDY" "$CHROMOSOME_COUNT" \
    "$SITE" "$BATCH" "$SITE_BATCH" "$TOTAL" "$MAPPED" "$UNMAPPED" \
    "$MAPPING_RATE" "$OUT_BAM" >> "$SUMMARY"

samtools flagstat "$OUT_BAM" > "$STATS_DIR/${SAMPLE_ID}_flagstat.txt"
