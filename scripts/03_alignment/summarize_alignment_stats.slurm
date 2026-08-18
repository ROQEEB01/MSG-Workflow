#!/bin/bash
#SBATCH --job-name=msg_alignment_stats
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=alignment_stats_%j.log

set -euo pipefail

BAM_DIR="$1"
OUTPUT="$2"

command -v samtools >/dev/null || { echo "ERROR: samtools not found."; exit 1; }

mkdir -p "$(dirname "$OUTPUT")"

printf "sample_id\ttotal_reads\tmapped_reads\tunmapped_reads\tmapping_rate\tmean_covered_depth\n" > "$OUTPUT"

mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)

for BAM in "${BAMS[@]}"; do
    SAMPLE_ID=$(basename "$BAM" .bam)

    TOTAL=$(samtools view -c "$BAM")
    MAPPED=$(samtools view -c -F 4 "$BAM")
    UNMAPPED=$(samtools view -c -f 4 "$BAM")

    if [[ "$TOTAL" -gt 0 ]]; then
        RATE=$(awk -v m="$MAPPED" -v t="$TOTAL" 'BEGIN {printf "%.4f", 100*m/t}')
    else
        RATE="0.0000"
    fi

    MEAN_DEPTH=$(samtools depth "$BAM" |
        awk '{sum += $3; n++} END {if (n>0) printf "%.6f", sum/n; else print "0"}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE_ID" "$TOTAL" "$MAPPED" "$UNMAPPED" "$RATE" "$MEAN_DEPTH" \
        >> "$OUTPUT"
done
