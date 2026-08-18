#!/bin/bash
#SBATCH --job-name=msg_mpileup
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=mpileup_%j.log

set -euo pipefail

REFERENCE="$1"
BAM_LIST="$2"
OUTPUT="$3"

[[ -f "$REFERENCE" ]] || { echo "ERROR: Reference not found."; exit 1; }
[[ -f "$BAM_LIST" ]] || { echo "ERROR: BAM list not found."; exit 1; }

command -v samtools >/dev/null || { echo "ERROR: samtools not found."; exit 1; }

mapfile -t BAMS < <(grep -v '^[[:space:]]*$' "$BAM_LIST")
[[ ${#BAMS[@]} -gt 0 ]] || { echo "ERROR: BAM list is empty."; exit 1; }

for BAM in "${BAMS[@]}"; do
    [[ -f "$BAM" ]] || { echo "ERROR: BAM not found: $BAM"; exit 1; }
    samtools quickcheck "$BAM"
done

[[ -f "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"

mkdir -p "$(dirname "$OUTPUT")"

samtools mpileup \
    -f "$REFERENCE" \
    -q 20 \
    -Q 20 \
    -B \
    -d 1000 \
    -b "$BAM_LIST" \
    > "$OUTPUT"

echo "Pools: ${#BAMS[@]}"
echo "Reference positions: $(wc -l < "$OUTPUT")"
