#!/bin/bash
#SBATCH --job-name=msg_mpileup2sync
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=mpileup2sync_%j.log

set -euo pipefail

POPOOLATION2="$1"
MPILEUP="$2"
OUTPUT="$3"

SCRIPT="$POPOOLATION2/mpileup2sync.pl"

[[ -f "$SCRIPT" ]] || { echo "ERROR: mpileup2sync.pl not found."; exit 1; }
[[ -f "$MPILEUP" ]] || { echo "ERROR: mpileup not found."; exit 1; }

mkdir -p "$(dirname "$OUTPUT")"

perl "$SCRIPT" \
    --fastq-type sanger \
    --min-qual 20 \
    --input "$MPILEUP" \
    --output "$OUTPUT"

POOLS=$(awk 'NR==1 {print NF-3; exit}' "$OUTPUT")
echo "Pools in sync: $POOLS"
