#!/bin/bash
#SBATCH --job-name=msg_snpfreq
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --output=snp_frequency_%j.log

set -euo pipefail

POPOOLATION2="$1"
SYNC="$2"
OUTPUT_PREFIX="$3"

SCRIPT="$POPOOLATION2/snp-frequency-diff.pl"

[[ -f "$SCRIPT" ]] || { echo "ERROR: snp-frequency-diff.pl not found."; exit 1; }
[[ -f "$SYNC" ]] || { echo "ERROR: Sync file not found."; exit 1; }

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

perl "$SCRIPT" \
    --input "$SYNC" \
    --output-prefix "$OUTPUT_PREFIX" \
    --min-count 4 \
    --min-coverage 5 \
    --max-coverage 500
