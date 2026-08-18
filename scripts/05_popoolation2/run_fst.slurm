#!/bin/bash
#SBATCH --job-name=msg_fst
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=fst_%j.log

set -euo pipefail

POPOOLATION2="$1"
SYNC="$2"
SAMPLE_ORDER="$3"
OUTPUT="$4"

SCRIPT="$POPOOLATION2/fst-sliding.pl"

[[ -f "$SCRIPT" ]] || { echo "ERROR: fst-sliding.pl not found."; exit 1; }
[[ -f "$SYNC" ]] || { echo "ERROR: Sync file not found."; exit 1; }
[[ -f "$SAMPLE_ORDER" ]] || { echo "ERROR: sample_order.tsv not found."; exit 1; }

SYNC_POOLS=$(awk 'NR==1 {print NF-3; exit}' "$SYNC")
ORDER_POOLS=$(awk 'NR>1 && NF>0 {n++} END {print n+0}' "$SAMPLE_ORDER")

if [[ "$SYNC_POOLS" -ne "$ORDER_POOLS" ]]; then
    echo "ERROR: Pool-count mismatch between sync and sample_order.tsv."
    exit 1
fi

CHROM_COL=$(awk -F'\t' '
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="chromosome_count") {print i; exit}
    }
}' "$SAMPLE_ORDER")

[[ -n "$CHROM_COL" ]] || {
    echo "ERROR: chromosome_count column not found."
    exit 1
}

POOL_SIZES=$(awk -F'\t' -v c="$CHROM_COL" '
NR>1 {
    if ($c !~ /^[0-9]+$/ || $c < 1) {
        printf "ERROR: Invalid chromosome_count on line %d\n", NR > "/dev/stderr"
        exit 1
    }
    values = values (values ? ":" : "") $c
}
END {print values}
' "$SAMPLE_ORDER")

mkdir -p "$(dirname "$OUTPUT")"

perl "$SCRIPT" \
    --input "$SYNC" \
    --output "$OUTPUT" \
    --suppress-noninformative \
    --min-count 4 \
    --min-coverage 5 \
    --max-coverage 500 \
    --min-covered-fraction 1 \
    --window-size 1 \
    --step-size 1 \
    --pool-size "$POOL_SIZES"
