#!/bin/bash
#SBATCH --job-name=msg_mm2_index
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=minimap2_index_%j.log

set -euo pipefail

REFERENCE_FASTA="$1"
OUTPUT_INDEX="$2"

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found."; exit 1; }
[[ -f "$REFERENCE_FASTA" ]] || { echo "ERROR: Reference not found."; exit 1; }

mkdir -p "$(dirname "$OUTPUT_INDEX")"
minimap2 -d "$OUTPUT_INDEX" "$REFERENCE_FASTA"
