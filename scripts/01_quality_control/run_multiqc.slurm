#!/bin/bash
#SBATCH --job-name=msg_multiqc
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=multiqc_%j.log

set -euo pipefail

INPUT_DIR="$1"
OUTPUT_DIR="$2"

command -v multiqc >/dev/null || { echo "ERROR: multiqc not found."; exit 1; }

mkdir -p "$OUTPUT_DIR"
multiqc "$INPUT_DIR" --outdir "$OUTPUT_DIR"
