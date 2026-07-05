#!/bin/bash
#SBATCH -p workq
#SBATCH -J SSCS_array
#SBATCH --mem=20G
#SBATCH -e logs/error_%A_%a.err
#SBATCH -o logs/output_%A_%a.out
#SBATCH --array=1-4

set -euo pipefail

TABLE="automaticrun.tsv"
PIPELINE="DSPipeline2.sh"

mkdir -p logs

line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TABLE")

alignRef=$(echo "$line" | cut -f1)
Sample=$(echo "$line" | cut -f2)
read1in=$(echo "$line" | cut -f3)
read2in=$(echo "$line" | cut -f4)

echo "Running sample: $Sample"
echo "Read1: $read1in"
echo "Read2: $read2in"
echo "Reference: $alignRef"

mkdir -p "$Sample"

cd "$Sample"

bash "../$PIPELINE" "$alignRef" "$Sample" "$read1in" "$read2in"
