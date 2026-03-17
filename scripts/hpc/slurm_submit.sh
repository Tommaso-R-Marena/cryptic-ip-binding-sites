#!/usr/bin/env bash
set -euo pipefail

PROTEOME_FILE=${1:-data/proteomes.txt}
OUTPUT_DIR=${2:-results/hpc}

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=cryptic-ip
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=${OUTPUT_DIR}/slurm-%j.out

mkdir -p ${OUTPUT_DIR}
python scripts/phase3_screen_proteome.py --proteome-list ${PROTEOME_FILE} --workers 16 --output-dir ${OUTPUT_DIR}
EOF
