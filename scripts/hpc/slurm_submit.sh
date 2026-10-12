#!/usr/bin/env bash
# SUPERSEDED by the harness in cryptic_ip/hpc/ and scripts/rockfish.py; see
# docs/ROCKFISH.md. This submits one proteome-screen job with hardcoded
# resources and no sharding, census or resume. Kept because it is what the
# proteome screen was run with, not because it is the way to run anything now.
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
