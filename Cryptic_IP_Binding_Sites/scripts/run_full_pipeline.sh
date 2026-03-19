#!/bin/bash
set -eo pipefail

echo "============================================="
echo "Running Full Cryptic IP Pipeline"
echo "============================================="

echo "Step 1: Download Validation Structures"
bash scripts/download_validation_structures.sh

echo "Step 2: Run Validation Checks"
python -m validation.validation_report

# Check status of validation
if [ $? -eq 0 ]; then
    echo "Validation PASSED. Proceeding to Screening."
else
    echo "Validation FAILED. Aborting full pipeline."
    exit 1
fi

echo "Step 3: Screening (Example: Yeast)"
# In a real local run, we would fetch yeast first, but typically this happens in CI.
# Example command for local execution:
# python -m screening.screen_proteome --organism=yeast --workers=4

echo "Pipeline script finished."
