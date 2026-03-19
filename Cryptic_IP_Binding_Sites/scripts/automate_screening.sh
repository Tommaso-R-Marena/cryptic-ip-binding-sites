#!/bin/bash
set -eo pipefail

echo "============================================================"
echo "    Cryptic IP Binding Sites: End-to-End Automation"
echo "============================================================"

if [ -z "$1" ]; then
    echo "Usage: ./scripts/automate_screening.sh <organism_name>"
    echo "Example: ./scripts/automate_screening.sh yeast"
    exit 1
fi

ORGANISM=$1
echo "-> Starting fully automated pipeline for organism: $ORGANISM"

# 1. Fetch data (assuming fetch script exists or structures are pre-downloaded)
# If fetching is required, the user should rely on .github workflows, 
# but we ensure directories exist.
mkdir -p data/structures/${ORGANISM}
mkdir -p data/results/${ORGANISM}
mkdir -p logs

# 2. Screening the entire proteome
echo "-> [1/2] Screening the $ORGANISM proteome..."
python -m screening.screen_proteome --organism ${ORGANISM} 2>&1 | tee logs/screen_${ORGANISM}.log

# 3. Triage results
echo "-> [2/2] Triaging results to find hit candidates..."
python -m screening.triage --organism ${ORGANISM} 2>&1 | tee logs/triage_${ORGANISM}.log

# 4. Success Output
echo "============================================================"
echo "-> Pipeline Complete. Top candidates are available in:"
echo "   data/results/${ORGANISM}/${ORGANISM}_top_candidates.csv"
echo "-> Visualizations maps are available at:"
echo "   data/results/${ORGANISM}/${ORGANISM}_hits_scatter.png"
echo "============================================================"
