#!/bin/bash
set -eo pipefail

VALIDATION_DIR="data/structures/validation"
mkdir -p "$VALIDATION_DIR"

echo "Downloading Validation Structures..."

# ADAR2 (AlphaFold)
wget -q -nc -O "$VALIDATION_DIR/AF-P78563-F1-model_v4.pdb" "https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb" || echo "AF-P78563-F1-model_v4.pdb already exists"

# Positives
wget -q -nc -O "$VALIDATION_DIR/1ZY7.pdb" "https://files.rcsb.org/download/1ZY7.pdb" || echo "1ZY7.pdb already exists"
wget -q -nc -O "$VALIDATION_DIR/5HDT.pdb" "https://files.rcsb.org/download/5HDT.pdb" || echo "5HDT.pdb already exists"
wget -q -nc -O "$VALIDATION_DIR/5ICN.pdb" "https://files.rcsb.org/download/5ICN.pdb" || echo "5ICN.pdb already exists"
wget -q -nc -O "$VALIDATION_DIR/4A69.pdb" "https://files.rcsb.org/download/4A69.pdb" || echo "4A69.pdb already exists"

# Negatives
wget -q -nc -O "$VALIDATION_DIR/1MAI.pdb" "https://files.rcsb.org/download/1MAI.pdb" || echo "1MAI.pdb already exists"
wget -q -nc -O "$VALIDATION_DIR/1BTK.pdb" "https://files.rcsb.org/download/1BTK.pdb" || echo "1BTK.pdb already exists"
wget -q -nc -O "$VALIDATION_DIR/1FAO.pdb" "https://files.rcsb.org/download/1FAO.pdb" || echo "1FAO.pdb already exists"
wget -q -nc -O "$VALIDATION_DIR/1FGY.pdb" "https://files.rcsb.org/download/1FGY.pdb" || echo "1FGY.pdb already exists"

echo "Validation structures downloaded to $VALIDATION_DIR"
