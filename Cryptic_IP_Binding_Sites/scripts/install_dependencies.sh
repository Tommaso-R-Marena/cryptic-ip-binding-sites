#!/bin/bash
set -eo pipefail

echo "========================================================="
echo "Installing Cryptic IP Binding Sites Dependencies (Ubuntu)"
echo "========================================================="

sudo apt-get update
sudo apt-get install -y \
    python3 \
    python3-pip \
    fpocket \
    apbs \
    pdb2pqr \
    wget \
    curl \
    git \
    build-essential

echo "System dependencies installed."

pip3 install -r requirements.txt

echo "Python dependencies installed."
echo "Dependencies setup complete!"
echo "Check installation with: python -m pipeline --help"
