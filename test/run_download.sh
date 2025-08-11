#!/bin/bash
set -e  # Exit on error

# Inform user about runtime
echo "Please be patient, this process may take a while..."

# Load conda into the shell session
source "$(conda info --base)/etc/profile.d/conda.sh"

# Activate the 'unionloops-nf' conda environment
conda activate unionloops-nf

# Run the Python script
python ./download_data.py
