#!/bin/bash
set -e  # Exit on error

# Load conda (ensure conda command is available)
source ~/.bashrc

# Activate the 'unionloops-nf' conda environment
conda activate unionloops-nf

# Run the Python script
./download_data.py
