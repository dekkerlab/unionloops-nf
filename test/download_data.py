#!/usr/bin/env python3

import os
import pandas as pd
import cooltools

# Get the absolute path of the directory containing this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Path to the data directory inside the script directory
data_dir = os.path.join(script_dir, 'data')

# Create the data directory if it doesn't exist
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

# Download files to the test/data/ folder
cool_file_HFF = cooltools.download_data("HFF_MicroC", cache=True, data_dir=data_dir)
cool_file_hESC = cooltools.download_data("hESC_MicroC", cache=True, data_dir=data_dir)

# Prepare the TSV contents with absolute paths
coolers = {
    'name': ['HFF', 'hESC'],
    'path': [cool_file_HFF, cool_file_hESC]
}

# Save to test/test_mcool_paths.tsv
output_tsv = os.path.join(script_dir, 'test_mcool_paths.tsv')
pd.DataFrame(coolers).to_csv(output_tsv, sep='\t', index=False)
