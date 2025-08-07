#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import cooler
import bioframe
import cooltools

clr_path=sys.argv[1]
# Available assembly names in bioframe package (e.g., hg38 and mm10)
# Genomic features from the UCSC
# See details in https://bioframe.readthedocs.io/en/latest/guide-io.html#curated-genome-assembly-build-information
assembly_name=sys.argv[2]
resolution=int(sys.argv[3])

clr_weight_name=sys.argv[4]
max_loci_separation=int(sys.argv[5])
max_nans_tolerated=int(sys.argv[6])
lambda_bin_fdr=float(sys.argv[7])
tile_size=int(sys.argv[8])
nproc=int(sys.argv[9])
dataset_name=sys.argv[10]


# Use bioframe to fetch the genomic features from the UCSC.
chromsizes = bioframe.fetch_chromsizes(assembly_name)
cens = bioframe.fetch_centromeres(assembly_name)
arms = bioframe.make_chromarms(chromsizes, cens)

clr=cooler.Cooler(clr_path+f'::resolutions/{resolution}')

# Select only chromosomes that are present in the cooler.
arms = arms.set_index("chrom").loc[clr.chromnames].reset_index()

expected = cooltools.expected_cis(
    clr,
    view_df=arms,
    clr_weight_name=clr_weight_name,
    nproc=nproc,)

# Generate enriched pixels for each matrix
enriched_pixels=cooltools.dots(
    clr,
    expected=expected,
    view_df=arms,
    lambda_bin_fdr=lambda_bin_fdr,
    max_loci_separation=max_loci_separation,
    tile_size=tile_size,
    max_nans_tolerated=max_nans_tolerated,
    # Only output enriched pixels, no clustering is performed.
    clustering_radius=None,
    clr_weight_name=clr_weight_name,
    nproc=nproc,)

filename= dataset_name + f'.enriched.pixels.resolution.{int(resolution/1000)}kb.tsv'
enriched_pixels.to_csv(filename, sep='\t', index=False, header=True)


