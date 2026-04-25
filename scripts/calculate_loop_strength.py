#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import cooler
import bioframe
import cooltools
import os.path

input_cooler_paths=pd.read_csv(sys.argv[1], sep='\t')
input_cooler_path=sys.argv[2]
dataset_name=sys.argv[3]
dots=pd.read_csv(sys.argv[4], sep='\t').loc[:,['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
loop_strengths=pd.read_csv(sys.argv[4], sep='\t')
resolution=int(sys.argv[5])
flank=int(sys.argv[6])
assembly_name=sys.argv[7]
clr_weight_name=sys.argv[8]
nproc=int(sys.argv[9])


def quantify_loops(mtx):
    sq_size = mtx.shape[0]
    midpoint = int(np.floor(sq_size/2))
    mid_9pixels_mean = np.nanmean(mtx[midpoint-1:midpoint+2,midpoint-1:midpoint+2])
    neighboring_size = int(np.ceil(0.3*sq_size) // 2 * 2 + 1) # get the closest odd number of 30% sq_size
    upper_left_mean = np.nanmean(mtx[:neighboring_size,:neighboring_size])
    upper_right_mean = np.nanmean(mtx[:neighboring_size, (sq_size-neighboring_size):])
    lower_right_mean = np.nanmean(mtx[(sq_size-neighboring_size):, (sq_size-neighboring_size):])
    return mid_9pixels_mean/np.nanmean(np.array([upper_left_mean,upper_right_mean,lower_right_mean]))

def quantify_individual_loops_in_stack(stack):
    loop_strength=np.array([])
    for loop_idx in np.arange(stack.shape[-1]):
        loop_strength=np.append(loop_strength, quantify_loops(stack[:,:,loop_idx]))
    return loop_strength


clr = cooler.Cooler(input_cooler_path+f'::resolutions/{resolution}')
# Use bioframe to fetch the genomic features from the UCSC.
chromsizes = bioframe.fetch_chromsizes(assembly_name)
cens = bioframe.fetch_centromeres(assembly_name)
arms = bioframe.make_chromarms(chromsizes, cens)

# Select only chromosomes that are present in the cooler.
arms = arms.set_index("chrom").loc[clr.chromnames].reset_index()

expected = cooltools.expected_cis(
    clr,
    view_df=arms,
    clr_weight_name=clr_weight_name,
    nproc=nproc,)
pileup_mtx=cooltools.pileup(clr, dots, view_df=arms, expected_df=expected, flank=flank, clr_weight_name=clr_weight_name, nproc=nproc)
loop_strengths['loop_strength'] = quantify_individual_loops_in_stack(pileup_mtx)

# Sort loop_strengths by chrom1 and then start1
## Make 'chrom1' a categorical with custom order
loop_strengths['chrom1'] = pd.Categorical(loop_strengths['chrom1'], categories=chromsizes.index, ordered=True)
loop_strengths = loop_strengths.sort_values(by=['chrom1', 'start1']).reset_index(drop=True)

clr_path=input_cooler_paths.loc[input_cooler_paths.name == dataset_name,'path'].values[0]
filename=os.path.splitext(os.path.basename(clr_path))[0] + f'.loop.strength.resolution.{int(resolution/1000)}kb.tsv'
loop_strengths.to_csv(filename, sep='\t', index=False, header=True)

