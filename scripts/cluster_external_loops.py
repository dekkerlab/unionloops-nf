#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import bioframe
import cooler
import cooltools
from cooltools.api.dotfinder import *

# Annotate with bin1_id, bin2_id, and balanced.count for all pixels in all conditions
# bin_id are same for all conditions with same resolution
def annotate_balanced(enriched_pixels, clr_dict):
    # bin1_id, bin2_id
    df=enriched_pixels.copy()
    clr_for_bin_id=cooler.Cooler(list(clr_dict.values())[0])
    df['bin1_id']=df.apply(lambda x: clr_for_bin_id.offset((x.chrom1, x.start1, x.end1)), axis=1)
    df['bin2_id']=df.apply(lambda x: clr_for_bin_id.offset((x.chrom2, x.start2, x.end2)), axis=1)
    
    for condition, clr_file in clr_dict.items():
        clr=cooler.Cooler(clr_file)
        df[condition+'_balanced']=df.apply(lambda x: clr.matrix(balance=True)[x.bin1_id, x.bin2_id].item(), axis=1)
    return df



def clustering_step_across_conditions(
    scored_df,
    dots_clustering_radius=20_000,
    assigned_regions_name="region",
):

    # cluster within each regions separately and accumulate the result:
    pixel_clust_list = []
    scored_pixels_by_region = scored_df.groupby(assigned_regions_name, observed=True)
    for region, _df in scored_pixels_by_region:
        logging.info(f"clustering enriched pixels in region: {region}")
        # Using genomic corrdinated for clustering, not bin_id
        pixel_clust = clust_2D_pixels(
            _df,
            threshold_cluster=dots_clustering_radius,
            bin1_id_name="start1",
            bin2_id_name="start2",
        )
        pixel_clust_list.append(pixel_clust)
    logging.info("Clustering is complete")

    # concatenate clustering results ...
    # indexing information persists here ...
    if not pixel_clust_list:
        logging.warning("No clusters found for any regions! Output will be empty")
        empty_output = pd.DataFrame(
            [],
            columns=list(scored_df.columns)
            + [
                assigned_regions_name + "1",
                assigned_regions_name + "2",
                "c_label",
                "c_size",
                "cstart1",
                "cstart2",
            ],
        )
        return empty_output  # Empty dataframe with the same columns as anticipated
    else:
        pixel_clust_df = pd.concat(
            pixel_clust_list, ignore_index=False
        )  # Concatenate the clustering results for different regions

    # now merge pixel_clust_df and scored_df DataFrame ...
    # TODO make a more robust merge here
    df = pd.merge(
        scored_df, pixel_clust_df, how="left", left_index=True, right_index=True
    )
    # TODO check if next str-cast is neccessary
    df[assigned_regions_name + "1"] = df[assigned_regions_name].astype(str)
    df[assigned_regions_name + "2"] = df[assigned_regions_name].astype(str)
    # report only centroids with highest Observed based on overlapped libraries in each cluster c_label:
    centroid_list = []
    chrom_clust_group = df.groupby(
        [assigned_regions_name + "1", assigned_regions_name + "2", "c_label"],
        observed=True,
    )
    for chrom_clust, _df in chrom_clust_group:
        condition_names=np.unique('&'.join(list(_df.condition)).split("&"))
        _df['sum_balanced']=_df[[condition+'_balanced' for condition in condition_names]].copy().sum(axis=1)
        centroid = _df.loc[_df['sum_balanced'].idxmax()]
        centroid['sample_name']='&'.join(list(condition_names))
        centroid_list.append(centroid)
        
    centroids=pd.concat(centroid_list, ignore_index=False, axis=1).T
    centroids = centroids[[
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "sample_name",
    ]]
    
    return centroids, df

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



input_cooler_paths=pd.read_csv(sys.argv[1], sep='\t')
input_loop_paths=pd.read_csv(sys.argv[2], sep='\t')
resolution=int(sys.argv[3])
flank=int(sys.argv[4])
clr_weight_name=sys.argv[5]
# Available assembly names in bioframe package (e.g., hg38 and mm10)
# See details in https://bioframe.readthedocs.io/en/latest/guide-io.html#curated-genome-assembly-build-information
assembly_name=sys.argv[6]
dots_clustering_radius=int(sys.argv[7])
output_filename=sys.argv[8]
nproc=int(sys.argv[9])
output_clusters=f'clusters_of_external_loops.resolution.{int(resolution/1000)}kb.tsv'

# Merge loops called from each individual library (by exactly overlapping)
conditions = input_loop_paths["name"].tolist()
loop_files = dict(zip(input_loop_paths["name"], input_loop_paths["path"]))

loop_list=[]
for cond in conditions:
    loop_df = pd.read_csv(
        loop_files[cond], # bedpe format
        usecols=range(6),
        header=None,
        sep='\t',
        names=["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    )
    loop_df['sample_name']=cond
    loop_list.append(loop_df)
loops=pd.concat(loop_list).reset_index(drop=True)

# Join those conditions with & if these loops are exactly overlapped at both anchors.
merged_loops=loops.groupby(['chrom1','start1','end1','chrom2','start2','end2'])['sample_name'].apply(lambda x: '&'.join(list(np.unique(x)))).reset_index()

merged_loops=merged_loops.rename({'sample_name': 'condition'}, axis=1)

cooler_files = {}
for i in range(len(input_cooler_paths)):
    cooler_files[input_cooler_paths.name[i]] = input_cooler_paths.path[i]+f'::resolutions/{resolution}'

loops_annotated=annotate_balanced(merged_loops, cooler_files)

# Use bioframe to fetch the genomic features from the UCSC.
chromsizes = bioframe.fetch_chromsizes(assembly_name)
cens = bioframe.fetch_centromeres(assembly_name)
full_arms = bioframe.make_chromarms(chromsizes, cens)

# Select only chromosomes that are present in the cooler.
chrom_list = merged_loops["chrom1"].unique().tolist()
arms = full_arms.loc[full_arms.chrom.isin(chrom_list)].reset_index(drop=True)

loops_annotated = assign_regions(loops_annotated, arms)

centroids, clusters_of_loops = clustering_step_across_conditions(loops_annotated)

coord_cols = ["start1", "end1", "start2", "end2"]

centroids = centroids.copy()

for col in coord_cols:
    centroids[col] = pd.to_numeric(centroids[col], errors="raise").astype("int64")

dots = centroids.loc[:,['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]


loop_strengths=pd.DataFrame()
for cond in conditions:
    clr = cooler.Cooler(input_cooler_paths.loc[input_cooler_paths.name == cond,'path'].reset_index(drop=True)[0]+'::resolutions/'+str(resolution))
    expected = cooltools.expected_cis(
        clr,
        view_df=arms,
        clr_weight_name=clr_weight_name,
        nproc=nproc,)
    pileup_mtx=cooltools.pileup(clr, dots, view_df=arms, expected_df=expected, flank=flank, clr_weight_name=clr_weight_name, nproc=nproc)
    loop_strengths[cond] = quantify_individual_loops_in_stack(pileup_mtx)

loop_strengths_formatted=pd.concat([dots, loop_strengths], axis=1)
loop_strengths_formatted=pd.merge(centroids, loop_strengths_formatted, on=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'], how="inner")

loop_strengths_formatted.to_csv(output_filename, sep='\t', index=False, header=True)
clusters_of_loops.to_csv(output_clusters, sep='\t', index=False, header=True)



