#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import bioframe
import cooler
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




input_cooler_paths=pd.read_csv(sys.argv[1], sep='\t')
resolution=int(sys.argv[2])
# Available assembly names in bioframe package (e.g., hg38 and mm10)
# See details in https://bioframe.readthedocs.io/en/latest/guide-io.html#curated-genome-assembly-build-information
assembly_name=sys.argv[3]
dots_clustering_radius=int(sys.argv[4])
output_centroids=f'centroids_of_clusters_of_enriched_pixels.resolution.{int(resolution/1000)}kb.tsv'
output_clusters=f'clusters_of_enriched_pixels.resolution.{int(resolution/1000)}kb.tsv'

enriched_pixels_filenames=[]
for i in range(len(input_cooler_paths)):
    filename=input_cooler_paths.name[i] + f'.enriched.pixels.resolution.{int(resolution/1000)}kb.tsv'
    enriched_pixels_filenames.append(filename)

meta_enriched_pixels_files = {
    'name' : list(input_cooler_paths.name),
    'path' : enriched_pixels_filenames
}
enrich_pixels_paths = pd.DataFrame(meta_enriched_pixels_files)
# save the paths of enriched pixels in enriched_pixels_meta.tsv as the input for singleton_filtering.R
enrich_pixels_paths.to_csv('enriched_pixels_meta.tsv', sep='\t', index=False, header=True)

# Use bioframe to fetch the genomic features from the UCSC.
chromsizes = bioframe.fetch_chromsizes(assembly_name)
cens = bioframe.fetch_centromeres(assembly_name)
full_arms = bioframe.make_chromarms(chromsizes, cens)

cooler_files = {}
arms_list=[]
for i in range(len(input_cooler_paths)):
    cooler_files[input_cooler_paths.name[i]] = input_cooler_paths.path[i]+f'::resolutions/{resolution}'
    clr=cooler.Cooler(input_cooler_paths.path[i]+f'::resolutions/{resolution}')
    # Select only chromosomes that are present in the cooler.
    arms_list.append(full_arms.set_index("chrom").loc[clr.chromnames].reset_index())

# Select arms shared in all samples.
arms_df=pd.concat(arms_list, ignore_index=True)
## Step 1: Count occurrences of each row
row_counts = arms_df.value_counts().reset_index(name='count')
## Step 2: Filter rows that appear exactly len(input_cooler_paths)
rows_len_times = row_counts[row_counts['count'] == len(input_cooler_paths)].drop(columns='count')
## Step 3: Merge back to original to select matching rows
df_len_times = arms_df.merge(rows_len_times, on=list(arms_df.columns))
## Step 4: Get unique rows among those
intersection_arms = df_len_times.drop_duplicates()
## Step 5: Sorted the intersection arms based on the order of full arms (left)
arms=pd.merge(full_arms, intersection_arms)

# Concatennate enriched pixels from all samples
list_enriched_pixels=[]
for i in range(len(enrich_pixels_paths)):
    df=pd.read_csv(enrich_pixels_paths.path[i], sep='\t')
    # only select loops whose chromosomes are shared in all samples
    df=df[df['chrom1'].isin(arms['chrom'])]
    df['condition']=enrich_pixels_paths.name[i]
    list_enriched_pixels.append(df)

df_enriched_pixels=pd.concat(list_enriched_pixels, ignore_index=True)
df_enriched_pixels=df_enriched_pixels[['chrom1','start1','end1','chrom2','start2','end2','condition']].copy().groupby(['chrom1','start1','end1','chrom2','start2','end2'])['condition'].apply(lambda x: '&'.join(list(np.unique(x)))).reset_index()

filtered_pixels_annotated=annotate_balanced(df_enriched_pixels, cooler_files)
filtered_pixels_annotated = assign_regions(filtered_pixels_annotated, arms)

centroids, clusters_of_enriched_pixels = clustering_step_across_conditions(filtered_pixels_annotated, dots_clustering_radius=dots_clustering_radius)

centroids.to_csv(output_centroids, sep='\t', index=False, header=True)
clusters_of_enriched_pixels.to_csv(output_clusters, sep='\t', index=False, header=True)
