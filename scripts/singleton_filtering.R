#!/usr/bin/env Rscript

library(dplyr)
library(stringr)

# compare each singleton with other nonsingleton clusters one by one
# Remove those singletons without shared anchors with other nonsingleton clusters (at least two pixels in that cluster)
.check_shared_singletons = function(singleton, nonsingletons){
    nonsingletons=nonsingletons[nonsingletons['chrom1'] == singleton['chrom1'],]
    share_left_anchor=any((singleton['start1'] >= nonsingletons['start1']) & (singleton['end1'] <= nonsingletons['end1']))
    if (share_left_anchor){
        return(share_left_anchor)
    }
    else{
        share_right_anchor=any((singleton['start2'] >= nonsingletons['start2']) & (singleton['end2'] <= nonsingletons['end2']))
        if (share_right_anchor){
            return(share_right_anchor)
        }
        else{
            return(FALSE)
        }
    }
}


remove_singletons_without_shared_anchors = function(df, resolution){
    df_singletons=df[(df['end1']-df['start1'] == resolution) & (df['end2']-df['start2'] == resolution),]
    df_nonsingletons=df[(df['end1']-df['start1'] != resolution) | (df['end2']-df['start2'] != resolution),]
    df_singletons=cbind(df_singletons,is_shared_anchor=apply(df_singletons,1,.check_shared_singletons,nonsingletons=df_nonsingletons))
    df_singletons_filtered=df_singletons[df_singletons[['is_shared_anchor']],c('chrom1','start1','end1','chrom2','start2','end2')]
    return(df_singletons_filtered)
}


# Remove those library-specific (HiC matrix) singletons that do not pass FDR_orphan_threshold (>0.02) according to Rao's HICCUPS paper.
remove_library_specific_singletons_hiccups=function(df, resolution, enriched_pixels_meta, FDR_orphan_threshold=0.02){
    df_library_specific_singletons=df[(df['end1']-df['start1'] == resolution) & (df['end2']-df['start2'] == resolution) & (df[['sample_name']] %in% enriched_pixels_meta$name),]
    
    df_library_specific_singletons_filtered=list()
    for (i in 1:length(enriched_pixels_meta)){
        enriched_pixels=read.table(file = enriched_pixels_meta$path[i], sep = '\t', header = TRUE)
        df=merge(df_library_specific_singletons[df_library_specific_singletons['sample_name'] == enriched_pixels_meta$name[i],], enriched_pixels, by=c('chrom1','start1','end1','chrom2','start2','end2'))
        df_library_specific_singletons_filtered[[i]]=df[rowSums(df[c("la_exp.lowleft.qval","la_exp.donut.qval","la_exp.vertical.qval","la_exp.horizontal.qval")]) <= FDR_orphan_threshold, c('chrom1','start1','end1','chrom2','start2','end2','sample_name')]
    }
    return(do.call(rbind, df_library_specific_singletons_filtered))
}




args <- commandArgs(trailingOnly = TRUE)
enriched_pixels_meta=read.table(file = 'enriched_pixels_meta.tsv', sep = '\t', header = TRUE)
resolution=as.numeric(args[1])
centroids=read.table(file = paste0('centroids_of_clusters_of_enriched_pixels.resolution.',resolution/1000,'kb.tsv'), sep = '\t', header = TRUE)
clusters_of_enriched_pixels=read.table(file = paste0('clusters_of_enriched_pixels.resolution.',resolution/1000,'kb.tsv'), sep = '\t', header = TRUE)
output='union_list_of_loops_tmp.tsv'


clusters=clusters_of_enriched_pixels %>%
           group_by(region, c_label) %>%
           summarise(start1 = min(start1), end1 = max(end1), start2 = min(start2), end2 = max(end2))
clusters[c('chrom1', 'arm')] <- str_split_fixed(clusters$region, '_', 2)
clusters['chrom2'] <- clusters['chrom1']
clusters <- clusters[c('chrom1','start1','end1','chrom2','start2','end2')]

singletons_with_shared_anchors=remove_singletons_without_shared_anchors(clusters, resolution)
singletons_with_shared_anchors <- merge(singletons_with_shared_anchors, clusters_of_enriched_pixels, by=c('chrom1','start1','end1','chrom2','start2','end2'))[c('chrom1','start1','end1','chrom2','start2','end2','condition')]
colnames(singletons_with_shared_anchors)[7] <- 'sample_name'

singletons_hiccups_filtered=remove_library_specific_singletons_hiccups(singletons_with_shared_anchors, resolution, enriched_pixels_meta, FDR_orphan_threshold=0.02)

df_singletons=clusters[(clusters['end1']-clusters['start1'] == resolution) & (clusters['end2']-clusters['start2'] == resolution),c('chrom1','start1','end1','chrom2','start2','end2')]
# get all centroids of nonsingleton clusters
centroids_of_nonsingleton_clusters=anti_join(centroids, df_singletons, by=c('chrom1','start1','end1','chrom2','start2','end2'))[,c('chrom1','start1','end1','chrom2','start2','end2','sample_name')]

# combine the rest of singletons passed filering and centriods of nonsingleton clusters
dots=rbind(singletons_hiccups_filtered, centroids_of_nonsingleton_clusters)

write.table(dots, output, quote = FALSE, row.names = FALSE, sep = "\t")
