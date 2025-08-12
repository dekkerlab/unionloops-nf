#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

// === Check for required parameters BEFORE using them ===
def required_params = [
    'input_cooler_paths',
    'outfilename',
    'conda_env'
]

required_params.each { p ->
    if (!params."$p") {
        error "Missing required parameter: --${p}"
    }
}

// Check input cooler resolution
if ((params.resolution as int) < 4000) {
    error "Provided cooler has resolution ${params.resolution} bases. Current version only supports >= 4000."
}


// === Now it's safe to use them ===
log.info """\
    U N I O N L O O P S - N F   P I P E L I N E
    ===========================================
    active_profile         : ${workflow.profile}
    cooler_paths           : ${params.input_cooler_paths}
    assembly               : ${params.assembly_name}
    resolution             : ${params.resolution}
    flank                  : ${params.flank}
    dots_clustering_radius : ${params.dots_clustering_radius}
    nf_path                : ${projectDir}
    conda_env              : ${params.conda_env}
    outdir                 : ${params.outdir}
    outfilename            : ${params.outfilename}
    clr_weight_name        : ${params.clr_weight_name}
    max_loci_separation    : ${params.max_loci_separation}
    max_nans_tolerated     : ${params.max_nans_tolerated}
    lambda_bin_fdr         : ${params.lambda_bin_fdr}
    tile_size              : ${params.tile_size}
    nproc                  : ${params.nproc}
    """
    .stripIndent().trim()


process generate_enriched_pixels {
    tag "$name"
    publishDir "${params.outdir}/enriched_pixels", mode:'copy'
    
    input:
    val name
    path cooler_path
    val assembly_name
    val resolution
    val clr_weight_name
    val max_loci_separation
    val max_nans_tolerated
    val lambda_bin_fdr
    val tile_size
    val nproc

    output:
    path '*'

    script:
    """
    python ${projectDir}/scripts/generate_enriched_pixels.py $cooler_path $assembly_name $resolution $clr_weight_name $max_loci_separation $max_nans_tolerated $lambda_bin_fdr $tile_size $nproc $name
    """
}

process cluster_enriched_pixels {
    publishDir "${params.outdir}/clusters", mode:'copy'
    
    input:
    path input_cooler_paths
    path input_cooler_files
    path enriched_pixels_files
    val assembly_name
    val resolution
    val dots_clustering_radius

    output:
    path '*'

    script:
    """
    python ${projectDir}/scripts/cluster_enriched_pixels.py $input_cooler_paths $resolution $assembly_name $dots_clustering_radius
    """
}

process singleton_filtering {

    input:
    path enriched_pixels_files
    path cluster_files
    val resolution

    output:
    path 'union_list_of_loops_tmp.tsv'

    script:
    """
    Rscript ${projectDir}/scripts/singleton_filtering.R $resolution
    """
}

process calculate_loop_strength {
    tag "$name"
    
    input:
    val name
    path input_cooler_paths
    path union_list_of_loops_file
    path input_cooler_file
    val resolution
    val flank
    val assembly_name
    val clr_weight_name
    val nproc

    output:
    path '*'

    script:
    """
    python ${projectDir}/scripts/calculate_loop_strength.py $input_cooler_paths $input_cooler_file $name $union_list_of_loops_file $resolution $flank $assembly_name $clr_weight_name $nproc
    """
}

process combine_loop_strength {
    publishDir params.outdir, mode:'copy'
    
    input:
    path input_cooler_paths
    path loop_strength_files
    val resolution
    val outfilename

    output:
    path "$outfilename"

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import os.path
    resolution=int("$resolution")
    cooler_paths=pd.read_csv("$input_cooler_paths", sep='\t')
    df=pd.read_csv(os.path.splitext(os.path.basename(cooler_paths.path[0]))[0] + f'.loop.strength.resolution.{int(resolution/1000)}kb.tsv', sep='\t').loc[:,['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'sample_name']]
    loop_strength=pd.DataFrame({cooler_paths.name[i]: pd.read_csv(os.path.splitext(os.path.basename(cooler_paths.path[i]))[0] + f'.loop.strength.resolution.{int(resolution/1000)}kb.tsv', sep='\t')['loop_strength'] for i in range(len(cooler_paths))})
    loop_strength=pd.concat([df, loop_strength], axis=1)
    loop_strength.to_csv("$outfilename", sep='\t', index=False, header=True)
    """
}

referenceFile = new FileReader(params.input_cooler_paths)
coolers = referenceFile.collect { it.split(/\t/) }.inject([:]) { map, val -> map[val[0]] = val[1]; map }
coolers.remove('name')
def name = coolers.collect{entry -> entry.key}
def path = coolers.collect{entry -> entry.value}

workflow {
    name_ch = Channel.fromList(name)
    path_ch = Channel.fromList(path)
    enriched_pixels_ch = generate_enriched_pixels(name_ch, path_ch, params.assembly_name, params.resolution, params.clr_weight_name, params.max_loci_separation, params.max_nans_tolerated, params.lambda_bin_fdr, params.tile_size, params.nproc)
    clusters_ch = cluster_enriched_pixels(params.input_cooler_paths, path_ch.collect(), enriched_pixels_ch.collect(), params.assembly_name, params.resolution, params.dots_clustering_radius)
    union_list_ch = singleton_filtering(enriched_pixels_ch.collect(), clusters_ch, params.resolution)
    loop_strength_ch=calculate_loop_strength(name_ch, params.input_cooler_paths, union_list_ch, path_ch, params.resolution, params.flank, params.assembly_name, params.clr_weight_name, params.nproc)
    combine_loop_strength(params.input_cooler_paths, loop_strength_ch.collect(), params.resolution, params.outfilename)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the final tsv file of union loops --> $params.outdir/$params.outfilename\n" : "Oops .. something went wrong" )
}


