#BSUB -q short
#BSUB -W 4:00
#BSUB -n 2
#BSUB -J unionloops-nf
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -eo dis.err
#BSUB -oo dis.out

# Load environment
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
elif [ -f "$HOME/.bash_profile" ]; then
    source "$HOME/.bash_profile"
fi

# Activate nextflow conda environment
conda activate nextflow

# Run Nextflow
nextflow run /full/path/to/unionloops-nf/unionloops.nf \
        -profile cluster \
        -ansi-log false \
        --input_cooler_paths /full/path/to/unionloops-nf/test/test_mcool_paths.tsv \
        --outfilename test_union_loop_list_10kb.tsv \
        --conda_env ~/miniconda3/envs/unionloops-nf \
        --nproc 2