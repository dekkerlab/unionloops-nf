# unionloops-nf

**Nextflow Version:** 22.10.6

### A HiCCUPS-based pipeline for simultaneous chromatin loop calling and cross-sample comparison, where both steps reinforce each other.

The `unionloops` pipeline provides:

- Enhanced sensitivity of loop detection using cross-sample evidence.
- Improved loop positional precision relative to CTCF/RAD21 sites.
- Loop annotation via clustering: shared vs. sample-specific.
- Loop strength quantification across all samples.

---

## Setup

### Step 1: Check Conda version and solver
#### 1. Show the configured solver
```bash
# Check Conda version
conda --version

# Show current solver configuration ("classic" or "libmamba")
conda config --show solver

# Note: On older Conda versions, if nothing is printed, it defaults to "classic".
# Make sure your current solver configuration is "libmamba", which is much faster than "classic".
```

#### 2. Update Conda to the latest version
```bash
# Example: Update to Conda v25.7.0 and force reinstall in the base environment
conda install -n base -c defaults conda=25.7.0 --force-reinstall
```

#### 3. Double-check versions and solver
```bash
# Verify Conda version after update
conda --version

# Check solver again (default on modern Conda is "libmamba")
conda config --show solver
```

### Step 2: Create conda environment for Nextflow

```bash
conda env create -f nextflow_env.yml
```

### Step 3: Create conda environment for `unionloops` pipeline

```bash
conda env create -f unionloops_env.yml
```

---

## Input configuration

Prepare a TSV file (e.g., `mcool_paths.tsv`) with sample names and `.mcool` file paths:

```tsv
name	path
sample1	/full/path/to/sample1.mcool
sample2	/full/path/to/sample2.mcool
sample3	/full/path/to/sample3.mcool
MEGA	/full/path/to/MEGA.mcool    # Optional: merged high-resolution map
```

For optimal performance, also consider including a single merged `MEGA.mcool` file that combines all samples. The high signal-to-noise ratio of the MEGA map can:

- Improve the rescue of **sample-specific** loops.
- Enhance the **precision** of loop detection.

> Use `.mcool` files generated from **distiller v0.3.3** ([distiller-nf](https://github.com/open2c/distiller-nf)) for best compatibility.

---

## Running the pipeline

You can launch `unionloops` using different hardware profiles:

1. Default hardware profile (`configs/local.config`) with your `mcool_paths.tsv` and conda env `unionloops-nf`:

```bash
nextflow run /full/path/to/unionloops-nf/unionloops.nf \
    -ansi-log false \
    --input_cooler_paths /full/path/to/mcool_paths.tsv \
    --outfilename union_loop_list.tsv \
    --conda_env /full/path/to/miniconda3/envs/unionloops-nf
```

2. `cluster` hardware profile (`configs/cluster.config`) with your `mcool_paths.tsv` and conda env `unionloops-nf`:

```bash
nextflow run /full/path/to/unionloops-nf/unionloops.nf \
    -profile cluster \
    -ansi-log false \
    --input_cooler_paths /full/path/to/mcool_paths.tsv \
    --outfilename union_loop_list.tsv \
    --conda_env /full/path/to/miniconda3/envs/unionloops-nf
```

3. `custom` hardware profile with your own configuration file with your `mcool_paths.tsv` and conda env `unionloops-nf`:

```bash
nextflow run /full/path/to/unionloops-nf/unionloops.nf \
    -profile custom --custom_config /full/path/to/your.config \
    -ansi-log false \
    --input_cooler_paths /full/path/to/mcool_paths.tsv \
    --outfilename union_loop_list.tsv \
    --conda_env /full/path/to/miniconda3/envs/unionloops-nf
```

You may override default parameters defined in `nextflow.config` as needed (see parameters section).

---

## Example output

By default, output files will be saved in the `results/` directory relative to your working directory.

### File structure

```
results/
├── enriched_pixels/               # Enriched pixels per sample
│   ├── sample1.enriched.pixels.resolution.10kb.tsv
│   ├── sample2.enriched.pixels.resolution.10kb.tsv
│   ├── sample3.enriched.pixels.resolution.10kb.tsv
│   └── MEGA.enriched.pixels.resolution.10kb.tsv
│
├── clusters/                      # Clustering results of pooled enriched pixels across all samples
│   ├── centroids_of_clusters_of_enriched_pixels.resolution.10kb.tsv # Without additional filtering
│   ├── clusters_of_enriched_pixels.resolution.10kb.tsv
│   └── enriched_pixels_meta.tsv
│
└── union_loop_list_10kb.tsv       # Final union list of loops
```

### Columns in the dataframe of the final union list of loops

| Column       | Description                                     |
| ------------ | ----------------------------------------------- |
| chr1         | Chromosome of anchor 1                          |
| start1       | Start position of anchor 1                      |
| end1         | End position of anchor 1                        |
| chr2         | Chromosome of anchor 2                          |
| start2       | Start position of anchor 2                      |
| end2         | End position of anchor 2                        |
| sample\_name | Detected sample(s); joined with `&` if multiple |
| sample1      | Loop strength in sample1                        |
| sample2      | Loop strength in sample2                        |
| sample3      | Loop strength in sample3                        |
| MEGA         | Loop strength in MEGA                           |

---

## Parameters

### Required parameters

| Parameter            | Description                              |
| -------------------- | ---------------------------------------- |
| input\_cooler\_paths | TSV with sample names and `.mcool` paths |
| outfilename          | Output filename for union loop list      |
| conda\_env           | Path to conda environment for pipeline   |

### Optional parameters

| Parameter                | Default       | Description                                                       |
| ------------------------ | ------------- | ----------------------------------------------------------------- |
| assembly\_name           | hg38          | Genome assembly name from UCSC database                           |
| resolution               | 10000         | Resolution for loop detection (bp)                                |
| outdir                   | results       | Output directory                                                  |
| custom\_config           | custom.config | Custom Nextflow config                                            |
| clr\_weight\_name        | weight        | Used by cooltools functions                                       |
| max\_loci\_separation    | 10000000      | Maximum loci separation for loop-calling (bp)                     |
| max\_nans\_tolerated     | 1             | Used in `cooltools.dots()`                                        |
| lambda\_bin\_fdr         | 0.1           | Used in `cooltools.dots()`                                        |
| tile\_size               | 5000000       | Used in `cooltools.dots()`                                        |
| nproc                    | 1             | Number of processes used                                          |
| dots\_clustering\_radius | 20000         | Clustering radius for enriched pixels (bp)                        |
| flank                    | 100000        | Flanking region for strength estimation (typically 10×resolution) |

---

## Test example

### Test data details

1. **HFF_MicroC**  
   - Description: Micro-C data from HFF human cells for two chromosomes (hg38) in a multi-resolution mcool format.  
   - Source: Krietenstein et al. 2021  
   - Downloaded from: https://osf.io/3h9js/download  
   - Stored as: `test.mcool`  
   - Original MD5 checksum: `e4a0fc25c8dc3d38e9065fd74c565dd1`

2. **hESC_MicroC**  
   - Description: Micro-C data from human ES cells for two chromosomes (hg38) in a multi-resolution mcool format.  
   - Source: Krietenstein et al. 2021  
   - Downloaded from: https://osf.io/3kdyj/download  
   - Stored as: `test_hESC.mcool`  
   - Original MD5 checksum: `ac0e636605505fb76fac25fa08784d5b`

### Run the pipeline using test data

#### Step 1: Clone this repository and change to the test directory
```bash
$ git clone https://github.com/dekkerlab/unionloops-nf.git
$ cd unionloops-nf/test/
```

#### Step 2: Download two test `.mcool` files to `test/data/` and generate a `test_mcool_paths.tsv` file under `test/`
```bash
$ bash ./run_download.sh
```

#### Step 3: Run the pipeline with the downloaded test data
```bash
$ nextflow run ../unionloops.nf \
>  -ansi-log false \
>  --input_cooler_paths /full/path/to/test/test_mcool_paths.tsv \
>  --outfilename test_union_loop_list.tsv \
>  --conda_env ~/miniconda3/envs/unionloops-nf
```
**Note:** You might need to replace `~/miniconda3/envs/unionloops-nf` with the path to your `unionloops-nf` conda environment. You can find it by running:
```bash
$ conda env list | grep 'unionloops-nf'
```

#### Step 4: Take a look at the final union list of loops
```bash
$ head results/test_union_loop_list.tsv
```

---

## Citations

- **distiller**: A modular Hi-C mapping pipeline for reproducible data analysis. [https://github.com/open2c/distiller-nf](https://github.com/open2c/distiller-nf)

- **HiCCUPS**: Rao et al. *A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping*. Cell, 2014. [https://doi.org/10.1016/j.cell.2014.11.021](https://doi.org/10.1016/j.cell.2014.11.021)

- **Nextflow**: Tommaso et al. *Nextflow enables reproducible computational workflows*. Nature Biotechnology, 2017. [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

- **cooltools**: Open2C. *Cooltools: scalable analysis tools for Hi-C and other genome-wide contact maps*. PLOS Computational Biology, 2024. [https://doi.org/10.1371/journal.pcbi.1012067](https://doi.org/10.1371/journal.pcbi.1012067)

- **cooler**: Abdennur and Mirny. *Cooler: scalable storage for Hi-C data and other genomically labeled arrays*. Bioinformatics, 2020. [https://doi.org/10.1093/bioinformatics/btz540](https://doi.org/10.1093/bioinformatics/btz540)

- **bioframe**: Open2C et al. Bioframe: Operations on Genomic Intervals in Pandas Dataframes. Bioinformatics, 2024. [https://doi.org/10.1093/bioinformatics/btae088](https://doi.org/10.1093/bioinformatics/btae088)

---

## Methods paper (To Do)

Title

