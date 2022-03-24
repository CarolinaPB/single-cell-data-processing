# Single-cell preprocessing pipeline

## First follow the instructions here

[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT

This pipeline includes the first steps in the analysis of Single-cell data.  
The first step is getting the reference package for your species. This will be used for read alignment and gene expression quantification. If you're working with human or mouse, you can download the reference from the Cellranger website, if not, the pipeline can create the reference for you. (details below)

Once you have the reference package, the pipeline starts by running [Cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count). Cellranger count performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis.    
If the fastq files are not named in the format accepted by Cellranger count: `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`, you can specify in the config file that these need to be renamed: option `RENAME: y` or you can rename them yourself to follow this naming convention.

The metrics from Cellranger count for all samples are combined into one file `cellranger_count_metrics_allsamples.tsv`. This will have information such as "estimated number of cells", and "mean reads per cell".

After the Cellranger count step, it's important to remove the ambient RNA, which is RNA that has been released from degraded or dying cells and is now in the cell suspension. The R package [SoupX](https://github.com/constantAmateur/SoupX) is used to correct for ambient RNA. In addition to the output files with the corrected data, one html document is created per sample processed (`2_ambient_RNA_correction/Ambient_RNA_correction_<sample>.html`). This html file shows the code used to perform the ambient RNA correction, as well as a few plots that illustrate this process - for the 5 most affected genes and for 5 random genes:

- Plot 1: in which cells the gene is expressed
- Plot 2: ratio of observed to expected counts
- Plot 3: change in expression due to correction

Once the data has been corrected for ambient RNA, it's time for quality control filtering. This is a step that depends on the cell type, library preparation method used, etc, so you should always check if the default parameters make sense, use your own, or even run several times with different ones.

QC is run for every sample separately. First [Scanpy](https://scanpy.readthedocs.io/en/stable/) calculates some general QC metrics for genes and cells. It will also calculate the proportion of counts for mitochondrial genes. Several plots will be created to help assess the quality of the data:  
Before filtering:

- Violin plots showing:
  - n_genes_by_counts: number of genes with positive counts in a cell
  - total_counts: total number of counts for a cell
  - pct_counts_mt: proportion of mitochondrial counts for a cell
- Scatter plot showing :
  - total_counts vs pct_counts_mt
  - total counts vs n_genes_by_counts  

After filtering:

- Percentage of counts per gene for the top 20 genes after filtering
- Violin plots showing:
  - n_genes_by_counts: number of genes with positive counts in a cell
  - total_counts: total number of counts for a cell
  - pct_counts_mt: proportion of mitochondrial counts for a cell

The final preprocessing step is doublet removal with [Scrublet](https://github.com/swolock/scrublet). This step may be run more than once to determine the ideal doublet score threshold. The histogram shown in `4_Doublets/<sample>/histogram_<sample>_doublets.pdf` should show a biomodal distribution, and the threshold shown in the "simulated doublets" plot should be at the minium between the two modes. The first run should be with parameter `SCRUB_TRESHOLD:`. Once the first run has finished, you should look at the histograms of all samples and see if you need to change the treshold. If you do, set it in the config file as explained further below and run the `remove_doublets` step again.

There will be an `h5ad` object containing the preprocessed data for each sample - ambient RNA removed, QC filtered and doublets removed - in the `4_Doublets` directory.  
This file can be used for further analysis.

#### Tools used

- Cellranger:
  - [mkgtf](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf) - filter GTF. [default: off]
  - [mkref](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkref) - create reference
  - [count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) - create feature counts
- R 
  - combine Cellranger count sample metrics
  - [SoupX](https://github.com/constantAmateur/SoupX) - remove ambient RNA
- Python
  - [Scanpy](https://scanpy.readthedocs.io/en/stable/index.html) - QC filtering
  - [Scrublet](https://github.com/swolock/scrublet) - Doublet removal

| ![DAG](https://github.com/CarolinaPB/single-cell-data-processing/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |

### Edit config.yaml with the paths to your files and set parameters

```yaml
DATA: /path/to/data/dir
OUTDIR: /path/to/outdir

# mkref options
MKREF: <y/n>
FASTA: /path/to/fasta
GTF: /path/to/gtf # if creating own reference
REF_VERSION: 
  - "--ref-version=<version>"
CR_MKREF_EXTRA: ""

# Filter GTF
FILTER_GTF: y
# # see here for available biotypes https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf
ATTRIBUTES:
  - "--attribute=<biotype>"

PREFIX: <species>

# rename fastq files
RENAME: <y/n>

# Cell ranger count options 
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count
CR_COUNT_extra: ""

# QC parameters
MITO_PERCENTAGE: 10 # keep cells with less than X% mitochondrial read fraction
NUMBER_GENES_PER_CELL: 500 # keep cells with more than X genes
NUMBER_UMI_PER_CELL: 1000 # keep cells with more than X UMIs
ENSEMBLE_BIOMART_SPECIES: "<species>" # ensembl biomart species used to get the mitochondrial genes for that species

# threshold doublet score (should be at the minimum between two modes of the simulated doublet histogram)
SCRUB_THRESHOLD: 
  <sample1>: <value>
  <sample2>: <empty>
```

- **DATA** - path to directory containing fastq files. Preferrably, the files should be named in the format accepted by Cellranger Count `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`. If they are, set `RENAME: n`. If not, they should be in the format `<sample>_R1.fastq.gz`. In this case, you should set `RENAME: y` so that the pipeline will rename the files according to the necessary format for Cellranger Count.
The path can be to a directory that contains a subdirectory per sample. For example:

    ```text
    DATA
    ├── SAMPLE_1
    │   ├──  sample_1_R1.fastq.gz
    │   └──  sample_1_R2.fastq.gz
    └── SAMPLE_2
       ├──  sample_2_R1.fastq.gz
       └──  sample_2_R2.fastq.gz
    ```

    Or to a directory with fastqs for all samples:

    ```text
    DATA
    ├── sample_1_R1.fastq.gz
    ├── sample_1_R2.fastq.gz
    ├── sample_2_R1.fastq.gz
    └── sample_2_R2.fastq.gz
    ```

- **OUTDIR** - directory where snakemake will run and where the results will be written to.  
If you don't want the results to be written to a new directory, open config.yaml and comment out `OUTDIR: /path/to/outdir`
- **MKREF** `y` if Cell Ranger doesn't provide a reference package for your species (currently, the species available are Human and Mouse). `n` if you're using an existing reference package. In this case, you should create a directory in the pipeline directory named `<prefix>_genome`, and download the reference from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).  

---  
**You only need to set the following if `MKREF: y`:**

- **FASTA** - path to fasta file
- **GTF**: path to gft file
- **PREFIX**: name of your species. Used to name the directory containing the reference package
- **REF_VERSION**: Reference version string to include with reference
- **CR_MKREF_EXTRA**: any other options for cellranger mkref. [see here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#singl)

----
- **PREFIX** - The name of your organism. The reference package used for cellranger count will be in the `<prefix>_genome` directory
- **RENAME** - `y` if your input fastqs are not named in this format `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`. Use `n` if they are.
- Options for Cellrange count:
  - **CR_COUNT_extra** - any other options for cellranger count. [Find other options here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count). [Default: ""]
- QC parameters
  - **MITO_PERCENTAGE** - Keep cells with less than X% mitochondrial read fraction. [Default: 10]
  - **NUMBER_GENES_PER_CELL** - keep cells with more than X genes. [Default: 500]
  - **NUMBER_UMI_PER_CELL** - keep cells with more than X UMIs. [Default: 1000]
  - **ENSEMBLE_BIOMART_SPECIES** -  ensembl biomart species used to get the mitochondrial genes
- Doublet removal score
  - **SCRUB_THRESHOLD** - threshold doublet score. It should be at the minimum between two modes of the simulated doublet histogram.
  In the first run it should be run as `SCRUB_TRESHOLD:` (with no parameters).
  After that is done, for each sample you should then look at the `4_Doublets/<sample>/histogram_<sample>_doublets.pdf` plot and see if the vertical line on the "simulated doublets" plot is at the minimum between the two modes. If it's not, you should manually set it in the config file as:

    ```yaml
    SCRUB_THRESHOLD: 
      <sample 1>: <value>
      <sample 2>: <empty>
    ```

    In `SCRUB_THRESHOLD` there should be a line for each sample, even if you don't need to set the threshold for that sample. If you need to change the treshold, set `<sample>: <value>`, if not, set `<sample>:` (without value).

## Additional set up

### Install Cellranger
Follow the instructions [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

1. First download the package from the [downloads page](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). For example:
```sh
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1648081234&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDgwODEyMzR9fX1dfQ__&Signature=HqCwx6eBEj~Lyw7C7UvsMAHzUH9aiPSM5yFcyflZiL2JRIwqzY2VWz1COtDQHNoJ48Ve41LZ5Q3eGv1yaAEf88SGhtxRUb2wJhFvvixBoR550bQ2wK7qfL6buLL9~u7MPw4q0-c1adXaSCm6otd6Xn0x2FIpZimOGJMYI9QEvNStN1Hi6MH4ZUOHGFFRBAvxlRxHmYBk-Vr~6qdc7nFXJW0C8OBWTn2g~XSKZRD50B5G5StMis0lLmgXZbRS0htQu8LPuUp8ZxqxQv20m9-HV9jEDVYEUP1sNJzAHGhAtq1FajN572Lptq0cWES8fheMexht1l-wRbQA-yOKAp7Bzg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

Attention: the `cellranger-6.1.2.tar.gz` name is an example, when you use this pipeline there might be a more recent version available.

2. Unpack Cellranger:
```sh
tar -xzvf cellranger-6.1.2.tar.gz
```
This will create a new directory, `cellranger-6.1.2`, that contains cellranger and its dependencies. 

3. Place the path to the `cellranger-6.1.2` (or the version you installed) in the config.yaml file. It should look like this

```yaml
CELLRANGER_PATH: /path/to/cellranger-6.1.2
```
> don't add the backslash ("\\") after the directory name

### Reference package
If you're working with human or mouse data, download the reference from [here](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) and place it in a folder in the pipeline directory called `<prefix>_genome`.  
If you're working with another organism, download the fasta file and gtf file for your organism and place them in a directory called `<prefix>_genome` directory (should be in the pipeline main directory - where the Snakefile is). You can download these from [Ensembl](https://www.ensembl.org/index.html).

### How to run

1. Dry run to check if everything was correctly set up and if the pipeline is ready to run

    ```sh
    snakemake -np
    ```

2. If all looks good, run the pipeline with

      ```sh
      snakemake --profile <name of hpc profile>
      ```

3. Once you have the `remove_doublets` results for each sample, you should look at the histogram, either in the jupyter notebook `4_Doublets/processed_notebook_<sample>.ipynb` or in the saved plot `4_Doublets/<sample>/histogram_<sample>_doublets.pdf`. The vertical bar (threshold) in the "simulated doublets" plot should be at the lowest point between the two modes. If not, you'll need to set it in the config file. Only change the threshold for the samples that need it. For some, the automatically set threshold may be good. To do this, edit the config file. Add your sample name (the same as the directory names in `4_Doublets`) and the new threshold or, in case you don't need to add a threshold, add the sample name and nothing in front of it.

    ```yaml
    SCRUB_THRESHOLD: 
      <sample 1>: <value>
      <sample 2>:
      <sample 3>: <value>
    ```

    In this case, the scrublet step for `sample 1` and `sample 3` will run with the user defined threshold and `sample 2` will run with the automatically defined treshold.

4. After changing the thresholds, you'll need to run the `remove_doublets` step again. Before you do this, you need to copy the previous results (`4_Doublets`) to another directory, or you need to delete those results. Example to move the results to another directory:

    ```sh
    mkdir -p 4_Doublets_first_run
    mv -r 4_Doublets/* 4_Doublets_first_run
    ```

5. Once that's done you can rerun the `remove_doublets` step with

    ```sh
    snakemake -np --forcerun remove_doublets
    snakemake --profile <profile name> --forcerun remove_doublets
    ```

## RESULTS

You will have results for each step of the pipeline.

- 1_renamed: directory with renamed fastq files softlinked
- 2_ambient_RNA_correction: directory containing results from ambient RNA correction
  - Ambient_RNA_correction_<sample>.html: shows the code used for the ambient RNA correction, as well as a few plots that illustrate this process - for the 5 most affected genes and for 5 random genes:

    Plot 1: in which cells the gene is expressed  
    Plot 2: ratio of observed to expected counts  
    Plot 3: change in expression due to correction  

- 2_ambient_RNA_correction_data: inside, the ambient RNA corrected data for each sample is its corresponding directory.

- 3_QC: directory containing results from QC:
  - processed_notebook_<sample>.ipynp - Jupyter notebooks used to calculate QC for each sample. These are interactive and can be used to do further QC.
  - <sample>.h5ad - Filtered data.
  - <sample> - Directory with QC plots

      ##### Description of the QC plots

      Before filtering:

      - Violin plots showing:
        - n_genes_by_counts: number of genes with positive counts in a cell
        - total_counts: total number of counts for a cell
        - pct_counts_mt: proportion of mitochondrial counts for a cell
      - Scatter plot showing :
        - total_counts vs pct_counts_mt
        - total counts vs n_genes_by_counts  

      After filtering:

      - Percentage of counts per gene for the top 20 genes after filtering
      - Violin plots showing:
        - n_genes_by_counts: number of genes with positive counts in a cell
        - total_counts: total number of counts for a cell
        - pct_counts_mt: proportion of mitochondrial counts for a cell
      These jupyter notebooks are interactive and can be used to do further QC control.
- 4_Doublets: directory containing results from doublet removal
  - processed_notebook_<sample>.ipynb - Jupyter notebooks used to do the doublet removal step for each sample. These are interactive and can be used to test different scrublet tresholds
  - <sample> - directory with saved plots
  - <sample>_doublets.h5ad - filtered data
- \<sample> - directory containing output from Cellranger count. See [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/overview) for more information
  - outs
    - web_summary.html - summary metrics and automated secondary analysis results. If an issue was detected during the pipeline run, an alert appears on this page.
    - metrics_summary.csv - metrics like "estimated number of cells"
    - possorted_genome_bam.bam - BAM file containing position-sorted reads aligned to the genome and transcriptome, as well as unaligned reads. Each read in this BAM file has Chromium cellular and molecular barcode information attached.
    - raw_feature_bc_matrix.h5 - Contains every barcode from the fixed list of known-good barcode sequences that has at least one read. This includes background and cell associated barcodes
    - filtered_feature_bc_matrix.h5 - Contains only detected cell-associated barcodes. For Targeted Gene Expression samples, non-targeted genes are removed from the filtered matrix.
    - analysis - directory containing secondary analysis results: clustering, differential expression analysis, PCA, t-SNE, UMAP
    - molecule_info.h5 - contains per-molecule information for all molecules that contain a valid barcode, a valid UMI, and were assigned with high confidence to a gene or Feature Barcode.
    - cloupe.cloup - file to be used with [Loupe Browser](https://support.10xgenomics.com/single-cell-gene-expression/software/visualization/latest/what-is-loupe-cell-browser)


## Common issues

### Only see `create_file_log` and `combine_cellranger_counter_metrics` when running `snakemake -np`

This can be because your fastq files are not named correctly and the `rename` option is not set to `y`.
