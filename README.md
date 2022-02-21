# Single-cell preprocessing pipeline

## First follow the instructions here

[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT

This pipeline includes the first steps in the analysis of Single-cell data.
It starts by running [Cellranger count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count). Cellranger count performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis.  
If the fastq files are not named in the format accepted by Cellranger count: `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`, you can specify in the config file that these need to be renamed: option `RENAME: y` or you can rename them yourself to follow this naming convention.

The metrics from Cellranger count for all samples are combined into one file `cellranger_count_metrics_allsamples.tsv`. This will have information such as "estimated number of cells", and "mean reads per cell".

After the Cellranger count step, it's important to remove the ambient RNA. This is done through a custom R notebook. The R package [SoupX](https://github.com/constantAmateur/SoupX) is used to correct for ambient RNA. In addition to the output files with the corrected data, one html document is created per sample processed (`ambient_RNA_correction/Ambient_RNA_correction_<sample>.html`). This html file shows the code used to perform the ambient RNA correction, as well as a few plots that illustrate this process - for the 5 most affected genes and for 5 random genes:

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



#### Tools used

- Cellranger:
  - [mkgtf](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf) - filter GTF. [default: off]
  - [mkref](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkref) - create reference
  - [count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) - create feature counts

| ![DAG](https://github.com/CarolinaPB/single-cell-data-processing/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |

### Edit config.yaml with the paths to your files

```yaml
DATA: /path/to/directory/with/fastqs
GTF: /path/to/gtf.gtf
PREFIX: <prefix>

RENAME: <y/n>

CR_COUNT_extra: ""

# QC parameters
MITO_PERCENTAGE: 10
NUMBER_GENES_PER_CELL: 500 
NUMBER_UMI_PER_CELL: 1000
ENSEMBLE_BIOMART_SPECIES: "<species>"
# Doublet removal parameters
# threshold doublet score (should be at the minimum between two modes of the simulated doublet histogram)
# Add a line for every sample, even if you don't add a value.
SCRUB_THRESHOLD: 
  <sample 1>: <value>
  <sample 2>: <empty>
```

- DATA - path to directory containing fastq files. Preferrably, the files should be named in the format accepted by Cellranger Count `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`. If they are, set `RENAME: n`. If not, they should be in the format `<sample>_R1.fastq.gz`. In this case, you should set `RENAME: y` so that the pipeline will rename the files according to the necessary format for Cellranger Count.
- PREFIX - The name of your organism. The reference package used for cellranger count will be in the `<prefix>_genome` directory
- RENAME - `y` if your input fastqs are not named in this format `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`. `n` if they are.
- Options for Cellrange counter:
  - CR_COUNT_extra - any other options for cellranger count. [Find other options here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count). [Default: ""]
- QC parameters
  - MITO_PERCENTAGE - Keep cells with less than X% mitochondrial read fraction. [Default: 10]
  - NUMBER_GENES_PER_CELL - keep cells with more than X genes. [Default: 500]
  - NUMBER_UMI_PER_CELL - keep cells with more than X UMIs. [Default: 1000]
  - ENSEMBLE_BIOMART_SPECIES -  ensembl biomart species used to get the
- Doublet removal score
  - SCRUB_THRESHOLD - threshold doublet score. It should be at the minimum between two modes of the simulated doublet histogram.   
  In the first run it should be run as `SCRUB_TRESHOLD: `.   
  After that is done, for each sample you should then look at the `4_Doublets/<sample>/histogram_<sample>_doublets.pdf` plot and see if the vertical line on the "simulated doublets" plot is at the minimum between the two modes. If it's not, you should manually set it in the config file as:

```
SCRUB_THRESHOLD: 
  <sample 1>: <value>
  <sample 2>: <empty>
```
There should be a line for each sample, even if you don't need to set the threshold for that sample. If you need to change the treshold, set `<sample>: <value>`, if not, set `<sample>: `.




### Other set up
If you're working with human or mouse data, download the reference from here:
<https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest> and place it in a folder in the pipeline directory called `<prefix>_genome`

## RESULTS

The most important files are and directories are:  

# TODO

First run until QC_and_remove_doublets. Look at the histogram of simulated doublets. The plot should have a bimodal distribution. The threshold should be in the dip between the two modes. The first run of the pipeline with use a automatically detected threshold. After the first run you should look at the histogram and see if the treshold is in the dip between the modes. If not, go to the config.yaml file, and add a list with your sample names (you can look at the name of the samples in the QC directory - name of directories.)
the list should be in this format:
Add all the samples, even if you want to keep the automatic thresholds. If you want to use the automatic thresholds don't add anything after the colon (":")

```yaml
SCRUB_THRESHOLD: 
  <sample name>: <threshold value>
  <sample name>: <threshold value>
  <sample name>: 
```

After changing the thresholds, you'll need to rerun this step again. Before you do this, you need to copy the previous results (QC directory) to another directory, or you need to delete those results

Once that's done you can rerun the QC_and_remove_doublets step

The QC_and_remove_doublets step creates several files:
.h5ad file with filtered results
jupyter notebook with the analysis steps and plots. This is is an interactive notebook. it can be opened and run, changed, etc. Use it to explore your data
one directory per sample containing the plots created during this step (also present in the jupyter notebook)

# TODO

ADD MKREF OPTION
