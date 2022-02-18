# Single-cell preprocessing pipeline

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT
This pipeline includes the first steps in the analysis of Single-cell data. It starts doing QC. Then the reads are corrected and trimmed for 3' poly A tail.  
It uses a reference genome and an annotation file (GTF) to build the reference genome index. Cellranger count [performs alignment, filtering, barcode counting, and UMI counting. It uses the Chromium cellular barcodes to generate feature-barcode matrices, determine clusters, and perform gene expression analysis](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). The annotation file was modified to include both gene symbol (if available) and Ensembl ID as gene reference.

#### Tools used:
- Cellranger:
    - [mkgtf](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf) - filter GTF. [default: off]
    - [mkref](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkref) - create reference
    - [count](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count) - create feature counts


| ![DAG](https://github.com/CarolinaPB/nanopore-assembly/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |


### Edit config.yaml with the paths to your files
```
DATA: /path/to/directory/with/fastqs
FASTA: /path/to/fasta.fa
GTF: /path/to/gtf.gtf

# Filter GTF
FILTER_GTF: <y/n> 
ATTRIBUTES:
  - "--attribute=gene_biotype:protein_coding"
  - "--attribute=<attribute>"

# mkref options
PREFIX: <prefix>
REF_VERSION: 
    - "--ref-version=<version>"
CR_MKREF_EXTRA: ""

# rcorrector options
KMER: <kmer length>
RCORRECTOR_EXTRA: "" 
RCORRECTOR: <se/pe>
# if se:
RCORRECTOR_R: <R1/R2>

# polyA trimming options
TRIM_EXTRA: 
  - "<option>"
  - "<option>"

# Cell ranger counter
CR_COUNT_extra: ""
```

- DATA - path to directory containing fastq files (fastq files should have *R1_fastq.gz and *R2.fastq ending)
- FASTA - path to reference genome fasta file
- GTF - path to GTF file
- Options for cellranger mkref
    - PREFIX - prefix to name output folder. The name will be <prefix>_genome
    - REF_VERSION - version of reference genome
    - CR_MKREF_EXTRA - any other parameters for mkref. [Default: ""]
- Options for filtering
    - FILTER_GTF - `Y` for filtering the GTF for the attributes in ATTRIBUTES. `N` for no filtering. [Default: N]
    - ATTRIBUTES - filters to be applied. Check options [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf) [Default: ""]
- Options for fastq correcting
    - KMER - kmer length [Default: 25]
    - RCORRECTOR_EXTRA - any other Rcorrect parameters. [Default: ""]
    - RCORRECTOR - `se` if only want to correct R1 or R2 reads. `pe` if you want to correct both R1 and R2 reads. [Default: se]
    - RCORRECTOR_R: `R1` or `R2` for correcting R1 or R2 reads, respectively. [Default: R2]
- Options for trimming:
    - TRIM_EXTRA - parameters to use. Check options [here](https://cutadapt.readthedocs.io/en/stable/guide.html). [Default: '-A "A{10}"', "--minimum-length :25"]
- Options for Cellrange counter:
    - CR_COUNT_extra - any other cellranger count parameters. [Default: ""]


## RESULTS
The most important files are and directories are:  



# TODO
First run until QC_and_remove_doublets. Look at the histogram of simulated doublets. The plot should have a bimodal distribution. The threshold should be in the dip between the two modes. The first run of the pipeline with use a automatically detected threshold. After the first run you should look at the histogram and see if the treshold is in the dip between the modes. If not, go to the config.yaml file, and add a list with your sample names (you can look at the name of the samples in the QC directory - name of directories.)
the list should be in this format:
Add all the samples, even if you want to keep the automatic thresholds. If you want to use the automatic thresholds don't add anything after the colon (":")
```
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