configfile: "config.yaml"

from snakemake.shell import shell
from snakemake.utils import makedirs
# from snakemake_wrapper_utils.java import get_java_opts
from pathlib import Path
pipeline = "single-cell-data-processing" # replace with your pipeline's name


include: "rules/create_file_log.smk"

if "OUTDIR" in config:
    workdir: config["OUTDIR"]

makedirs("logs_slurm")

DATA_DIR = config["DATA"]
PREFIX = config["PREFIX"]
CR_COUNT_EXTRA = config["CR_COUNT_extra"]
RENAME = config["RENAME"]

MITO_PERCENTAGE = config["MITO_PERCENTAGE"] # keep cells with less than X% mitochondrial read fraction
NUMBER_GENES_PER_CELL = config["NUMBER_GENES_PER_CELL"] # keep cells with more than X genes
NUMBER_UMI_PER_CELL = config["NUMBER_UMI_PER_CELL"] # keep cells with more than X UMIs
ENSEMBLE_BIOMART_SPECIES = config["ENSEMBLE_BIOMART_SPECIES"] # ensembl biomart species used to get the mitochondrial genes for that species
SCRUB_THRESHOLD = config['SCRUB_THRESHOLD']





localrules:  create_file_log, remove_ambient_RNA, combine_cellrange_counter_metrics


if RENAME.lower() =="n":
    SAMPLES, REST, = glob_wildcards(os.path.join(DATA_DIR, "{test}_S{rest}_R1_001.fastq.gz"))
elif RENAME.lower() == "y":
    SAMPLES, = glob_wildcards(os.path.join(DATA_DIR, "{test}_R1.fastq.gz"))


cellranger_count_outfiles = ["web_summary.html",
                            "metrics_summary.csv", 
                            "possorted_genome_bam.bam",
                            "possorted_genome_bam.bam.bai",
                            "filtered_feature_bc_matrix.h5",
                            "raw_feature_bc_matrix.h5",
                            "molecule_info.h5",
                            "cloupe.cloupe"]
cellranger_count_outdirs = ["filtered_feature_bc_matrix",
                            "raw_feature_bc_matrix",
                            "analysis"]
# Create dictionary with sample names 
samples_dict = {}
for el in SAMPLES:
    if "/" in el:
        op = el.split("/")
        if op[0] in samples_dict:
            samples_dict[op[0]].append(op[1])
        else:
            samples_dict[op[0]] = [op[1]]
    else:
        samples_dict[el] = [el]

print(samples_dict)
rule all:
    input:
        files_log,
        'cellranger_count_metrics_allsamples.tsv',
        expand('4_Doublets/{samples}_QC_doublets.h5ad', samples = samples_dict.keys()),

        # expand('2_ambient_RNA_correction/Ambient_RNA_correction_{samples}.html', samples = samples_dict.keys()),
        # expand("3_QC/{samples}_QC.h5ad", samples = samples_dict.keys())
        # expand("SoupX/{samples}/features.tsv", samples = samples_dict.keys()), # remove
        # expand("remove_ambient_RNA_{samples}.done", samples = samples_dict.keys()),
        # expand("{samples}/outs/{counts_out}",counts_out = cellranger_count_outfiles, samples = samples_dict.keys()),
        # expand("cellranger_count_{samples}.done", samples = samples_dict.keys()),




FASTQS_DIR = DATA_DIR
if RENAME.lower() == "y":
    subworkflow rename_files:
        workdir:
            config["OUTDIR"]
        snakefile:
            os.path.join(workflow.basedir,"subworkflow/rename_fastqs/Snakefile")
        configfile:
            os.path.join(workflow.basedir,"config.yaml")

    FASTQS_DIR = os.path.join(config["OUTDIR"], "renamed")

def cellranger_count_input(wildcards):
    if RENAME.lower() == "y":
        return(rename_files(f"renamed_{wildcards.samples}.done"))
    elif RENAME.lower() == "n":
        return([])


rule cellranger_count:
    input:
        cellranger_count_input
    #     directory(os.path.join(DATA_DIR, "{samples_short}"))
    output:
        expand("{{samples}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
        directory(expand("{{samples}}/outs/{counts_outdirs}",counts_outdirs = cellranger_count_outdirs)),
        # touch("steps_done/cellranger_count_{samples}.done")
    message:
        'Rule {rule} processing'
    params:
        extra = CR_COUNT_EXTRA,
        transcriptome = os.path.join(workflow.basedir, PREFIX + "_genome"),
        # fastqs = lambda wildcards: os.path.join(DATA_DIR, samples_dict[wildcards.samples]),
        fastqs = FASTQS_DIR,
        samples = lambda wildcards: ",".join(samples_dict[wildcards.samples]),
    shell:
        """
/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2/cellranger count \
{params.extra} \
--id={wildcards.samples} \
--transcriptome={params.transcriptome} \
--fastqs={params.fastqs} \
--sample={params.samples}
        """


rule remove_ambient_RNA:
    input:
        # "{samples}/outs/"
        # "cellranger_count_{samples}.done"
        expand("{{samples}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
    output:
        '2_ambient_RNA_correction/Ambient_RNA_correction_{samples}.html',
        # expand("SoupX/{{samples}}/{soupxfile}", soupxfile = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]), 
    message:
        'Rule {rule} processing'
    params:
        wd = os.getcwd(),
        input = "{samples}/outs/",
        sample = "{samples}"
    script:
        'scripts/remove_ambient_RNA.Rmd'


rule combine_cellrange_counter_metrics:
    input:
        # expand("count_{sample}.done", sample_stem = sample_stem, sample_long = sample_long)
        # expand("count_{sample_stem}.done", sample_stem = sample_stem)
        expand("{samples}/outs/{counts_out}",samples = samples_dict.keys(),counts_out = cellranger_count_outfiles ),
        # expand("cellranger_count_{samples}.done", samples = samples_dict.keys())
    output:
        'cellranger_count_metrics_allsamples.tsv'
    message:
        'Rule {rule} processing'
    script:
        'scripts/combine_cellrange_counter_metrics.R'

rule get_mito_genes:
    output:
        f"{ENSEMBLE_BIOMART_SPECIES}_mito_genes.csv"
    message:
        'Rule {rule} processing'
    params:
        species =  ENSEMBLE_BIOMART_SPECIES
    run:
        import pandas as pd
        import scanpy as sc
        mito_ensembl_ids = sc.queries.mitochondrial_genes(params.species, attrname="ensembl_gene_id")
        mito_ensembl_ids.to_csv(output[0])  

def set_scrub_treshold(wildcards):
    if not SCRUB_THRESHOLD:
        return("")
    else:
        return(SCRUB_THRESHOLD[wildcards.samples])


rule QC:
    input:
        amibient_RNA = rules.remove_ambient_RNA.output,
        mito_genes = rules.get_mito_genes.output
    output:
        "3_QC/{samples}_QC.h5ad"
    log:
        notebook = "3_QC/processed_notebook_{samples}.ipynb"
    params:
        mito_percentage = 10,
        number_genes_per_cell = 500,
        number_UMI_per_cell = 1000,
        ensemble_biomart_species = "sscrofa",
        sample = "{samples}",
        # scrub_threshold = lambda wildcards: SCRUB_THRESHOLD[wildcards.samples]
    message:
        'Rule {rule} processing'
    notebook:
        'scripts/QC_Scanpy.py.ipynb'


rule remove_doublets:
    input:
        rules.QC.output
    output:
        '4_Doublets/{samples}_QC_doublets.h5ad'
    log:
        notebook = "4_Doublets/processed_notebook_{samples}.ipynb"
    message:
        'Rule {rule} processing'
    params:
        sample = "{samples}",
        scrub_threshold = lambda wildcards: set_scrub_treshold(wildcards)
    notebook:
        'scripts/Doublet_removal.py.ipynb'