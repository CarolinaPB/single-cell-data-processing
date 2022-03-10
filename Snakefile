# configfile: "config.yaml"

import sys
from snakemake.shell import shell
from snakemake.utils import makedirs
# from snakemake_wrapper_utils.java import get_java_opts
from pathlib import Path
pipeline = "single-cell-data-processing" # replace with your pipeline's name


include: "rules/create_file_log.smk"

if "OUTDIR" in config:
    workdir: config["OUTDIR"]


# Get name of config file (to use in the subworkflow)
args = sys.argv
print(args[args.index("--configfile") + 1])
config_path = args[args.index("--configfile") + 1]


makedirs("logs_slurm")

DATA_DIR = config["DATA"]
MKREF = config["MKREF"].lower()
if MKREF == "y":
    REF_VERSION = config["REF_VERSION"]
    CR_MKREF_EXTRA = config["CR_MKREF_EXTRA"]
    PREFIX = config["PREFIX"]
    FASTA = config["FASTA"]
    GTF = config["GTF"]
    ATTRIBUTES = config["ATTRIBUTES"]
    FILTER_GTF = config["FILTER_GTF"].lower()

RENAME = config["RENAME"].lower()

CR_COUNT_EXTRA = config["CR_COUNT_extra"]

# QC parameters
MITO_PERCENTAGE = config["MITO_PERCENTAGE"] # keep cells with less than X% mitochondrial read fraction
NUMBER_GENES_PER_CELL = config["NUMBER_GENES_PER_CELL"] # keep cells with more than X genes
NUMBER_UMI_PER_CELL = config["NUMBER_UMI_PER_CELL"] # keep cells with more than X UMIs
ENSEMBLE_BIOMART_SPECIES = config["ENSEMBLE_BIOMART_SPECIES"] # ensembl biomart species used to get the mitochondrial genes for that species
# Doublet removal - threshold doublet score
SCRUB_THRESHOLD = config['SCRUB_THRESHOLD']


localrules:  create_file_log, remove_ambient_RNA, combine_cellranger_counter_metrics, get_mito_genes, edit_gtf, filter_GTF, QC, remove_doublets

if RENAME.lower() =="n":
    SAMPLES, REST, = glob_wildcards(os.path.join(DATA_DIR, "{samp}_S{rest}_R1_001.fastq.gz"))
elif RENAME.lower() == "y":
    SAMPLES, = glob_wildcards(os.path.join(DATA_DIR, "{samp}_R1.fastq.gz"))


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
# If strings in SAMPLES are in a <sample>/<sample>.fastq format, the name of the directory ("sample") will be used as the sample name
samples_dict = {}
for el in SAMPLES:
    print(el)
    if "/" in el: 
    # if the string corresponds to directory/file.fastq.gz, split the string at the separator "/"
        op = el.split("/")
        # if the sample name op[0], the first element after the split, is already in the dictionary, add the 
        if op[0] in samples_dict:
            print(op[0])
            print(op[1])
            samples_dict[op[0]].append(op[1])
        else:
            print(op[0])
            samples_dict[op[0]] = [op[1]]
    else:
        samples_dict[el] = [el]

rule all:
    input:
        files_log,
        'cellranger_count_metrics_allsamples.tsv',
        expand('4_Doublets/{samples}_QC_doublets.h5ad', samples = samples_dict.keys()),

rule edit_gtf:
    '''
    combine gene symbol and ensembl ID
    '''
    input:
        GTF
    output:
        f"{Path(GTF).stem}.edited.gtf"
    message:
        'Rule {rule} processing'
    shell:
        """
perl -p -e 'if (/gene_name/) {{s{{(gene_id\s+"([^"]+).+?gene_name\s+")([^"]+)}}{{$1$3_$2}}}} \
elsif (!/^#/ && /gene_id/) {{s/(gene_id\s+"([^"]+)";\s+)/$1gene_name "$2"; /}}' {input} > {output}
        """


rule filter_GTF:
    input:
        rules.edit_gtf.output
    output:
        f"{Path(GTF).stem}.filtered.gtf"
    message:
        'Rule {rule} processing'
    params:
        attributes = ATTRIBUTES
    shell:
        """
/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2/cellranger mkgtf \
{input} {output} {params.attributes}
        """

def input_cellranger_mkref():
    """
    Function that returns the file to be used as input for cellranger_mkref.
    If the option to filter the GTF file was set to "y", it will return the filtered GTF. 
    If the option to filter the GTF file was set to "n", it will return the unfiltered GTF.
    """
    if FILTER_GTF.lower() == "y":
        return(f"{Path(GTF).stem}.filtered.gtf")
    elif FILTER_GTF.lower() == "n":
        return(f"{Path(GTF).stem}.edited.gtf")


rule cellranger_mkref: # short run time. around 10 min
    input:
        fasta = FASTA,
        gtf = input_cellranger_mkref()
    output:
        outdir = directory(f"{PREFIX}_genome"),
        done = touch("mkref_done.txt")
    message:
        'Rule {rule} processing'
    params:
        ref_version = REF_VERSION,
        extra = CR_MKREF_EXTRA
    shell:
        """
/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2/cellranger mkref \
--genome={output.outdir} \
--fasta={input.fasta} \
--genes={input.gtf} \
--nthreads=16 \
--ref-version={params.ref_version} \
{params.extra}
        """

"""
Sets the path of the fastqs to be used in the rest of the pipeline
If they are not in the format accepted by cellranger, the option "RENAME: y" should be set
and the fastq files will be renamed (a symbolic link will be created with the new name).
If the fastq files have not been renamed, then their directory is DATA_DIR, if they have, 
the directory will be 1_renamed. If an output directory has been specified in the config file, 
1_renamed will be appended to the outdir path
"""
FASTQS_DIR = DATA_DIR
if RENAME.lower() == "y":
    subworkflow rename_files:
        snakefile:
            os.path.join(workflow.basedir,"subworkflow/rename_fastqs/Snakefile")
        configfile:
            os.path.join(workflow.basedir, config_path) # use same config as used for the main snakemake file

    # Path to renamed fastq files
    if "OUTDIR" in config:
        FASTQS_DIR = os.path.join(config["OUTDIR"], "1_renamed")
    else:
        FASTQS_DIR = "1_renamed"

def get_cellranger_count_input_rename(wildcards):
    """
    If the fastq files have to be renamed, this function makes sure that 
    cellranger_count only runs after they have been renamed by giving the 
    result of the renaming rule as input
    """
    if RENAME.lower() == "y":
        return(rename_files(f"renamed_{wildcards.samples}.done"))
    elif RENAME.lower() == "n":
        return([])

def get_cellranger_count_input_mkref():
    """
    If the mkref step is used, this function will make sure cellranger_count only
    runs after mkref is finished
    """
    if MKREF == "y":
        return("mkref_done.txt")
    elif MKREF == "n":
        return([])

"""
Cellranger count:
Can't specify output of the rule with the actual output files since cellranger count gives the error:
RuntimeError: <directory> is not a pipestance directory
if the cellranger count output directory was not created by cellranger count. 
By adding the output files to the output of the rule, snakemake would create the necessary
directories, and therefore trigger the error
One option would be to keep the complete snakemake output and do "rm -r {{samples}}" before the
cellranger_count rule
"""

rule cellranger_count:
    input:
        get_cellranger_count_input_rename,
        get_cellranger_count_input_mkref()
    output:
        # expand("{{samples}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
        # directory(expand("{{samples}}/outs/{counts_outdirs}",counts_outdirs = cellranger_count_outdirs)),
        touch("steps_done/cellranger_count_{samples}.done")
    message:
        'Rule {rule} processing'
    params:
        extra = CR_COUNT_EXTRA,
        transcriptome = os.path.join(workflow.basedir, f"{PREFIX}_genome"),
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
    """
    For each sample, creates a notebook, with the analysis steps and results. 
    Results are saved to SoupX/<sample>
    These files will be created: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
    """
    input:
        # "{samples}/outs/"
        # expand("{{samples}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
        rules.cellranger_count.output
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


rule combine_cellranger_counter_metrics:
    """
    Combines cellranger count metrics for all samples in a single table
    """
    input:
        # expand("{samples}/outs/{counts_out}",samples = samples_dict.keys(),counts_out = cellranger_count_outfiles ),
        expand("steps_done/cellranger_count_{samples}.done", samples = samples_dict.keys())
    output:
        'cellranger_count_metrics_allsamples.tsv'
    message:
        'Rule {rule} processing'
    script:
        'scripts/combine_cellrange_counter_metrics.R'

rule get_mito_genes:
    """
    Extract mitochondrial gene names for your species from Ensembl
    Will create two files: one with ensembl gene symbols and the other with gene IDs
    """
    output:
        ensembl = f"{ENSEMBLE_BIOMART_SPECIES}_mito_genes_ensembl.csv",
        symbol = f'{ENSEMBLE_BIOMART_SPECIES}_mito_genes_symbol.csv'
    message:
        'Rule {rule} processing'
    params:
        species =  ENSEMBLE_BIOMART_SPECIES
    run:
        import pandas as pd
        import scanpy as sc
        mito_ensembl_ids = sc.queries.mitochondrial_genes(params.species, attrname="ensembl_gene_id")
        mito_gene_ids = sc.queries.mitochondrial_genes(params.species, attrname="external_gene_name")
        mito_ensembl_ids.to_csv(output[0])  
        mito_gene_ids.to_csv(output[1])

def set_scrub_treshold(wildcards):
    """
    If SCRUB_THRESHOLD is set, return the treshold per sample
    """
    if not SCRUB_THRESHOLD:
        return("")
    else:
        return(SCRUB_THRESHOLD[wildcards.samples])


rule QC:
    """
    For each sample, creates a notebook with QC steps and results.
    It will also create a data object file for each sample (_QC.h5ad) and
    a directory with the plots created
    """
    input:
        ambient_RNA = rules.remove_ambient_RNA.output,
        mito_genes_ensembl = rules.get_mito_genes.output.ensembl,
        mito_genes_symbol = rules.get_mito_genes.output.symbol
    output:
        "3_QC/{samples}_QC.h5ad"
    log:
        notebook = "3_QC/processed_notebook_{samples}.ipynb"
    params:
        mito_percentage = MITO_PERCENTAGE,
        number_genes_per_cell = NUMBER_GENES_PER_CELL,
        number_UMI_per_cell = NUMBER_UMI_PER_CELL,
        ensemble_biomart_species = ENSEMBLE_BIOMART_SPECIES,
        sample = "{samples}",
    message:
        'Rule {rule} processing'
    notebook:
        'scripts/QC_Scanpy.py.ipynb'


rule remove_doublets:
    """
    For each sample creates a notebook with steps to remove doublets and results.
    It will create a data object file (_doublets.h5ad) and a directory with the plots created
    """
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
        scrub_threshold = lambda wildcards: set_scrub_treshold(wildcards) #treshold per sample defined in the configfile
    notebook:
        'scripts/Doublet_removal.py.ipynb'

onsuccess:
    print("Workflow finished")