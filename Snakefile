configfile: "config.yaml"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from pathlib import Path
pipeline = "single-cell-data-processing" # replace with your pipeline's name


## cell ranger count. seems like it will pick up where it left off

include: "rules/create_file_log.smk"

if "OUTDIR" in config:
    workdir: config["OUTDIR"]

DATA_DIR = config["DATA"]
GTF = config["GTF"]
FILTER_GTF = config["FILTER_GTF"]
FASTA = config["FASTA"]
PREFIX = config["PREFIX"]
KMER = config["KMER"]
RCORRECTOR = config["RCORRECTOR"]
RCORRECTOR_EXTRA = config["RCORRECTOR_EXTRA"]
if "RCORRECTOR" in config:
    RCORRECTOR_R = config["RCORRECTOR_R"]
TRIM_EXTRA = config["TRIM_EXTRA"]
REF_VERSION = config["REF_VERSION"]
CR_MKREF_EXTRA = config["CR_MKREF_EXTRA"]
CR_COUNT_EXTRA = config["CR_COUNT_extra"]

sample, = glob_wildcards(os.path.join(DATA_DIR, "{sample,[^/]+}.fastq.gz"))

# sample_stem, = glob_wildcards(os.path.join(DATA_DIR, "{sample_stem}_R1.fastq.gz"))
sample_stem, = glob_wildcards(os.path.join(DATA_DIR, "{sample_stem}_R1_001.fastq.gz"))

print(sample_stem)
localrules: edit_gtf, create_file_log, rename_trimmed, combine_cellrange_counter_metrics, remove_ambient_RNA

cellranger_count_outfiles = ["web_summary.html",
                            "metrics_summary.csv", 
                            "possorted_genome_bam.bam",
                            "possorted_genome_bam.bam.bai",
                            "filtered_feature_bc_matrix",
                            "filtered_feature_bc_matrix_h5.h5",
                            "raw_feature_bc_matrix",
                            "raw_feature_bc_matric_h5.h5",
                            "analysis",
                            "molecule_info.h5",
                            "cloupe.cloupe",
                            ""]

# number_of_samples = len(sample_stem)
rule all:
    input:
        "read_quality_assessment/multiqc_report.html",
        "read_quality_assessment/trimmed/multiqc_report.html",
        files_log,
        'cellranger_count_metrics_allsamples.tsv',
        expand('SoupX/Ambient_RNA_correction_{sample_stem}.html', sample_stem = sample_stem)

        # expand("renamed_trimmed/trimmed_renamed_{sample_stem}.done", sample_stem= sample_stem),
        # expand("count_{sample_stem}.done", sample_stem = sample_stem),
        # "polyA_trimmed/{{sample_stem}}_S{n}_L001_R1_001.fastq.gz"
        # expand("fastqc_results/trimmed/{sample}_fastqc.{ext}", sample=sample, ext = ["zip", "html"]),
        # expand(expand("polyA_trimmed/{sample_stem}_S{n}_L001_R1_001.fastq.gz", zip,sample_stem = sample_stem, n = range(1, int(number_samples)+1))),
        # expand(expand("polyA_trimmed/{sample_stem}_S{n}_L001_R2_001.fastq.gz", zip,sample_stem = sample_stem, n = range(1, int(number_samples)+1))),
        # "mkref_done.txt",
        # expand("{sample_stem}/outs/{counts_out}",counts_out = cellranger_count_outfiles, sample_stem=sample_stem),

rule fastqc:
    input:
        expand(os.path.join(DATA_DIR, "{sample}.fastq.gz"), sample=sample),
    output:
        expand(os.path.join("fastqc_results", "{sample}_fastqc.{ext}"), sample=sample, ext = ["zip", "html"]),
    message:
        'Rule {rule} processing'
    params:
        outdir = "fastqc_results", # change to "1_fastq_results/not_trimmed"
    group:
        'qc'
    shell:
        'fastqc --outdir {params.outdir} --threads 16 {input}'

rule multiqc:
    input:
        expand(os.path.join("fastqc_results", "{sample}_fastqc.{ext}"), sample=sample, ext = ["zip", "html"]),
    output:
        "read_quality_assessment/multiqc_report.html" # change to "1_fastq_results/not_trimmed/multiqc_report.html"
    message:
        'Rule {rule} processing'
    params:
        # outdir = "read_quality_assessment,
        fastqc_results = "fastqc_results" # change to "fastq_results/not_trimmed"
    group:
        'qc'
    shell:
        'multiqc {params.fastqc_results} --filename {output}'


if RCORRECTOR == "pe":
    rule rcorrector_pe:
        input:
            R1 = os.path.join(DATA_DIR, "{sample_stem}_R1.fastq.gz"),
            R2 = os.path.join(DATA_DIR, "{sample_stem}_R2.fastq.gz"),
        output:
            R1 = temp(os.path.join("corrected_reads", "{sample_stem}_R1.cor.fq.gz")), # change directory to "2_corrected_reads"
            R2 = temp(os.path.join("corrected_reads", "{sample_stem}_R2.cor.fq.gz")), # change directory to "2_corrected_reads"
        message:
            'Rule {rule} processing'
        params:
            outdir = "corrected_reads",
            kmer = KMER,
            extra = RCORRECTOR_EXTRA
        group:
            'group'
        shell:
            'run_rcorrector.pl {params.extra} -k {params.kmer} -1 {input.R1} -2 {input.R2} -od {params.outdir} -t 16'

elif RCORRECTOR == "se":
    rule rcorrector_se:
        input:
            os.path.join(DATA_DIR, "{sample_stem}_"+RCORRECTOR_R+".fastq.gz"),
        output:
            temp(os.path.join("corrected_reads", "{sample_stem}_"+RCORRECTOR_R+".cor.fq.gz")), # change directory to "2_corrected_reads"
        message:
            'Rule {rule} processing'
        params:
            outdir = "corrected_reads",
            kmer = KMER,
            extra = RCORRECTOR_EXTRA
        group:
            'group'
        shell:
            "run_rcorrector.pl {params.extra} -s {input} -k {params.kmer} -od {params.outdir} -t 16"



def input_trim_poly_A(wildcards):
    if RCORRECTOR == "pe": 
        return({"R1": os.path.join("corrected_reads", "{wildcards.sample_stem}_R1.cor.fq.gz".format(wildcards=wildcards)), 
        "R2": os.path.join("corrected_reads", "{wildcards.sample_stem}_R2.cor.fq.gz".format(wildcards=wildcards))})
    elif RCORRECTOR == "se":
        if RCORRECTOR_R == "R1":
            return({"R1": os.path.join("corrected_reads", "{wildcards.sample_stem}_R1.cor.fq.gz".format(wildcards=wildcards)),
            "R2": os.path.join(DATA_DIR, "{wildcards.sample_stem}_R2.fastq.gz".format(wildcards=wildcards))})
        elif RCORRECTOR_R == "R2":
            return({"R2": os.path.join("corrected_reads", "{wildcards.sample_stem}_R2.cor.fq.gz".format(wildcards=wildcards)),
            "R1": os.path.join(DATA_DIR, "{wildcards.sample_stem}_R1.fastq.gz".format(wildcards=wildcards))})


rule trim_poly_A:
    input:
        unpack(input_trim_poly_A)
    output:
        R1 = "polyA_trimmed/{sample_stem}_R1.fastq.gz", # rename dir to "3_polyA_trimmed"
        R2 = "polyA_trimmed/{sample_stem}_R2.fastq.gz", # rename dir to "3_polyA_trimmed"
        # done = touch("trim_polyA_{sample_stem}.done")
    message:
        'Rule {rule} processing'
    group:
        'group'
    params:
        extra = lambda wc: TRIM_EXTRA,
    shell:
        """
        cutadapt --cores 0 {params.extra} -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        """

def create_sample_numbered_dict(sample_list):
    sample_dict = {}
    i = 1
    for element in sample_list:
        sample_dict[element] = i
        i +=1
    return(sample_dict)

numbered_samples = create_sample_numbered_dict(sample_stem)

rule rename_trimmed:
    input:
        R1 = rules.trim_poly_A.output.R1,
        R2 = rules.trim_poly_A.output.R2
    output:
        # "polyA_trimmed/{{sample_stem}}_S{n}_L001_R1_001.fastq.gz".format(n=numbered_samples[wildcards.sample_stem]),
        # "polyA_trimmed/{{sample_stem}}_S{n}_L001_R2_001.fastq.gz".format(n=numbered_samples[wildcards.sample_stem]),
        touch("renamed_trimmed/trimmed_renamed_{sample_stem}.done") # change dir to "5_renamed_dir"
    message:
        'Rule {rule} processing'
    params: 
        number = lambda wildcards: numbered_samples["{}".format(wildcards.sample_stem)],
        outdir = "renamed_trimmed" # change dir to "5_renamed_dir"
    shell:
        """
        cp {input.R1} {params.outdir}/{wildcards.sample_stem}_S{params.number}_L001_R1_001.fastq.gz
        cp {input.R2} {params.outdir}/{wildcards.sample_stem}_S{params.number}_L001_R2_001.fastq.gz
        """

use rule fastqc as fastqc_trimmed with:
    input:
        expand("polyA_trimmed/{sample}.fastq.gz", sample=sample)
    output:
        expand("fastqc_results/trimmed/{sample}_fastqc.{ext}", sample=sample, ext = ["zip", "html"]), # change dir to 1_fastq_results
    params:
        outdir = "fastqc_results/trimmed" # change dir to 1_fastq_results

use rule multiqc as multiqc_trimmed with:
    input:
        expand("fastqc_results/trimmed/{sample}_fastqc.{ext}", sample=sample, ext = ["zip", "html"]),
    output:
        "read_quality_assessment/trimmed/multiqc_report.html" # change to "1_fastq_results/trimmed/multiqc_report.html"
    params:
        fastqc_results = "fastqc_results/trimmed" # change dir to 1_fastq_results

 
rule edit_gtf:
    '''
    combine gene symbol and ensembl ID
    '''
    input:
        GTF
    output:
        Path(GTF).stem + ".edited.gtf"
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
        Path(GTF).stem + ".filtered.gtf"
    message:
        'Rule {rule} processing'
    params:
        attributes = config["ATTRIBUTES"]
    shell:
        """
/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2/cellranger mkgtf \
{input} {output} {params.attributes}
        """

def input_cellranger_mkref():
    if FILTER_GTF.lower() == "y":
        return(Path(GTF).stem + ".filtered.gtf")
    elif FILTER_GTF.lower() == "n":
        return(Path(GTF).stem + ".edited.gtf")


rule cellranger_mkref: # short run time. around 10 min
    input:
        fasta = FASTA,
        # gtf = Path(GTF).stem + ".edited.gtf"
        gtf = input_cellranger_mkref()
    output:
        outdir = directory(PREFIX + "_genome"),
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

# output explanation: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/gex-outputs
rule cellranger_count:
    input:
        fastq = rules.rename_trimmed.output,
        mkref = rules.cellranger_mkref.output.done
    output:
        expand("{{sample_stem}}/outs/{counts_out}",counts_out = cellranger_count_outfiles ),
        touch("count_{sample_stem}.done")
    message:
        'Rule {rule} processing'
    params:
        transcriptome = rules.cellranger_mkref.output.outdir,
        extra = CR_COUNT_EXTRA
    shell:
        """
/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2/cellranger count \
{params.extra} \
--id={wildcards.sample_stem} \
--transcriptome={params.transcriptome} \
--fastqs=renamed_trimmed \
--sample={wildcards.sample_stem}
        """

rule combine_cellrange_counter_metrics:
    input:
        expand("count_{sample_stem}.done", sample_stem = sample_stem)
    output:
        'cellranger_count_metrics_allsamples.tsv'
    message:
        'Rule {rule} processing'
    script:
        'scripts/combine_cellrange_counter_metrics.R'

rule remove_ambient_RNA:
    input:
        "{sample_stem}/outs/"
    output:
        'SoupX/Ambient_RNA_correction_{sample_stem}.html'
    message:
        'Rule {rule} processing'
    params:
        wd = os.getcwd()
    script:
        'scripts/remove_ambient_RNA.Rmd'

