configfile: "config.yaml"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from pathlib import Path
pipeline = "single-cell-data-processing" # replace with your pipeline's name


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

sample_stem, = glob_wildcards(os.path.join(DATA_DIR, "{sample_stem}_R1.fastq.gz"))

localrules: edit_gtf, create_file_log

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
                            "cloupe.cloupe"]

rule all:
    input:
        "read_quality_assessment/multiqc_report.html",
        expand("polyA_trimmed/{sample_stem}_R2.Atrimmed.fastq.gz", sample_stem = sample_stem),
        expand("polyA_trimmed/{sample_stem}_R1.Atrimmed.fastq.gz", sample_stem = sample_stem),
        "mkref_done.txt",
        expand("{sample_stem}/outs/{counts_out}",counts_out = cellranger_count_outfiles, sample_stem=sample_stem)



rule fastqc:
    input:
        expand(os.path.join(DATA_DIR, "{sample}.fastq.gz"), sample=sample),
    output:
        expand(os.path.join("fastqc_results", "{sample}_fastqc.{ext}"), sample=sample, ext = ["zip", "html"]),
    message:
        'Rule {rule} processing'
    params:
        outdir = "fastqc_results",
        DATA = DATA_DIR
    group:
        'qc'
    shell:
        'fastqc --outdir {params.outdir} --threads 16 {input}'

rule multiqc:
    input:
        expand(os.path.join("fastqc_results", "{sample}_fastqc.{ext}"), sample=sample, ext = ["zip", "html"]),
    output:
        "read_quality_assessment/multiqc_report.html"
    message:
        'Rule {rule} processing'
    params:
        # outdir = "read_quality_assessment,
        fastqc_results = "fastqc_results"
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
            R1 = os.path.join("corrected_reads", "{sample_stem}_R1.cor.fq.gz"),
            R2 = os.path.join("corrected_reads", "{sample_stem}_R2.cor.fq.gz"),
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
            os.path.join("corrected_reads", "{sample_stem}_"+RCORRECTOR_R+".cor.fq.gz"),
        message:
            'Rule {rule} processing'
        params:
            outdir = "corrected_reads",
            kmer = KMER,
            extra = RCORRECTOR_EXTRA
        shell:
            "run_rcorrector.pl {params.extra} -s {input} -k {params.kmer} -od {params.outdir}"


# def rcorrector_input_se():
#     if RCORRECTOR == "se":
#         return(",".join(expand(os.path.join("corrected_reads", "{sample_stem}_R2.cor.fq.gz"), sample_stem =sample_stem)))
#     elif RCORRECTOR == "se":
#         return(",".join(expand(os.path.join("corrected_reads", "{sample_stem}_R2.cor.fq.gz"), sample_stem =sample_stem)))

def input_trim_poly_A(wildcards):
    if RCORRECTOR == "pe": 
        # print({"R1": os.path.join("corrected_reads", "{wildcards.sample_stem}_R1.cor.fq.gz".format(wildcards=wildcards)), 
        # "R2": os.path.join("corrected_reads", "{wildcards.sample_stem}_R2.cor.fq.gz".format(wildcards=wildcards))})
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
        R1 = "polyA_trimmed/{sample_stem}_R1.Atrimmed.fastq.gz",
        R2 = "polyA_trimmed/{sample_stem}_R2.Atrimmed.fastq.gz",
    message:
        'Rule {rule} processing'
    group:
        'group'
    params:
        extra = lambda wc: TRIM_EXTRA
    shell:
        """
        cutadapt --cores 0 {params.extra} -o {output.R1} -p {output.R2} {input.R1} {input.R2}
        """


# cell ranger installation 
# export PATH=/lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2:$  
rule edit_gtf:
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
        gtf = Path(GTF).stem + ".edited.gtf"
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


rule cellranger_count:
    input:
        fastq = rules.trim_poly_A.output,
        mkref = rules.cellranger_mkref.output.done
    output:
        expand("{{sample_stem}}/outs/{counts_out}",counts_out = cellranger_count_outfiles )
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
--fastqs=polyA_trimmed \
--sample={wildcards.sample_stem} \
--jobmode=slurm
        """

rule subset_count:
    input:
        'input'
    output:
        'file.out'
    message:
        'Rule {rule} processing'
    shell:
        ''