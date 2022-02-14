configfile: "config.yaml"

from snakemake.shell import shell
# from snakemake_wrapper_utils.java import get_java_opts
from pathlib import Path
pipeline = "single-cell-data-processing" # replace with your pipeline's name


include: "rules/create_file_log.smk"

if "OUTDIR" in config:
    workdir: config["OUTDIR"]

DATA_DIR = config["DATA"]
PREFIX = config["PREFIX"]
CR_COUNT_EXTRA = config["CR_COUNT_extra"]
RENAME = config["RENAME"]


localrules:  create_file_log, remove_ambient_RNA, combine_cellrange_counter_metrics, rename_files

# SAMPLES_SHORT,SAMPLES, REST, = glob_wildcards(os.path.join(DATA_DIR, "{sample_short}/{sample}_S{rest}_R1_001.fastq.gz"))

if RENAME.lower() =="n":
    TEST, REST, = glob_wildcards(os.path.join(DATA_DIR, "{test}_S{rest}_R1_001.fastq.gz"))
elif RENAME.lower() == "y":
    TEST, = glob_wildcards(os.path.join(DATA_DIR, "{test}_R1.fastq.gz"))


print(TEST)

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

test_dict = {}
for el in TEST:
    if "/" in el:
        op = el.split("/")
        if op[0] in test_dict:
            test_dict[op[0]].append(op[1])
        else:
            test_dict[op[0]] = [op[1]]
    else:
        test_dict[el] = [el]
print(test_dict)

rule all:
    input:
        files_log,
        # expand("{samples}/outs/{counts_out}",counts_out = cellranger_count_outfiles, samples = test_dict.keys()),
        # expand("cellranger_count_{samples}.done", samples = test_dict.keys()),
        expand('ambient_RNA_correction/Ambient_RNA_correction_{samples}.html', samples = test_dict.keys()),
        'cellranger_count_metrics_allsamples.tsv'

        # expand("remove_ambient_RNA_{samples}.done", samples = test_dict.keys()),


# def cellranger_count_input():
if RENAME.lower() == "y":
    subworkflow rename_files:
        # workdir:
        #     os.path.join(workflow.basedir,"subworkflow/rename_fastqs")
        snakefile:
            os.path.join(workflow.basedir,"subworkflow/rename_fastqs/Snakefile")
        configfile:
            os.path.join(workflow.basedir,"subworkflow/rename_fastqs/config.yaml")

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
        # expand("{{samples}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
        touch("cellranger_count_{samples}.done")
    message:
        'Rule {rule} processing'
    params:
        extra = CR_COUNT_EXTRA,
        transcriptome = PREFIX + "_genome",
        # fastqs = lambda wildcards: os.path.join(DATA_DIR, test_dict[wildcards.samples]),
        fastqs = DATA_DIR,
        samples = lambda wildcards: ",".join(test_dict[wildcards.samples]),
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
        "cellranger_count_{samples}.done"
    output:
        # "remove_ambient_RNA_{samples}.done"
        'ambient_RNA_correction/Ambient_RNA_correction_{samples}.html'
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
        # expand("{samples}/outs/{counts_out}",samples = test_dict.keys(),counts_out = cellranger_count_outfiles ),
        expand("cellranger_count_{samples}.done", samples = test_dict.keys())
    output:
        'cellranger_count_metrics_allsamples.tsv'
    message:
        'Rule {rule} processing'
    script:
        'scripts/combine_cellrange_counter_metrics.R'












# number_of_samples = len(sample_stem)
# rule all:
#     input:
#         # "read_quality_assessment/multiqc_report.html",
#         # "read_quality_assessment/trimmed/multiqc_report.html",
#         files_log,
#         # 'cellranger_count_metrics_allsamples.tsv',
#         expand("{sample}/outs/{counts_out}",counts_out = cellranger_count_outfiles, sample=sample),

#         # expand('SoupX/Ambient_RNA_correction_{sample_dir}.html', sample_dir = samples_dir)


# output explanation: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/gex-outputs
# rule cellranger_count:
#     input:
#         # fastq = os.path.join(DATA_DIR, "{sample}_R1_001.fastq.gz"),
#         ref_dir = PREFIX + "_genome" 
#     output:
#         expand("{{sample_dir}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
#         # touch("count_{sample_stem}_{sample_long}.done")
#     message:
#         'Rule {rule} processing'
#     params:
#         # transcriptome = rules.cellranger_mkref.output.outdir,
#         transcriptome = PREFIX + "_genome",
#         extra = CR_COUNT_EXTRA,
#         fastqs = os.path.join(DATA_DIR, "{sample}"),
#         sample = "{sample}" # prefix fastq
#     # wildcard_constraints:
#     #     sample = "[^/]+"
#     shell:
#         """
# /lustre/nobackup/WUR/ABGC/moiti001/TOOLS/cellranger-6.1.2/cellranger count \
# {params.extra} \
# --id={params.sample} \
# --transcriptome={params.transcriptome} \
# --fastqs={params.fastqs} \
# --sample={params.sample}
#         """


# rule remove_ambient_RNA:
#     input:
#         "{sample, [^/]+}/outs/"
#     output:
#         'SoupX/Ambient_RNA_correction_{sample_dir}.html'
#     message:
#         'Rule {rule} processing'
#     params:
#         wd = os.getcwd()
#     script:
#         'scripts/remove_ambient_RNA.Rmd'

# rule fastqc:
#     input:
#         expand(os.path.join(DATA_DIR, "{sample_stem}/{sample_long}.fastq.gz"), sample_stem=sample_stem, sample_long = sample_long),
#     output:
#         expand(os.path.join("1_fastqc_results", "{sample_stem}/{sample_long}_fastqc.{ext}"), sample_stem=sample_stem, sample_long = sample_long, ext = ["zip", "html"]),
#     message:
#         'Rule {rule} processing'
#     params:
#         outdir = "1_fastqc_results",
#     group:
#         'qc'
#     shell:
#         'fastqc --outdir {params.outdir} --threads 16 {input}'

# rule multiqc:
#     input:
#         expand(os.path.join("1_fastqc_results", "{sample_stem}/{sample_long}_fastqc.{ext}"), sample_stem=sample_stem, sample_long = sample_long, ext = ["zip", "html"]),
#     output:
#         "1_fastqc_results/multiqc_report.html" # change to "1_fastq_results/not_trimmed/multiqc_report.html"
#     message:
#         'Rule {rule} processing'
#     params:
#         # outdir = "read_quality_assessment,
#         fastqc_results = "fastqc_results" # change to "fastq_results/not_trimmed"
#     group:
#         'qc'
#     shell:
#         'multiqc {params.fastqc_results} --filename {output}'
