# FastQC of the raw reads and of the trimmed reads (what bowtie_align aligns:
# adapter-trimmed, collapsed, random nt removed). Reports only; no rule reads
# them.
#
# FastQC names its report after the input file, so each job links its fastq
# into its own tmp dir as <sample>.fastq.gz and moves the report out. The tmp
# dir also takes FastQC's temporary files (--dir): compute nodes' /tmp can be full.
FASTQC_SHELL = (
    "rm -rf {params.tmp} && mkdir -p {params.tmp} && "
    "ln -s $(readlink -f {input.fastq}) {params.tmp}/{wildcards.sample}.fastq.gz && "
    "fastqc --quiet --threads 1 --dir {params.tmp} --outdir {params.tmp} "
    "{params.tmp}/{wildcards.sample}.fastq.gz > {log} 2>&1 && "
    "mv {params.tmp}/{wildcards.sample}_fastqc.html {output.html} && "
    "mv {params.tmp}/{wildcards.sample}_fastqc.zip {output.zip} && "
    "rm -rf {params.tmp}"
)


rule fastqc_raw:
    input:
        fastq=lambda wc: samples.loc[wc.sample, "fastq_1"],
    output:
        html=f"{RESULTS_DIR}/fastqc/raw/{{sample}}_fastqc.html",
        zip=f"{RESULTS_DIR}/fastqc/raw/{{sample}}_fastqc.zip",
    params:
        tmp=f"{RESULTS_DIR}/fastqc/raw/{{sample}}_tmp",
    log:
        f"{LOG_DIR}/fastqc/raw/{{sample}}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        FASTQC_SHELL


rule fastqc_trimmed:
    input:
        fastq=f"{RESULTS_DIR}/trim_reads/{{sample}}.fastq.gz",
    output:
        html=f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_fastqc.html",
        zip=f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_fastqc.zip",
    params:
        tmp=f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_tmp",
    log:
        f"{LOG_DIR}/fastqc/trimmed/{{sample}}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        FASTQC_SHELL
