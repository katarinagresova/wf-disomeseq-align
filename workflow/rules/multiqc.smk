# One MultiQC report: FastQC of the raw reads, cutadapt, FastQC of the trimmed
# reads. Sections and sample names: workflow/multiqc_config.yaml.
MULTIQC_REPORTS = [
    f"{RESULTS_DIR}/fastqc/raw/{{sample}}_fastqc.zip",
    f"{LOG_DIR}/cutadapt/{{sample}}.log",
    f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_fastqc.zip",
]


rule multiqc:
    input:
        reports=expand(MULTIQC_REPORTS, sample=samples.index),
        config=workflow.source_path("../multiqc_config.yaml"),
    output:
        html=f"{RESULTS_DIR}/multiqc/multiqc_report.html",
        data=directory(f"{RESULTS_DIR}/multiqc/multiqc_data"),
    log:
        f"{LOG_DIR}/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --config {input.config} --outdir $(dirname {output.html}) "
        "{input.reports} > {log} 2>&1"
