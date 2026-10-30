# One MultiQC report: FastQC of the raw reads, cutadapt, the QC summary and the
# read fates (qc_summary), FastQC of the trimmed reads. Sections and sample
# names: workflow/multiqc_config.yaml.
MULTIQC_REPORTS = [
    f"{RESULTS_DIR}/fastqc/raw/{{sample}}_fastqc.zip",
    f"{LOG_DIR}/cutadapt/{{sample}}.log",
    f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_fastqc.zip",
]


rule multiqc:
    input:
        reports=expand(MULTIQC_REPORTS, sample=samples.index),
        qc=[f"{RESULTS_DIR}/qc/summary.tsv", f"{RESULTS_DIR}/qc/read_fates.tsv"],
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
        "{input.reports} {input.qc} > {log} 2>&1"
