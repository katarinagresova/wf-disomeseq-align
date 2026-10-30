# QC: footprint positions of each library (qc_library.py), then one row for
# each library (qc_summary.py).


rule qc_library:
    input:
        bam=f"{RESULTS_DIR}/bowtie/{{sample}}.bam",
        transcripts=f"{RESULTS_DIR}/reference/transcripts.tsv",
        script=workflow.source_path("../scripts/qc_library.py"),
    output:
        regions=f"{RESULTS_DIR}/qc/{{sample}}.regions.tsv",
        metagene=f"{RESULTS_DIR}/qc/{{sample}}.metagene.tsv",
        genes=f"{RESULTS_DIR}/qc/{{sample}}.genes.tsv",
        dups=f"{RESULTS_DIR}/qc/{{sample}}.dups.tsv",
    log:
        f"{LOG_DIR}/qc/{{sample}}.log",
    params:
        disome=config["disome_lengths"],
    conda:
        "../envs/pysam.yaml"
    shell:
        "python {input.script} {input.bam} {input.transcripts} {params.disome[0]} {params.disome[1]} "
        "{output.regions} {output.metagene} {output.genes} {output.dups} 2> {log}"


QC_LIBRARY_FILES = [
    f"{LOG_DIR}/cutadapt/{{sample}}.log",
    f"{LOG_DIR}/collapse_reads/{{sample}}.log",
    f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv",
    f"{RESULTS_DIR}/qc/{{sample}}.regions.tsv",
]


rule qc_summary:
    input:
        files=expand(QC_LIBRARY_FILES, sample=samples.index),
        samples=config["samples"],
        script=workflow.source_path("../scripts/qc_summary.py"),
    output:
        summary=f"{RESULTS_DIR}/qc/summary.tsv",
        fates=f"{RESULTS_DIR}/qc/read_fates.tsv",
    log:
        f"{LOG_DIR}/qc/qc_summary.log",
    params:
        disome=config["disome_lengths"],
        gel=config["gel_lengths"],
        # <sample>:<file>:<file>... in QC_LIBRARY_FILES order
        libraries=" ".join(":".join([s] + expand(QC_LIBRARY_FILES, sample=s)) for s in samples.index),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {output.summary} {output.fates} {input.samples} {params.disome[0]} "
        "{params.disome[1]} {params.gel[0]} {params.gel[1]} {params.libraries} 2> {log}"
