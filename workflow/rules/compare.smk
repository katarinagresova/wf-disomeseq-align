# Compare the libraries and the conditions of samples.tsv (compare.py): the
# metrics table, saturation curves, gene counts and a PDF report.

COMPARE_LIBRARY_FILES = [
    f"{RESULTS_DIR}/qc/{{sample}}.regions.tsv",
    f"{RESULTS_DIR}/qc/{{sample}}.metagene.tsv",
    f"{RESULTS_DIR}/qc/{{sample}}.genes.tsv",
    f"{RESULTS_DIR}/qc/{{sample}}.dups.tsv",
]


rule compare:
    input:
        files=expand(COMPARE_LIBRARY_FILES, sample=samples.index),
        summary=f"{RESULTS_DIR}/qc/summary.tsv",
        transcripts=f"{RESULTS_DIR}/reference/transcripts.tsv",
        script=workflow.source_path("../scripts/compare.py"),
    output:
        multiext(f"{RESULTS_DIR}/compare/", "metrics.tsv", "saturation.tsv", "gene_counts.tsv", "correlation.tsv",
                 "report.pdf"),
    log:
        f"{LOG_DIR}/compare.log",
    params:
        disome=config["disome_lengths"],
        gel=config["gel_lengths"],
        outdir=f"{RESULTS_DIR}/compare",
        # <sample>:<file>:<file>... in COMPARE_LIBRARY_FILES order
        libraries=" ".join(":".join([s] + expand(COMPARE_LIBRARY_FILES, sample=s)) for s in samples.index),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.summary} {input.transcripts} {params.disome[0]} {params.disome[1]} "
        "{params.gel[0]} {params.gel[1]} {params.outdir} {params.libraries} 2> {log}"
