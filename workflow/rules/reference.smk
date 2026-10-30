# The transcripts to align to (one for each protein-coding gene) with their
# CDS positions, and the ncRNAs of the annotation (contaminants.smk).


rule extract_transcripts:
    input:
        genome=config["genome_fa"],
        gtf=config["annotation_gtf"],
        script=workflow.source_path("../scripts/extract_transcripts.py"),
    output:
        fa=f"{RESULTS_DIR}/reference/transcripts.fa",
        tsv=f"{RESULTS_DIR}/reference/transcripts.tsv",
        ncrna=f"{RESULTS_DIR}/reference/ncrna.fa",
    # the counts of transcripts and genes left out
    log:
        f"{LOG_DIR}/extract_transcripts.log",
    conda:
        "../envs/pysam.yaml"
    resources:
        mem_mb=8000,
    shell:
        "python {input.script} {input.genome} {input.gtf} {output.fa} {output.tsv} {output.ncrna} 2> {log}"
