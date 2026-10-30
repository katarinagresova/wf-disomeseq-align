# The contaminant set: the ncRNAs of the annotation (extract_transcripts) and
# the fastas of config contaminant_fastas, by type. bowtie_index indexes it
# with the transcripts.


rule build_contaminants:
    input:
        ncrna=f"{RESULTS_DIR}/reference/ncrna.fa",
        fastas=[path for paths in config["contaminant_fastas"].values() for path in paths],
        script=workflow.source_path("../scripts/build_contaminants.py"),
    output:
        f"{RESULTS_DIR}/reference/contaminants.fa",
    log:
        f"{LOG_DIR}/build_contaminants.log",
    params:
        typed=" ".join(f"{rtype}={path}" for rtype, paths in config["contaminant_fastas"].items() for path in paths),
    conda:
        "../envs/pysam.yaml"
    shell:
        "python {input.script} {output} {input.ncrna} {params.typed} 2> {log}"
