# The contaminant set: the ncRNAs of the annotation (extract_transcripts), the
# fastas of config contaminant_fastas, by type, and the ncRNAs of the spike-in
# annotation (type spikein<type>, not its rRNA). bowtie_index indexes it with
# the transcripts.


rule build_contaminants:
    input:
        ncrna=f"{RESULTS_DIR}/reference/main/ncrna.fa",
        spike_in_ncrna=expand(f"{RESULTS_DIR}/reference/{{species}}/ncrna.fa", species=SPECIES[1:]),
        fastas=[path for paths in config["contaminant_fastas"].values() for path in paths],
        script=workflow.source_path("../scripts/build_contaminants.py"),
    output:
        f"{RESULTS_DIR}/reference/contaminants.fa",
    log:
        f"{LOG_DIR}/build_contaminants.log",
    params:
        typed=" ".join(f"{rtype}={path}" for rtype, paths in config["contaminant_fastas"].items() for path in paths),
        spike_in=lambda wildcards, input: " ".join(f"spike_in={path}" for path in input.spike_in_ncrna),
    conda:
        "../envs/pysam.yaml"
    shell:
        "python {input.script} {output} {input.ncrna} {params.typed} {params.spike_in} 2> {log}"
