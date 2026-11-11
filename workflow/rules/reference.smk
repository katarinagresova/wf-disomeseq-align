# The transcripts to align to (one for each protein-coding gene) with their
# CDS positions, for the main species and the optional spike-in, and the
# ncRNAs of the annotation (contaminants.smk uses those of the main species).


def species_reference(wildcards):
    reference = config if wildcards.species == "main" else config["spike_in"]
    return {"genome": reference["genome_fa"], "gtf": reference["annotation_gtf"]}


rule extract_transcripts:
    input:
        unpack(species_reference),
        script=workflow.source_path("../scripts/extract_transcripts.py"),
    output:
        fa=f"{RESULTS_DIR}/reference/{{species}}/transcripts.fa",
        tsv=f"{RESULTS_DIR}/reference/{{species}}/transcripts.tsv",
        ncrna=f"{RESULTS_DIR}/reference/{{species}}/ncrna.fa",
    # the counts of transcripts and genes left out
    log:
        f"{LOG_DIR}/extract_transcripts/{{species}}.log",
    conda:
        "../envs/pysam.yaml"
    resources:
        mem_mb=8000,
    shell:
        "python {input.script} {input.genome} {input.gtf} {output.fa} {output.tsv} {output.ncrna} 2> {log}"
