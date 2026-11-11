# The transcripts to align to (one for each protein-coding gene, without
# retrocopies) with their CDS positions, for the main species and the optional
# spike-in, and the ncRNAs of the annotation (contaminants.smk uses those of
# the main species).


def species_reference(wildcards):
    reference = config if wildcards.species == "main" else config["spike_in"]
    return {"genome": reference["genome_fa"], "gtf": reference["annotation_gtf"]}


rule extract_transcripts:
    input:
        unpack(species_reference),
        script=workflow.source_path("../scripts/extract_transcripts.py"),
    output:
        fa=f"{RESULTS_DIR}/reference/{{species}}/all_transcripts.fa",
        tsv=f"{RESULTS_DIR}/reference/{{species}}/all_transcripts.tsv",
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


# Retrocopies (processed pseudogenes with a protein-coding model) take the
# footprints of their parent gene as multi-mapping reads: remove them
# (retrocopies.py). The CDS tiles of the transcripts with few exons are aligned
# to the transcripts with many exons. The spike-in can have no such pairs.
rule retrocopies:
    input:
        fa=f"{RESULTS_DIR}/reference/{{species}}/all_transcripts.fa",
        tsv=f"{RESULTS_DIR}/reference/{{species}}/all_transcripts.tsv",
        script=workflow.source_path("../scripts/retrocopies.py"),
    output:
        fa=f"{RESULTS_DIR}/reference/{{species}}/transcripts.fa",
        tsv=f"{RESULTS_DIR}/reference/{{species}}/transcripts.tsv",
        retrocopies=f"{RESULTS_DIR}/reference/{{species}}/retrocopies.tsv",
    log:
        f"{LOG_DIR}/retrocopies/{{species}}.log",
    params:
        tmp=f"{RESULTS_DIR}/reference/{{species}}/retrocopies_tmp",
    conda:
        "../envs/bowtie.yaml"
    threads: 8
    shell:
        r"""
        rm -rf {params.tmp} && mkdir -p {params.tmp}
        python {input.script} tiles {input.fa} {input.tsv} {params.tmp}/partners.fa {params.tmp}/tiles.fa 2> {log}
        if [ -s {params.tmp}/partners.fa ] && [ -s {params.tmp}/tiles.fa ]; then
            bowtie-build --threads {threads} -q {params.tmp}/partners.fa {params.tmp}/partners >> {log} 2>&1
            # tile name and partner transcript of the first alignment (up to 2 mismatches)
            bowtie -p {threads} -f -v 2 --norc -k 1 --suppress 2,4,5,6,7,8 {params.tmp}/partners \
                {params.tmp}/tiles.fa > {params.tmp}/hits.txt 2>> {log}
        else
            : > {params.tmp}/hits.txt
        fi
        python {input.script} filter {input.fa} {input.tsv} {params.tmp}/hits.txt {output.fa} {output.tsv} \
            {output.retrocopies} 2>> {log}
        rm -rf {params.tmp}
        """
