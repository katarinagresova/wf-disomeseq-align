# Figures (plot_*.py, style in plotting.py, computations on footprints in
# footprints.py), in results/figures/: <name>.pdf, .png, .tsv (the numbers in
# the figure) and .txt (the caption). The libraries are in the order of
# design.tsv (config design, common.smk).

import json
import shlex

FIGURE = (".pdf", ".png", ".tsv", ".txt")


# The design table of common.smk, for the plot scripts. Written again when the
# design changes (params).
rule design_table:
    output:
        f"{RESULTS_DIR}/figures/design.tsv",
    params:
        table=DESIGN.to_csv(sep="\t", index=False),
    run:
        with open(output[0], "w") as f:
            f.write(params.table)


TRIM_LOGS = [
    f"{LOG_DIR}/cutadapt/{{sample}}.log",
    f"{LOG_DIR}/collapse_reads/{{sample}}.log",
    f"{LOG_DIR}/trim_reads/{{sample}}.log",
]


rule plot_read_depth:
    input:
        logs=expand(TRIM_LOGS, sample=DESIGN["sample_id"]),
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_read_depth.py"),
        style=workflow.source_path("../scripts/plotting.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/read_depth", *FIGURE),
    log:
        f"{LOG_DIR}/figures/read_depth.log",
    params:
        # the shortest insert (footprint + random nt) that trim_reads keeps
        min_insert=config["minlength"] + config["randomer_5p"] + config["randomer_3p"],
        minlength=config["minlength"],
        prefix=f"{RESULTS_DIR}/figures/read_depth",
        # <sample>:<file>:<file>... in TRIM_LOGS order
        libraries=" ".join(":".join([s] + expand(TRIM_LOGS, sample=s)) for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {params.min_insert} {params.minlength} {params.prefix} "
        "{params.libraries} 2> {log}"


# The read classes of the aligned reads (counts.tsv), in the groups of config
# read_classes: as % of each library's aligned reads, and as reads.
rule plot_read_classes:
    input:
        counts=expand(f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv", sample=DESIGN["sample_id"]),
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_read_classes.py"),
        style=workflow.source_path("../scripts/plotting.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/read_classes_percent", *FIGURE),
        multiext(f"{RESULTS_DIR}/figures/read_classes_reads", *FIGURE),
    log:
        f"{LOG_DIR}/figures/read_classes.log",
    params:
        groups=shlex.quote(json.dumps(READ_CLASSES)),
        minlength=config["minlength"],
        prefix=f"{RESULTS_DIR}/figures/read_classes",
        libraries=" ".join(f"{s}:{RESULTS_DIR}/bowtie/{s}.counts.tsv" for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {params.groups} {params.minlength} {params.prefix} "
        "{params.libraries} 2> {log}"


def class_label(read_class):
    """The label of the config read_classes group that has only this class, else the class."""
    return next((label for label, classes in READ_CLASSES.items() if classes == [read_class]), read_class)


# Footprint lengths of one read class (plot_lengths.py): figure -> class,
# monosome lengths, disome lengths ("-": not marked). The mRNA figure marks the
# monosome lengths and, as disome lengths, the gel cut.
LENGTH_FIGURES = {
    "rrna_lengths": ("rRNA", "-", "-"),
    "mrna_lengths": ("transcript", "{}-{}".format(*config["monosome_lengths"]), "{}-{}".format(*config["gel_lengths"])),
}


rule plot_lengths:
    input:
        counts=expand(f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv", sample=DESIGN["sample_id"]),
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_lengths.py"),
        style=workflow.source_path("../scripts/plotting.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/{{figure}}", *FIGURE),
    wildcard_constraints:
        figure="|".join(LENGTH_FIGURES),
    log:
        f"{LOG_DIR}/figures/{{figure}}.log",
    params:
        read_class=lambda wildcards: LENGTH_FIGURES[wildcards.figure][0],
        label=lambda wildcards: shlex.quote(class_label(LENGTH_FIGURES[wildcards.figure][0])),
        monosome=lambda wildcards: LENGTH_FIGURES[wildcards.figure][1],
        disome=lambda wildcards: LENGTH_FIGURES[wildcards.figure][2],
        minlength=config["minlength"],
        prefix=lambda wildcards: f"{RESULTS_DIR}/figures/{wildcards.figure}",
        libraries=" ".join(f"{s}:{RESULTS_DIR}/bowtie/{s}.counts.tsv" for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {params.read_class} {params.label} {params.monosome} {params.disome} "
        "{params.minlength} {params.prefix} {params.libraries} 2> {log}"


# Disome yield, step by step from the sequenced reads to the disome molecules
# (gel_lengths) on the main transcripts.
YIELD_FILES = [
    f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv",
    f"{LOG_DIR}/cutadapt/{{sample}}.log",
]


rule plot_yield:
    input:
        files=expand(YIELD_FILES, sample=DESIGN["sample_id"]),
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_yield.py"),
        style=workflow.source_path("../scripts/plotting.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/disome_yield", *FIGURE),
    log:
        f"{LOG_DIR}/figures/disome_yield.log",
    params:
        label=shlex.quote(class_label("transcript")),
        disome="{}-{}".format(*config["gel_lengths"]),
        minlength=config["minlength"],
        prefix=f"{RESULTS_DIR}/figures/disome_yield",
        # <sample>:<file>:<file> in YIELD_FILES order
        libraries=" ".join(":".join([s] + expand(YIELD_FILES, sample=s)) for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {params.label} {params.disome} {params.minlength} {params.prefix} "
        "{params.libraries} 2> {log}"


# Saturation of the disome footprints (gel_lengths) on the main transcripts:
# molecules and covered genes in subsamples of the sequenced reads.
SATURATION_FILES = [
    f"{RESULTS_DIR}/bowtie/main/{{sample}}.bam",
    f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv",
    f"{LOG_DIR}/cutadapt/{{sample}}.log",
]


rule plot_saturation:
    input:
        files=expand(SATURATION_FILES, sample=DESIGN["sample_id"]),
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_saturation.py"),
        style=workflow.source_path("../scripts/plotting.py"),
        footprints=workflow.source_path("../scripts/footprints.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/saturation", *FIGURE),
    log:
        f"{LOG_DIR}/figures/saturation.log",
    params:
        disome="{}-{}".format(*config["gel_lengths"]),
        min_molecules=config["saturation_gene_molecules"],
        prefix=f"{RESULTS_DIR}/figures/saturation",
        # <sample>:<file>:<file>... in SATURATION_FILES order
        libraries=" ".join(":".join([s] + expand(SATURATION_FILES, sample=s)) for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {params.disome} {params.min_molecules} {params.prefix} "
        "{params.libraries} 2> {log}"


# Where the disome footprints (gel_lengths) are on the main transcripts: the
# 5' ends around the start and stop codons, the 5' UTR / CDS / 3' UTR, the
# frame of the 5' ends in the CDS by length, and a QC by length (share, CDS,
# frame 0).
POSITION_FILES = [
    f"{RESULTS_DIR}/bowtie/main/{{sample}}.bam",
    f"{RESULTS_DIR}/bowtie/{{sample}}.counts.tsv",
]


rule plot_positions:
    input:
        files=expand(POSITION_FILES, sample=DESIGN["sample_id"]),
        transcripts=f"{RESULTS_DIR}/reference/main/transcripts.tsv",
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_positions.py"),
        style=workflow.source_path("../scripts/plotting.py"),
        footprints=workflow.source_path("../scripts/footprints.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/disome_metagene", *FIGURE),
        multiext(f"{RESULTS_DIR}/figures/disome_regions", *FIGURE),
        multiext(f"{RESULTS_DIR}/figures/disome_frames", *FIGURE),
        multiext(f"{RESULTS_DIR}/figures/disome_length_qc", *FIGURE),
    log:
        f"{LOG_DIR}/figures/disome_positions.log",
    params:
        disome="{}-{}".format(*config["gel_lengths"]),
        prefix=f"{RESULTS_DIR}/figures/disome",
        # <sample>:<file>:<file> in POSITION_FILES order
        libraries=" ".join(":".join([s] + expand(POSITION_FILES, sample=s)) for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {input.transcripts} {params.disome} {params.prefix} "
        "{params.libraries} 2> {log}"


# The disome footprints (gel_lengths) on the main transcripts around the sites
# of the peptide motifs of config motifs, as observed over expected 5' ends.
rule plot_motifs:
    input:
        files=expand(POSITION_FILES, sample=DESIGN["sample_id"]),
        transcripts=f"{RESULTS_DIR}/reference/main/transcripts.tsv",
        sequences=f"{RESULTS_DIR}/reference/main/transcripts.fa",
        design=f"{RESULTS_DIR}/figures/design.tsv",
        script=workflow.source_path("../scripts/plot_motifs.py"),
        style=workflow.source_path("../scripts/plotting.py"),
        footprints=workflow.source_path("../scripts/footprints.py"),
    output:
        multiext(f"{RESULTS_DIR}/figures/disome_motifs", *FIGURE),
    log:
        f"{LOG_DIR}/figures/disome_motifs.log",
    params:
        motifs=shlex.quote(json.dumps(config["motifs"])),
        disome="{}-{}".format(*config["gel_lengths"]),
        prefix=f"{RESULTS_DIR}/figures/disome_motifs",
        # <sample>:<file>:<file> in POSITION_FILES order
        libraries=" ".join(":".join([s] + expand(POSITION_FILES, sample=s)) for s in DESIGN["sample_id"]),
    conda:
        "../envs/python.yaml"
    shell:
        "python {input.script} {input.design} {input.transcripts} {input.sequences} {params.motifs} "
        "{params.disome} {params.prefix} {params.libraries} 2> {log}"
