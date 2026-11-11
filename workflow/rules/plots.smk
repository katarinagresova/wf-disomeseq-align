# Figures (plot_*.py, style in plotting.py), in results/figures/: <name>.pdf,
# .png, .tsv (the numbers in the figure) and .txt (the caption). The libraries
# are in the order of design.tsv (config design, common.smk).

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
