# Config and samples.tsv -> the globals the rules use.

import re

import pandas as pd
from snakemake.utils import validate

# also fills in the defaults (RESULTS_DIR, LOG_DIR)
validate(config, "../schemas/config.schema.yaml")

RESULTS_DIR = config["RESULTS_DIR"]
LOG_DIR = config["LOG_DIR"]

samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample_id", drop=False)
validate(samples, "../schemas/samples.schema.yaml")
if samples.index.duplicated().any():
    raise ValueError(f"samples.tsv: duplicate sample_id {sorted(set(samples.index[samples.index.duplicated()]))}")

# {sample} only ever matches a whole sample_id
wildcard_constraints:
    sample="|".join(re.escape(s) for s in samples.index),

def design_table():
    """config design -> one row for each library, in plot order (replicate, then condition, then samples.tsv
    order): sample_id, condition, replicate, the factors, condition_label and condition_order. Without factors,
    the conditions keep the order of samples.tsv and are their own labels."""
    design = config.get("design", {})
    table = samples[["sample_id", "condition"]].reset_index(drop=True)
    replicate = design.get("replicate")
    if replicate and replicate not in samples:
        raise ValueError(f"design replicate: samples.tsv has no column {replicate!r}")
    table["replicate"] = samples[replicate].to_numpy() if replicate else ""
    if table["replicate"].isna().any():
        raise ValueError(f"samples.tsv: empty {replicate!r}")
    # YAML reads a level such as 1 as a number; samples.tsv values are strings
    factors = {f: {str(k): v for k, v in levels.items()} for f, levels in design.get("factors", {}).items()}
    for factor, levels in factors.items():
        if factor not in samples:
            raise ValueError(f"design factor {factor!r}: samples.tsv has no such column")
        unknown = set(samples[factor]) - set(levels)
        if unknown:
            raise ValueError(f"samples.tsv {factor!r}: {sorted(map(str, unknown))} are not levels of the design")
        table[factor] = samples[factor].to_numpy()
    if factors:
        if (table.groupby("condition")[list(factors)].nunique() > 1).to_numpy().any():
            raise ValueError("samples.tsv: a condition has more than one level of a design factor")
        if (table.groupby(list(factors))["condition"].nunique() > 1).any():
            raise ValueError("samples.tsv: two conditions have the same levels of the design factors")
        first = table.drop_duplicates("condition").set_index("condition")
        conditions = sorted(first.index, key=lambda c: [list(factors[f]).index(first.at[c, f]) for f in factors])
        labels = {c: "; ".join(factors[f][first.at[c, f]] for f in factors) for c in conditions}
    else:
        conditions = list(dict.fromkeys(table["condition"]))
        labels = {c: c for c in conditions}
    table["condition_label"] = table["condition"].map(labels)
    table["condition_order"] = table["condition"].map({c: i for i, c in enumerate(conditions)})
    replicates = {r: i for i, r in enumerate(dict.fromkeys(table["replicate"]))}
    table["replicate_order"] = table["replicate"].map(replicates)
    return table.sort_values(["replicate_order", "condition_order"], kind="stable").drop(columns="replicate_order")


DESIGN = design_table()

# Optional spike-in species (config spike_in): its transcripts are aligned to
# with the main ones; its reads get their own BAM.
SPECIES = ["main", "spike_in"] if "spike_in" in config else ["main"]

# Read classes of assign_reads.py: transcript, spike_in, multi, unaligned and
# the contaminant types (build_contaminants.py: the annotation's rRNA, tRNA and
# ncRNA, the types of contaminant_fastas, spikein<type> for the spike-in
# ncRNAs but not its rRNA). config read_classes groups them for the plots.
CLASSES = {"transcript", "multi", "unaligned", "rRNA", "tRNA", "ncRNA", *config["contaminant_fastas"]}
if "spike_in" in config:
    CLASSES |= {"spike_in", "spikeintRNA", "spikeinncRNA"}
READ_CLASSES = config.get("read_classes", {})
_listed = [c for classes in READ_CLASSES.values() for c in classes]
if set(_listed) - CLASSES:
    raise ValueError(f"read_classes: {sorted(set(_listed) - CLASSES)} are not read classes; "
                     f"the classes are {sorted(CLASSES)}")
if len(_listed) != len(set(_listed)):
    raise ValueError("read_classes: a class is in more than one group")

wildcard_constraints:
    species="main|spike_in",
