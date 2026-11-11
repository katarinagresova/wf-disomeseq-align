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

# Optional spike-in species (config spike_in): its transcripts are aligned to
# with the main ones; its reads get their own BAM.
SPECIES = ["main", "spike_in"] if "spike_in" in config else ["main"]

wildcard_constraints:
    species="main|spike_in",
