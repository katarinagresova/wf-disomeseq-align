"""One QC row for each library, from the logs and tables of the read processing and alignment.

Read fates (reads; they add up to the sequenced reads):
  no_adapter     no adapter: not a footprint (cutadapt_reads)
  adapter_dimer  insert of random nt only (cutadapt_reads)
  too_short      footprint shorter than minlength (trim_reads)
  <type>, multi, unaligned, transcript: the classes of assign_reads.py
Molecules (collapsed reads, i.e. without the PCR duplicates):
  molecules, duplicate_fraction: after cutadapt_reads (collapse_reads)
  transcript_molecules, then by the region of the approximate A-site (qc_library.py): utr5, cds, utr3
  usable: cds, disome lengths; gel_cds: cds, gel lengths; monosome_cds: cds, MONOSOME_LENGTHS

Writes summary.tsv (all columns) and read_fates.tsv (sample and the read fates, for MultiQC).

Usage: qc_summary.py <summary.tsv> <read_fates.tsv> <samples.tsv> <disome_min> <disome_max> <gel_min>
       <gel_max> <sample>:<cutadapt.log>:<collapse.log>:<counts.tsv>:<regions.tsv> ...
"""
import re
import sys

import pandas as pd

MONOSOME_LENGTHS = (26, 34)
CUTADAPT = {"Total reads processed": "sequenced", "Reads discarded as untrimmed": "no_adapter",
            "Reads that were too short": "adapter_dimer", "Reads written (passing filters)": "written"}


def cutadapt_stats(path):
    stats = {}
    with open(path) as f:
        for line in f:
            m = re.match(r"(.+?):\s+([\d,]+)", line)
            if m and m.group(1) in CUTADAPT:
                stats[CUTADAPT[m.group(1)]] = int(m.group(2).replace(",", ""))
    stats.setdefault("no_adapter", 0)
    return stats


def collapse_stats(path):
    with open(path) as f:
        stats = dict(line.rstrip("\n").split("\t") for line in f if line.startswith(("input\t", "uniq\t")))
    return int(stats["input"]), int(stats["uniq"])


def in_range(lengths, lo_hi):
    return (lengths >= lo_hi[0]) & (lengths <= lo_hi[1])


def main():
    if len(sys.argv) < 9:
        sys.exit(__doc__)
    out_summary, out_fates, samples_tsv, dmin, dmax, gmin, gmax, *libraries = sys.argv[1:]
    disome, gel = (int(dmin), int(dmax)), (int(gmin), int(gmax))
    condition = pd.read_csv(samples_tsv, sep="\t", dtype=str).set_index("sample_id")["condition"]

    rows, classes = [], set()
    for library in libraries:
        sample, cutadapt_log, collapse_log, counts_tsv, regions_tsv = library.split(":")
        row = {"sample": sample, "condition": condition[sample], **cutadapt_stats(cutadapt_log)}
        collapsed_input, uniq = collapse_stats(collapse_log)
        if collapsed_input != row["written"]:
            sys.exit(f"{sample}: cutadapt wrote {row['written']} reads, collapse_reads read {collapsed_input}")

        counts = pd.read_csv(counts_tsv, sep="\t")
        by_class = counts.groupby("class")[["molecules", "reads"]].sum()
        classes |= set(by_class.index)
        row.update(by_class["reads"])
        row["too_short"] = row.pop("written") - by_class["reads"].sum()

        row["molecules"] = uniq
        row["duplicate_fraction"] = 1 - uniq / collapsed_input
        regions = pd.read_csv(regions_tsv, sep="\t")
        row["transcript_molecules"] = regions["molecules"].sum()
        for region in ("utr5", "cds", "utr3"):
            row[f"{region}_molecules"] = regions.loc[regions["region"] == region, "molecules"].sum()
        cds = regions[regions["region"] == "cds"]
        for name, lo_hi in (("usable", disome), ("gel_cds", gel), ("monosome_cds", MONOSOME_LENGTHS)):
            row[name] = cds.loc[in_range(cds["length"], lo_hi), "molecules"].sum()
        rows.append(row)

    fates = (["no_adapter", "adapter_dimer", "too_short"]
             + sorted(classes - {"multi", "unaligned", "transcript"}) + ["multi", "unaligned", "transcript"])
    summary = pd.DataFrame(rows).fillna({c: 0 for c in fates})
    if (summary[fates].sum(axis=1) != summary["sequenced"]).any():
        sys.exit("the read fates do not add up to the sequenced reads")
    first = ["sample", "condition", "sequenced"] + fates
    summary = summary[first + [c for c in summary.columns if c not in first]]
    summary.to_csv(out_summary, sep="\t", index=False)
    summary[["sample"] + fates].to_csv(out_fates, sep="\t", index=False)


if __name__ == "__main__":
    main()
