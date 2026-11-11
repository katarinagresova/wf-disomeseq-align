"""Figure: sequencing depth of each library, and the reads left for the alignment after the read filters.

One bar for each library (design order, grouped by replicate): its sequenced reads, split into the reads that
are aligned (footprint of minlength or more) and the reads that the filters of trim.smk remove: no adapter and
adapter dimer (cutadapt_reads), footprint shorter than minlength (trim_reads).

From the logs of trim.smk. cutadapt_reads gives the sequenced reads, no adapter and adapter dimers; the length
histogram of collapse_reads (inserts with the random nt) gives the reads with a footprint shorter than
minlength. Checks (exits with a message if one fails): the cutadapt counts add up; collapse_reads read the reads
that cutadapt_reads wrote, and its histogram has them all; trim_reads read the molecules that collapse_reads
wrote, and the molecules it removed as too short are those of the histogram of distinct reads.

Usage: plot_read_depth.py <design.tsv> <min_insert> <minlength> <out prefix>
       <sample>:<cutadapt.log>:<collapse.log>:<trim.log> ...
"""
import re
import sys

import numpy as np
import pandas as pd

import plotting as P

M = 1e6
CUTADAPT = {"Total reads processed": "processed", "Reads discarded as untrimmed": "no_adapter",
            "Reads that were too short": "too_short", "Reads written (passing filters)": "written"}


def cutadapt_stats(path):
    stats = {"no_adapter": 0, "too_short": 0}
    with open(path) as f:
        for line in f:
            m = re.match(r"(.+?):\s+([\d,]+)", line)
            if m and m.group(1) in CUTADAPT:
                stats[CUTADAPT[m.group(1)]] = int(m.group(2).replace(",", ""))
    return stats


def collapse_stats(path):
    """Reads in, distinct reads out, and the reads and the distinct reads by insert length."""
    stats, lengths, block = {}, {"read_lengths_raw": {}, "read_lengths_unique": {}}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                block = line[1:]
            elif block == "read statistics" and line.startswith(("input\t", "uniq\t")):
                key, value = line.split("\t")
                stats[key] = int(value)
            elif block in lengths and line and not line.startswith("#"):
                length, n = line.split("\t")
                lengths[block][int(length)] = int(n)
    return stats["input"], stats["uniq"], lengths["read_lengths_raw"], lengths["read_lengths_unique"]


def library_counts(sample, cutadapt_log, collapse_log, trim_log, min_insert):
    cut = cutadapt_stats(cutadapt_log)
    if cut["processed"] != cut["no_adapter"] + cut["too_short"] + cut["written"]:
        sys.exit(f"{sample}: the cutadapt_reads counts do not add up")
    if not cut["written"]:
        sys.exit(f"{sample}: no reads passed cutadapt_reads (check the adapter)")
    reads_in, uniq, raw, unique = collapse_stats(collapse_log)
    if reads_in != cut["written"] or sum(raw.values()) != reads_in or sum(unique.values()) != uniq:
        sys.exit(f"{sample}: cutadapt_reads wrote {cut['written']} reads, collapse_reads read {reads_in}")
    trim = cutadapt_stats(trim_log)
    short_molecules = sum(n for length, n in unique.items() if length < min_insert)
    if trim["processed"] != uniq or trim["too_short"] != short_molecules:
        sys.exit(f"{sample}: trim_reads read {trim['processed']} molecules and removed {trim['too_short']}; "
                 f"collapse_reads wrote {uniq}, {short_molecules} of them shorter than {min_insert} nt")
    short = sum(n for length, n in raw.items() if length < min_insert)
    return {"sequenced": cut["processed"], "no_adapter": cut["no_adapter"], "adapter_dimer": cut["too_short"],
            "too_short": short, "aligned": cut["written"] - short}


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    design_tsv, min_insert, minlength, prefix, *libraries = sys.argv[1:]
    design = P.read_design(design_tsv)
    counts = {}
    for library in libraries:
        sample, *logs = library.split(":")
        counts[sample] = library_counts(sample, *logs, int(min_insert))
    if set(design["sample_id"]) != set(counts):
        sys.exit("design.tsv and the logs have different libraries")
    s = pd.DataFrame(counts).T.loc[design["sample_id"]]

    short = f"footprint < {minlength} nt"
    parts = {  # stack order, left to right
        f"aligned (footprint ≥ {minlength} nt)": s["aligned"],
        short: s["too_short"],
        "adapter dimer": s["adapter_dimer"],
        "no adapter": s["no_adapter"],
    }
    colors = [P.CATEGORICAL[0], P.MUTED, "#b0aea6", "#d6d4cc"]
    data = pd.DataFrame({"replicate": design["replicate"].to_numpy(), "condition": design["condition"].to_numpy(),
                         "sequenced": s["sequenced"], "after_adapter_trimming":
                         s["sequenced"] - s["no_adapter"] - s["adapter_dimer"], **parts})
    data.index.name = "sample"
    aligned = data[list(parts)[0]]

    y, labels, pos = P.library_rows(design)
    fig, ax = P.plt.subplots(figsize=(10, 0.42 * pos + 1.6), layout="constrained")
    left = np.zeros(len(s))
    for (name, values), color in zip(parts.items(), colors):
        ax.barh(y, values / M, left=left, height=0.62, color=color, edgecolor="white", linewidth=1.2, label=name)
        left += values.to_numpy() / M
    xmax = data["sequenced"].max() / M
    for yi, total, kept in zip(y, data["sequenced"], aligned):
        ax.text(total / M + xmax * 0.01, yi, f"{total / M:.1f} M", va="center", fontsize=9.5, color=P.INK2)
        text = f"{kept / M:.1f} M ({100 * kept / total:.0f}%)"
        if kept / M > xmax * 0.16:  # fits inside the bar
            ax.text(xmax * 0.008, yi, text, va="center", fontsize=9.5, color=P.LABEL_ON[P.CATEGORICAL[0]])
    P.library_axis(ax, design, y, labels, pos)
    ax.set_xlim(0, xmax * 1.1)
    P.round_ticks(ax.xaxis)
    ax.set_xlabel("reads (millions)")
    handles, names = ax.get_legend_handles_labels()
    ax.legend(handles, names, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, handlelength=1, handleheight=1,
              columnspacing=1.6)
    P.save(fig, prefix, data,
           "Sequencing depth of each library (bar length, millions of reads) and the reads left for the alignment "
           f"(blue: footprint of {minlength} nt or more, after the adapter and the random nt are removed). Grey: "
           f"reads removed by the filters: {short}, adapter dimer (insert of random nt only) and no adapter. "
           "Text: sequenced reads; aligned reads and their share.")


if __name__ == "__main__":
    main()
