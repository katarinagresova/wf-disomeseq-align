"""Figures: the read classes (assign_reads.py) of the reads that were aligned, for each library.

The aligned reads are the reads with a footprint of minlength or more (plot_read_depth.py, blue); their class
is in counts.tsv. The classes are drawn in groups (config read_classes: label -> classes, in stack order);
classes in no group are drawn as "other". Without groups, each class is its own group (transcript, spike_in,
multi, the contaminant types, unaligned). The colours follow the palette order; a group of only unaligned reads
is grey; groups after the eighth coloured one are drawn as "other".

Two figures, one bar for each library (design order, grouped by replicate):
  <prefix>_percent  % of the library's aligned reads (adds up to 100)
  <prefix>_reads    millions of reads

Usage: plot_read_classes.py <design.tsv> <groups.json> <minlength> <out prefix> <sample>:<counts.tsv> ...
"""
import json
import sys

import numpy as np
import pandas as pd

import plotting as P

M = 1e6
OTHER, OTHER_COLOR = "other", "#d6d4cc"


def groups_of(groups, classes):
    """[(label, classes, colour)] in stack order, for the classes in the data."""
    if not groups:
        first = [c for c in ("transcript", "spike_in", "multi") if c in classes]
        rest = sorted(set(classes) - set(first) - {"unaligned"})
        groups = {c: [c] for c in first + rest + (["unaligned"] if "unaligned" in classes else [])}
    out, slot = [], 0
    for label, members in groups.items():
        if members == ["unaligned"]:
            out.append((label, members, P.MUTED))
        elif slot < len(P.CATEGORICAL):
            out.append((label, members, P.CATEGORICAL[slot]))
            slot += 1
        else:
            out.append((OTHER, members, OTHER_COLOR))
    listed = {c for _, members, _ in out for c in members}
    unlisted = sorted(set(classes) - listed)
    if unlisted:
        print(f"classes in no group (drawn as {OTHER}): {unlisted}", file=sys.stderr)
        out.append((OTHER, unlisted, OTHER_COLOR))
    merged = {}  # several "other" entries become one
    for label, members, color in out:
        merged.setdefault(label, ([], color))[0].extend(members)
    return [(label, members, color) for label, (members, color) in merged.items()]


def bars(design, data, groups, xlabel, percent):
    y, labels, height = P.library_rows(design)
    fig, ax = P.plt.subplots(figsize=(10, 0.42 * height + 1.9), layout="constrained")
    left = np.zeros(len(data))
    for label, _, color in groups:
        values = data[label].to_numpy()
        ax.barh(y, values, left=left, height=0.62, color=color, edgecolor="white", linewidth=1.2, label=label)
        if percent:  # share inside the segment, if it fits
            for yi, li, v in zip(y, left, values):
                if v >= 5:
                    ax.text(li + 0.8, yi, f"{v:.0f}%", va="center", fontsize=9,
                            color=P.LABEL_ON.get(color, P.INK))
        left += values
    if percent:
        ax.set_xlim(0, 100)
        P.percent_ticks(ax.xaxis)
    else:
        xmax = left.max()
        for yi, total in zip(y, left):
            ax.text(total + xmax * 0.01, yi, f"{total:.1f} M", va="center", fontsize=9.5, color=P.INK2)
        ax.set_xlim(0, xmax * 1.1)
        P.round_ticks(ax.xaxis)
    ax.set_xlabel(xlabel)
    P.library_axis(ax, design, y, labels, height)
    ncol = 5 if len(groups) > 4 else len(groups)
    handles, names = ax.get_legend_handles_labels()
    ax.legend(P.row_major(handles, ncol), P.row_major(names, ncol), loc="lower left", bbox_to_anchor=(0, 1.0),
              ncol=ncol, handlelength=1, handleheight=1, columnspacing=1.4, fontsize=9.5)
    return fig


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    design_tsv, groups_json, minlength, prefix, *libraries = sys.argv[1:]
    design = P.read_design(design_tsv)
    reads = {}
    for library in libraries:
        sample, counts_tsv = library.split(":")
        reads[sample] = pd.read_csv(counts_tsv, sep="\t").groupby("class")["reads"].sum()
    if set(design["sample_id"]) != set(reads):
        sys.exit("design.tsv and the counts.tsv files have different libraries")
    by_class = pd.DataFrame(reads).T.fillna(0).astype(int).loc[design["sample_id"]]

    groups = groups_of(json.loads(groups_json), list(by_class.columns))
    counts = pd.DataFrame({label: by_class[members].sum(axis=1) for label, members, _ in groups})
    if (counts.sum(axis=1) != by_class.sum(axis=1)).any():
        sys.exit("the groups do not add up to the aligned reads")
    if (by_class.sum(axis=1) == 0).any():
        sys.exit("a library has no aligned reads")
    percent = 100 * counts.div(counts.sum(axis=1), axis=0)
    info = pd.DataFrame({"replicate": design["replicate"].to_numpy(), "condition": design["condition"].to_numpy()},
                        index=counts.index)
    members = "; ".join(f"{label}: {', '.join(m)}" for label, m, _ in groups)

    fig = bars(design, percent, groups, "% of the aligned reads", percent=True)
    P.save(fig, f"{prefix}_percent", info.join(percent.assign(aligned_reads=counts.sum(axis=1))).rename_axis(
        "sample"),
           f"Read classes of the aligned reads (footprint of {minlength} nt or more) of each library, as % of its "
           f"aligned reads. Classes: {members}. Text: shares of 5% or more.")
    fig = bars(design, counts / M, groups, "aligned reads (millions)", percent=False)
    P.save(fig, f"{prefix}_reads", info.join(counts.assign(aligned_reads=counts.sum(axis=1))).rename_axis("sample"),
           f"Read classes of the aligned reads (footprint of {minlength} nt or more) of each library, in millions "
           f"of reads (bar length: all aligned reads). Classes: {members}.")


if __name__ == "__main__":
    main()
