"""Figure: the footprint length distribution of the reads of one read class (assign_reads.py), for each library.

One panel for each replicate (design order), one line for each library: the % of its reads of the class at
each footprint length (counts.tsv; reads, not molecules). The colour is the condition, in design order; a
second library of the same condition in a panel is dashed. The x axis ends at the length below which 99.9% of
the reads of the class are (all libraries); the TSV has all lengths.

With the monosome and disome lengths (min-max, nt), these are shaded, and a table next to each panel gives the
reads of each library at disome lengths: as % of its reads of the class, as % of its aligned reads (all classes,
unaligned included: the reads with a footprint of minlength or more, as in plot_read_classes.py) and in millions.

Usage: plot_lengths.py <design.tsv> <class> <label> <monosome min-max | -> <disome min-max | -> <minlength>
       <out prefix> <sample>:<counts.tsv> ...
"""
import sys

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory, offset_copy

import plotting as P

M = 1e6
BAND = "#efeee9"
SHOWN = 0.999  # the x axis shows the lengths of this fraction of the reads


def length_range(arg):
    if arg == "-":
        return None
    low, high = map(int, arg.split("-"))
    if low > high:
        sys.exit(f"length range {arg}: the first length is larger")
    return low, high


def table(ax, fig, libs, stats, color, style, label, disome):
    """The disome reads of the libraries of one panel, in an empty axes next to it."""
    ax.axis("off")
    columns = [(0.58, f"% of\n{label}"), (0.8, "% of\naligned"), (1.0, "million\nreads")]
    ax.text(0, 1, f"Reads at disome lengths ({disome[0]}–{disome[1]} nt)", transform=ax.transAxes, va="top",
            fontsize=9.5, fontweight="bold", color=P.INK)
    for x, text in columns:
        ax.annotate(text, (x, 1), xycoords="axes fraction", xytext=(0, -16), textcoords="offset points", ha="right",
                    va="top", fontsize=8.5, color=P.INK2, linespacing=1.1)
    for i, r in enumerate(libs.itertuples()):
        dy = -50 - 15 * i
        at = offset_copy(ax.transAxes, fig=fig, y=dy, units="points")
        ax.add_line(Line2D([0, 0.09], [1, 1], transform=at, color=color[r.condition], linestyle=style[r.sample_id],
                           linewidth=2, clip_on=False))
        s = stats.loc[r.sample_id]
        values = [f"{s['disome_pct_class']:.0f}%", f"{s['disome_pct_aligned']:.1f}%", f"{s['disome_reads'] / M:.2f}"]
        ax.annotate(f"{r.sample_id}  {r.condition}", (0.12, 1), xycoords="axes fraction", xytext=(0, dy),
                    textcoords="offset points", va="center", fontsize=9, color=P.INK)
        for (x, _), v in zip(columns, values):
            ax.annotate(v, (x, 1), xycoords="axes fraction", xytext=(0, dy), textcoords="offset points", ha="right",
                        va="center", fontsize=9, color=P.INK)


def main():
    if len(sys.argv) < 9:
        sys.exit(__doc__)
    design_tsv, read_class, label, monosome, disome, minlength, prefix, *libraries = sys.argv[1:]
    monosome, disome = length_range(monosome), length_range(disome)
    design = P.read_design(design_tsv)
    counts = {}
    for library in libraries:
        sample, counts_tsv = library.split(":")
        counts[sample] = pd.read_csv(counts_tsv, sep="\t")
    if set(design["sample_id"]) != set(counts):
        sys.exit("design.tsv and the counts.tsv files have different libraries")
    samples = list(design["sample_id"])

    reads = pd.DataFrame({s: counts[s][counts[s]["class"] == read_class].groupby("length")["reads"].sum()
                          for s in samples}).fillna(0)
    if reads.empty or (reads.sum() == 0).any():
        sys.exit(f"a library has no reads of class {read_class}")
    reads = reads.reindex(range(reads.index.min(), reads.index.max() + 1), fill_value=0).astype(int)
    percent = 100 * reads / reads.sum()
    pooled = reads.sum(axis=1).cumsum() / reads.to_numpy().sum()
    xmax = int(np.ceil(pooled.index[pooled.searchsorted(SHOWN)] / 10) * 10)

    stats = pd.DataFrame({"replicate": design["replicate"].to_numpy(), "condition": design["condition"].to_numpy(),
                          "class_reads": reads.sum().to_numpy(),
                          "aligned_reads": [counts[s]["reads"].sum() for s in samples]}, index=samples)
    if disome:
        stats["disome_reads"] = reads.loc[disome[0]:disome[1]].sum()
        stats["disome_pct_class"] = 100 * stats["disome_reads"] / stats["class_reads"]
        stats["disome_pct_aligned"] = 100 * stats["disome_reads"] / stats["aligned_reads"]

    color, style = P.condition_colors(design), P.line_styles(design)
    replicates = list(dict.fromkeys(design["replicate"]))

    n = len(replicates)
    if disome:
        fig, axes = P.plt.subplots(n, 2, figsize=(12, 2.4 * n + 1.2), width_ratios=[2.6, 1], squeeze=False,
                                   layout="constrained")
        plots = axes[:, 0]
        for ax in plots[1:]:
            ax.sharex(plots[0])
            ax.sharey(plots[0])
    else:
        fig, plots = P.plt.subplots(n, 1, figsize=(10, 2.4 * n + 1.2), sharex=True, sharey=True, squeeze=False,
                                    layout="constrained")
        plots = plots[:, 0]
    ymax = percent.loc[:xmax].to_numpy().max() * (1.3 if disome else 1.08)
    for i, (ax, rep) in enumerate(zip(plots, replicates)):
        libs = design[design["replicate"] == rep]
        for name, lengths in (("monosome", monosome), ("disome", disome)):
            if lengths:
                ax.axvspan(lengths[0] - 0.5, lengths[1] + 0.5, color=BAND, zorder=0, linewidth=0)
                ax.text((lengths[0] + lengths[1]) / 2, 0.97, f"{name}\n{lengths[0]}–{lengths[1]} nt",
                        transform=blended_transform_factory(ax.transData, ax.transAxes), ha="center", va="top",
                        fontsize=8.5, color=P.INK2, linespacing=1.1)
        for r in libs.itertuples():
            ax.plot(percent.index, percent[r.sample_id], color=color[r.condition], linestyle=style[r.sample_id],
                    linewidth=1.8)
        if rep:
            ax.set_title(f"Replicate {rep}", loc="left")
        ax.set_xlim(int(minlength) - 0.5, xmax + 0.5)
        ax.set_ylim(0, ymax)
        P.round_ticks(ax.xaxis, 12)
        P.percent_ticks(ax.yaxis, 4)
        ax.set_ylabel(f"% of {label} reads")
        if i < n - 1:
            ax.tick_params(axis="x", labelbottom=False)
        if disome:
            table(axes[i, 1], fig, libs, stats, color, style, label, disome)
    plots[-1].set_xlabel("footprint length (nt)")
    P.condition_legend(fig, design, color)

    data = stats.join(percent.T.rename(columns=str)).rename_axis("sample")
    ranges = []
    if monosome:
        ranges.append(f"monosome lengths ({monosome[0]}–{monosome[1]} nt)")
    if disome:
        ranges.append(f"disome lengths ({disome[0]}–{disome[1]} nt)")
    caption = (f"Footprint lengths of the {label} reads (class {read_class}) of each library, as % of its {label} "
               f"reads at each length; one panel for each replicate. The x axis ends at {xmax} nt "
               f"({SHOWN:.1%} of the reads); the TSV has all lengths.")
    if ranges:
        caption += f" Shaded: {' and '.join(ranges)}."
    if disome:
        caption += (f" Table: the reads at disome lengths as % of the {label} reads, as % of the aligned reads "
                    f"(all reads with a footprint of {minlength} nt or more) and in millions.")
    P.save(fig, prefix, data, caption)


if __name__ == "__main__":
    main()
