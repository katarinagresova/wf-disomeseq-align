"""Figures: where the disome footprints are on the main transcripts, for each library.

The disome footprints are the molecules (collapsed reads: footprint and random nt, so no PCR duplicates) of the
main BAM (class transcript of assign_reads.py: one best alignment, to a main transcript) with a disome length.
Positions are on the transcript (transcripts.tsv: cds_start, cds_end; the stop codon is in the CDS). No
P-site offsets are used.

  <prefix>_metagene  the 5' ends around the start codon (0: the A of the ATG) and the stop codon (0: the first
                     nt of the stop codon), per million disome footprints of the library. One panel for each
                     replicate and codon (design order), one line for each library: the colour is the condition,
                     a second library of the same condition in a panel is dashed.
  <prefix>_regions   the % of the disome footprints in the 5' UTR, the CDS and the 3' UTR, from the footprint
                     midpoint (5' end + length // 2). One bar for each library (design order).
  <prefix>_frames    the frame of the 5' end (0: the first nt of a codon, from the start codon) of the disome
                     footprints with the midpoint in the CDS: for each length, the % of its footprints in each
                     frame. One heatmap for each library (design order).
  <prefix>_length_qc by footprint length: the % of the disome footprints, the % of the footprints of the length in
                     the CDS, and the % of these in frame 0. One panel for each measure and replicate, one line for
                     each library.

Usage: plot_positions.py <design.tsv> <transcripts.tsv> <disome min-max> <out prefix>
       <sample>:<bam>:<counts.tsv> ...
"""
import sys

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize

import footprints as F
import plotting as P

START_WINDOW = (-60, 150)
STOP_WINDOW = (-150, 60)
SEQUENTIAL = LinearSegmentedColormap.from_list("blue", P.SEQUENTIAL)


def metagene_figure(design, metagene):
    color, style = P.condition_colors(design), P.line_styles(design)
    replicates = list(dict.fromkeys(design["replicate"]))
    n = len(replicates)
    fig, axes = P.plt.subplots(n, 2, figsize=(11, 2.4 * n + 1.3), sharex="col", sharey=True, squeeze=False,
                               layout="constrained")
    for i, rep in enumerate(replicates):
        for ax, codon in zip(axes[i], ("start", "stop")):
            ax.axvline(0, color=P.AXIS, linewidth=1, zorder=1)
            for r in design[design["replicate"] == rep].itertuples():
                y = metagene[codon][r.sample_id]
                ax.plot(y.index, y.to_numpy(), color=color[r.condition], linestyle=style[r.sample_id],
                        linewidth=1.1)
            ax.set_title(f"Replicate {rep}, {codon} codon" if rep else f"{codon.capitalize()} codon", loc="left")
            ax.set_xlim(*(START_WINDOW if codon == "start" else STOP_WINDOW))
            P.round_ticks(ax.xaxis, 10)
            P.round_ticks(ax.yaxis, 4)
            ax.tick_params(labelsize=P.TICK_LABEL_SIZE)
            if codon == "start":
                ax.set_ylabel("5' ends per million\ndisome footprints", fontsize=P.AXIS_LABEL_SIZE)
            if i == n - 1:
                ax.set_xlabel(f"5' end to the {codon} codon (nt)", fontsize=P.AXIS_LABEL_SIZE)
    for ax in axes.flat:  # after all lines: the shared limits include them all
        ax.set_ylim(bottom=0)
    P.condition_legend(fig, design, color)
    return fig


def frames_figure(design, percent, lengths):
    """One heatmap for each library (design order): frames (rows) by length (columns), % in each cell."""
    replicates = list(dict.fromkeys(design["replicate"]))
    fig, axes = P.plt.subplots(len(design), 1, figsize=(11, 0.82 * len(design) + 0.35 * len(replicates) + 1.2),
                               sharex=True, squeeze=False, layout="constrained")
    axes = axes[:, 0]
    norm = Normalize(0, 100)
    for i, (ax, r) in enumerate(zip(axes, design.itertuples())):
        values = percent[r.sample_id].T.to_numpy()  # frames x lengths
        image = ax.imshow(values, cmap=SEQUENTIAL, norm=norm, aspect="auto", interpolation="none",
                          extent=(lengths[0] - 0.5, lengths[1] + 0.5, 2.5, -0.5))
        for (frame, j), v in np.ndenumerate(values):
            if np.isfinite(v):
                dark = sum(SEQUENTIAL(norm(v))[:3] * np.array([0.2126, 0.7152, 0.0722])) < 0.45
                ax.text(lengths[0] + j, frame, f"{v:.0f}", ha="center", va="center", fontsize=8.5,
                        color="white" if dark else P.INK)
        ax.set_yticks(range(3), [f"frame {f}" for f in range(3)])
        ax.set_ylabel(f"{r.sample_id}\n{r.condition}", rotation=0, ha="right", va="center", color=P.INK,
                      fontsize=10.5)
        ax.tick_params(length=0)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if r.replicate and (i == 0 or r.replicate != design["replicate"].iloc[i - 1]):
            ax.set_title(f"Replicate {r.replicate}", loc="left")
    axes[-1].set_xticks(range(lengths[0], lengths[1] + 1))
    axes[-1].tick_params(axis="x", labelsize=P.TICK_LABEL_SIZE)
    axes[-1].set_xlabel("footprint length (nt)", fontsize=P.AXIS_LABEL_SIZE)
    bar = fig.colorbar(image, ax=axes, shrink=0.4, aspect=25, pad=0.02, anchor=(0, 1))
    bar.set_label("% of the footprints of the length", fontsize=10.5)
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=0)
    P.percent_ticks(bar.ax.yaxis, 5)
    return fig


LENGTH_ROWS = [("% of footprints", "% of the disome\nfootprints"), ("% CDS", "% of the length\nin the CDS"),
               ("% frame 0", "% of the length\nin frame 0")]


def length_qc_figure(design, qc):
    """One row for each measure of LENGTH_ROWS by footprint length, one column for each replicate (design
    order), one line for each library."""
    color, style = P.condition_colors(design), P.line_styles(design)
    replicates = list(dict.fromkeys(design["replicate"]))
    fig, axes = P.plt.subplots(len(LENGTH_ROWS), len(replicates), figsize=(4 * len(replicates) + 1.5, 8.2),
                               sharex=True, sharey="row", squeeze=False, layout="constrained")
    for j, rep in enumerate(replicates):
        for r in design[design["replicate"] == rep].itertuples():
            lib = qc.loc[[r.sample_id]]
            for ax, (column, _) in zip(axes[:, j], LENGTH_ROWS):
                ax.plot(lib["length"], lib[column], color=color[r.condition], linestyle=style[r.sample_id],
                        linewidth=1.6, marker="o", markersize=4)
        axes[-1, j].axhline(100 / 3, color=P.MUTED, linestyle=":", linewidth=1.2, zorder=1)
        if rep:
            axes[0, j].set_title(f"Replicate {rep}", loc="left")
        axes[-1, j].set_xlabel("footprint length (nt)", fontsize=P.AXIS_LABEL_SIZE)
    axes[-1, 0].annotate("no periodicity", (axes[-1, 0].get_xlim()[0], 100 / 3), xytext=(2, -3),
                         textcoords="offset points", va="top", fontsize=9, color=P.INK2)
    for (_, ylabel), ax in zip(LENGTH_ROWS, axes[:, 0]):
        ax.set_ylabel(ylabel, fontsize=P.AXIS_LABEL_SIZE)
    for ax in axes.flat:
        ax.tick_params(labelsize=P.TICK_LABEL_SIZE)
        P.round_ticks(ax.xaxis, 8)
        P.percent_ticks(ax.yaxis, 4)
    # after all lines: the shared limits include them all
    axes[0, 0].set_ylim(bottom=0)
    axes[1, 0].set_ylim(top=100)
    axes[2, 0].set_ylim(bottom=0)
    P.condition_legend(fig, design, color)
    return fig


def regions_figure(design, percent):
    y, labels, height = P.library_rows(design)
    fig, ax = P.plt.subplots(figsize=(10, 0.42 * height + 1.6), layout="constrained")
    left = np.zeros(len(percent))
    for slot, region in enumerate(F.REGIONS):
        values = percent[region].to_numpy()
        color = P.CATEGORICAL[slot]
        ax.barh(y, values, left=left, height=0.62, color=color, edgecolor="white", linewidth=1.2, label=region)
        for yi, li, v in zip(y, left, values):
            if v >= 5:
                ax.text(li + 0.8, yi, f"{v:.0f}%", va="center", fontsize=9, color=P.LABEL_ON[color])
        left += values
    ax.set_xlim(0, 100)
    P.percent_ticks(ax.xaxis)
    ax.set_xlabel("% of the disome footprints")
    P.library_axis(ax, design, y, labels, height)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, handlelength=1, handleheight=1, fontsize=9.5)
    return fig


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    design_tsv, transcripts_tsv, disome, prefix, *libraries = sys.argv[1:]
    disome = tuple(map(int, disome.split("-")))
    design = P.read_design(design_tsv)
    files = {}
    for library in libraries:
        sample, *paths = library.split(":")
        files[sample] = paths
    if set(design["sample_id"]) != set(files):
        sys.exit("design.tsv and the BAM files have different libraries")
    transcripts = pd.read_csv(transcripts_tsv, sep="\t").set_index("transcript_id")

    metagene = {"start": {}, "stop": {}}
    counts, frames, by_length = {}, {}, {}
    for sample in design["sample_id"]:
        bam, counts_tsv = files[sample]
        fp = F.read_bam(bam)
        P.check_footprints(fp, bam, counts_tsv)
        fp = F.annotate(fp[fp["length"].between(*disome)], transcripts)
        per_million = 1e6 / len(fp)
        metagene["start"][sample] = F.metagene(fp["to_start"], START_WINDOW) * per_million
        metagene["stop"][sample] = F.metagene(fp["to_stop"], STOP_WINDOW) * per_million
        counts[sample] = fp["region"].value_counts(sort=False)  # all regions, in REGIONS order
        frames[sample] = F.frames(fp[fp["region"] == "CDS"], disome)
        by_length[sample] = F.regions(fp, disome)

    counts = pd.DataFrame(counts, index=F.REGIONS).T
    percent = 100 * counts.div(counts.sum(axis=1), axis=0)
    info = design.set_index("sample_id")[["replicate", "condition"]]
    lengths = f"{disome[0]}–{disome[1]} nt"
    source = (f"disome footprints ({lengths}; molecules, so no PCR duplicates) on the main transcripts (reads "
              f"with one best alignment, to a main transcript)")

    data = pd.concat({codon: pd.DataFrame(metagene[codon]) for codon in ("start", "stop")}, names=["codon",
                     "position"])
    data = data.T.rename_axis("sample")
    data.columns = [f"{codon}_{position}" for codon, position in data.columns]
    P.save(metagene_figure(design, metagene), f"{prefix}_metagene", info.join(data).rename_axis("sample"),
           f"5' ends of the {source} around the start codon (0: the A of the ATG) and the stop codon (0: the "
           f"first nt of the stop codon), per million disome footprints of the library. No P-site offsets.")
    P.save(regions_figure(design, percent), f"{prefix}_regions",
           info.join(percent).join(counts.add_suffix(" footprints")).rename_axis("sample"),
           f"Region of the {source} of each library, as % of its disome footprints: the region of the footprint "
           f"midpoint (5' end + length / 2). Text: shares of 5% or more.")

    frames = pd.concat(frames, names=["sample"])  # (sample, length) x frame: footprints
    total = frames.sum(axis=1)
    frame_percent = 100 * frames.div(total, axis=0)
    frame_table = frame_percent.add_prefix("% frame ").assign(footprints=total).reset_index("length")
    P.save(frames_figure(design, {s: frame_percent.loc[s] for s in design["sample_id"]}, disome), f"{prefix}_frames",
           info.join(frame_table).rename_axis("sample"),
           f"Frame of the 5' end of the {source}, with the midpoint in the CDS: for each library and footprint "
           f"length, the % of its footprints with the 5' end at the first (frame 0), second (frame 1) or third nt "
           f"(frame 2) of a codon, from the start codon. No periodicity: 33% in each frame. No P-site offsets.")

    by_length = pd.concat(by_length, names=["sample"])  # (sample, length) x region: footprints
    n = by_length.sum(axis=1)
    qc = (pd.DataFrame({"footprints": n, "% of footprints": 100 * n / n.groupby(level="sample").transform("sum")})
          .join(100 * by_length.div(n, axis=0).add_prefix("% "))
          .join(frame_percent.add_prefix("% frame "))
          .reset_index("length"))
    P.save(length_qc_figure(design, qc), f"{prefix}_length_qc", info.join(qc).rename_axis("sample"),
           f"The {source}, by footprint length, for each library. Top: the % of its disome footprints at each "
           f"length. Middle: the % of the footprints of that length with the midpoint in the CDS. Bottom: the % of "
           f"the footprints of that length with the midpoint in the CDS that have the 5' end in frame 0 (the first "
           f"nt of a codon, from the start codon); dotted: 33%, no periodicity. The TSV has all regions and "
           f"frames. No P-site offsets.")


if __name__ == "__main__":
    main()
