"""Figure: saturation of the disome footprints on the main transcripts, for each library.

The disome footprints are the reads of the main BAM (class transcript of assign_reads.py: one best alignment, to
a main transcript, one transcript for each gene) with a disome length. Each BAM record is one molecule (a
collapsed read: footprint and random nt); its name has its read count (x<count>, collapse_reads.py).

Each library is subsampled to FRACTIONS of its sequenced reads: each read is kept with that probability (one
random subsample for each fraction, fixed seed; the full depth is exact). For each subsample:
  molecules  the disome molecules with at least one read left
  genes      the genes with min_molecules or more of these molecules

One panel for each replicate and measure (design order), one line for each library, x: the sequenced reads of
the subsample. The colour is the condition, in design order; a second library of the same condition in a
panel is dashed. The dot is the full depth.

Checks (exits with a message if one fails): the reads and the molecules of each BAM are those of class
transcript in its counts.tsv.

Usage: plot_saturation.py <design.tsv> <disome min-max> <min_molecules> <out prefix>
       <sample>:<bam>:<counts.tsv>:<cutadapt.log> ...
"""
import sys

import numpy as np
import pandas as pd

import footprints as F
import plotting as P

M = 1e6
FRACTIONS = np.linspace(0, 1, 21)
SEED = 1


def main():
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    design_tsv, disome, min_molecules, prefix, *libraries = sys.argv[1:]
    disome = tuple(map(int, disome.split("-")))
    min_molecules = int(min_molecules)
    design = P.read_design(design_tsv)
    files = {}
    for library in libraries:
        sample, *paths = library.split(":")
        files[sample] = paths
    if set(design["sample_id"]) != set(files):
        sys.exit("design.tsv and the BAM files have different libraries")

    rng = np.random.default_rng(SEED)
    curves = []
    for r in design.itertuples():
        bam, counts_tsv, cutadapt_log = files[r.sample_id]
        sequenced = P.sequenced_reads(cutadapt_log)
        fp = F.read_bam(bam)
        P.check_footprints(fp, bam, counts_tsv)
        curve = F.saturation(fp[fp["length"].between(*disome)], FRACTIONS, min_molecules, rng)
        # one transcript for each gene: the transcripts are the genes
        curves.append(pd.DataFrame({"sample": r.sample_id, "replicate": r.replicate, "condition": r.condition,
                                    "fraction": curve["fraction"].round(3),
                                    "sequenced_reads": (curve["fraction"] * sequenced).round().astype(int),
                                    "disome_reads": curve["reads"], "molecules": curve["molecules"],
                                    "genes": curve["transcripts"]}))
    data = pd.concat(curves).set_index("sample")

    color, style = P.condition_colors(design), P.line_styles(design)
    replicates = list(dict.fromkeys(design["replicate"]))
    n = len(replicates)
    fig, axes = P.plt.subplots(n, 2, figsize=(11, 2.5 * n + 1.3), sharex=True, sharey="col", squeeze=False,
                               layout="constrained")
    measures = [("molecules", M, "disome molecules\n(millions)"),
                ("genes", 1, f"genes with ≥ {min_molecules}\ndisome molecules")]
    for i, rep in enumerate(replicates):
        for r in design[design["replicate"] == rep].itertuples():
            curve = data.loc[[r.sample_id]]
            for ax, (column, scale, _) in zip(axes[i], measures):
                x, y = curve["sequenced_reads"] / M, curve[column] / scale
                ax.plot(x, y, color=color[r.condition], linestyle=style[r.sample_id], linewidth=1.8)
                ax.plot(x.iloc[-1], y.iloc[-1], "o", color=color[r.condition], markersize=6,
                        markeredgecolor="white", markeredgewidth=1)
        for j, (ax, (_, _, ylabel)) in enumerate(zip(axes[i], measures)):
            if rep:
                ax.set_title(f"Replicate {rep}", loc="left")
            ax.set_ylabel(ylabel, fontsize=P.AXIS_LABEL_SIZE)
            ax.tick_params(labelsize=P.TICK_LABEL_SIZE)
            P.round_ticks(ax.xaxis)
            P.round_ticks(ax.yaxis, 4)
            if i == n - 1:
                ax.set_xlabel("sequenced reads (millions)", fontsize=P.AXIS_LABEL_SIZE)
    for ax in axes.flat:  # after all lines: the shared limits include them all
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
    P.condition_legend(fig, design, color)

    P.save(fig, prefix, data,
           f"Saturation of the disome footprints ({disome[0]}–{disome[1]} nt) on the main transcripts (reads with "
           f"one best alignment, to a main transcript). Each library is subsampled to 5%, 10%, ..., 100% of its "
           f"sequenced reads (each read kept with that probability; one random subsample, seed {SEED}). Left: "
           f"the disome molecules (collapsed reads: footprint and random nt) with a read left. Right: the genes "
           f"with {min_molecules} or more of these molecules. Dot: the full depth.")


if __name__ == "__main__":
    main()
