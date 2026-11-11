"""Figure: the disome yield of each library, step by step from the sequenced reads to the disome molecules.

One panel for each step, one bar for each library (design order); the colour is the condition.
  sequenced   the sequenced reads (cutadapt: Total reads processed), in millions
  minlength   the reads with a footprint of minlength nt or more (the aligned reads of plot_read_classes.py: all
              classes, unaligned included), % of the sequenced reads
  mRNA        the reads of class transcript (assign_reads.py: one best alignment, to a main transcript), % of
              the reads of the step before
  disome      the mRNA reads at disome lengths, % of the mRNA reads
  molecules   the molecules of these reads (collapsed reads: footprint and random nt, so no PCR duplicates), %
              of these reads
  yield       the disome molecules per 100 sequenced reads: the product of the four steps

The reads and molecules are from counts.tsv (each class and length). Checks (exits with a message if one fails):
each library has reads of class transcript, and not more reads in counts.tsv than sequenced reads.

Usage: plot_yield.py <design.tsv> <mRNA label> <disome min-max> <minlength> <out prefix>
       <sample>:<counts.tsv>:<cutadapt.log> ...
"""
import sys

import pandas as pd

import plotting as P

M = 1e6


def percent(v):
    return f"{v:.1f}%" if v < 10 else f"{v:.0f}%"


def library_counts(sample, counts_tsv, cutadapt_log, disome):
    counts = pd.read_csv(counts_tsv, sep="\t")
    mrna = counts[counts["class"] == "transcript"]
    at_disome = mrna[mrna["length"].between(*disome)]
    row = {"sequenced": P.sequenced_reads(cutadapt_log), "min_reads": counts["reads"].sum(),
           "class_reads": mrna["reads"].sum(), "disome_reads": at_disome["reads"].sum(),
           "disome_molecules": at_disome["molecules"].sum()}
    if row["class_reads"] == 0:
        sys.exit(f"{sample}: no reads of class transcript in {counts_tsv}")
    if row["min_reads"] > row["sequenced"]:
        sys.exit(f"{sample}: more reads in {counts_tsv} ({row['min_reads']}) than sequenced ({row['sequenced']})")
    return row


def main():
    if len(sys.argv) < 7:
        sys.exit(__doc__)
    design_tsv, label, disome, minlength, prefix, *libraries = sys.argv[1:]
    disome = tuple(map(int, disome.split("-")))
    design = P.read_design(design_tsv)
    rows = {}
    for library in libraries:
        sample, counts_tsv, cutadapt_log = library.split(":")
        rows[sample] = library_counts(sample, counts_tsv, cutadapt_log, disome)
    if set(design["sample_id"]) != set(rows):
        sys.exit("design.tsv and the counts.tsv files have different libraries")
    s = pd.DataFrame(rows).T.loc[design["sample_id"]]

    lengths = f"{disome[0]}–{disome[1]} nt"
    steps = [  # TSV column, panel title, % of the step before
        ("pct_min_of_sequenced", f"footprint ≥ {minlength} nt\n% of sequenced", s["min_reads"] / s["sequenced"]),
        ("pct_class_of_min", f"{label}\n% of ≥ {minlength} nt", s["class_reads"] / s["min_reads"]),
        ("pct_disome_of_class", f"{lengths}\n% of {label}", s["disome_reads"] / s["class_reads"]),
        ("pct_molecules_of_disome", "molecules (no PCR\nduplicates), % of reads",
         s["disome_molecules"] / s["disome_reads"]),
    ]
    steps = [(column, title, 100 * v.astype(float)) for column, title, v in steps]
    yield_ = 100 * (s["disome_molecules"] / s["sequenced"]).astype(float)
    panels = ([("sequenced reads\n(millions)", s["sequenced"] / M, [f"{v / M:.1f}" for v in s["sequenced"]])]
              + [(title, v, [percent(x) for x in v]) for _, title, v in steps]
              + [("disome molecules per\n100 sequenced reads", yield_,
                  [f"{v:.1f}  ({n / 1e3:.0f}k)" for v, n in zip(yield_, s["disome_molecules"])])])

    color = P.condition_colors(design)
    colors = [color[c] for c in design["condition"]]
    y, labels, height = P.library_rows(design)
    fig, axes = P.plt.subplots(1, len(panels), figsize=(13.5, 0.42 * height + 1.9), sharey=True,
                               width_ratios=[1] * (len(panels) - 1) + [1.35], layout="constrained")
    for i, (ax, (title, values, texts)) in enumerate(zip(axes, panels)):
        ax.barh(y, values.to_numpy(dtype=float), height=0.62, color=colors)
        for yi, v, text in zip(y, values, texts):
            ax.annotate(text, (v, yi), xytext=(4, 0), textcoords="offset points", va="center", fontsize=9.5,
                        color=P.INK)
        ax.set_xlim(0, values.max() * (1.9 if i == len(panels) - 1 else 1.45))
        ax.set_title(title, loc="left", fontsize=10.5, pad=14)  # pad: room for the first replicate label
        ax.set_xticks([])
        ax.grid(False)
        ax.spines["bottom"].set_visible(False)
        if i:
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", length=0, labelleft=False)
    P.library_axis(axes[0], design, y, labels, height)
    P.condition_legend(fig, design, color)

    data = (design.set_index("sample_id")[["replicate", "condition"]].join(s)
            .assign(**{column: v for column, _, v in steps}, disome_molecules_per_100_sequenced=yield_)
            .rename_axis("sample"))
    P.save(fig, prefix, data,
           f"Disome yield of each library, step by step: the reads with a footprint of {minlength} nt or more (% of "
           f"the sequenced reads); of these, the {label} reads (one best alignment, to a main transcript); of "
           f"these, the reads at disome lengths ({lengths}); of these, the molecules (collapsed reads: footprint "
           f"and random nt, so no PCR duplicates). Right: the product of the four steps, the disome molecules per "
           f"100 sequenced reads (in brackets: the disome molecules).")


if __name__ == "__main__":
    main()
