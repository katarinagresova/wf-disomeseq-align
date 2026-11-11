"""Figure: the disome footprints around the sites of peptide motifs in the CDSs, for each library.

The disome footprints are the molecules (collapsed reads: footprint and random nt, so no PCR duplicates) of the
main BAM (class transcript of assign_reads.py: one best alignment, to a main transcript) with a disome length.
The sites of a motif are its matches in the CDSs (footprints.motif_sites): at least MARGINS nt from the start and
the stop codon, so that the window has no 5' ends of the start and stop codon peaks. No P-site offsets are used.

  <prefix>  the 5' ends around the sites (0: the first nt of the first codon of the motif), as observed over
            expected (footprints.motif_counts; 1: the mean density of the CDSs). One panel for each replicate and
            motif (design order, config order), one line for each library: the colour is the condition, a second
            library of the same condition in a panel is dashed.

Usage: plot_motifs.py <design.tsv> <transcripts.tsv> <transcripts.fa> <motifs JSON> <disome min-max> <out prefix>
       <sample>:<bam>:<counts.tsv> ...
"""
import json
import sys

import numpy as np
import pandas as pd
import pysam

import footprints as F
import plotting as P

WINDOW = (-90, 30)
MARGINS = (150, 90)  # nt after the start codon, before the stop codon: the window stays 60 nt from both


def motifs_figure(design, ratio, motifs):
    color, style = P.condition_colors(design), P.line_styles(design)
    replicates = list(dict.fromkeys(design["replicate"]))
    n = len(replicates)
    fig, axes = P.plt.subplots(n, len(motifs), figsize=(2.9 * len(motifs) + 1.5, 2.4 * n + 1.3), sharex=True,
                               sharey=True, squeeze=False, layout="constrained")
    for i, rep in enumerate(replicates):
        for j, (ax, label) in enumerate(zip(axes[i], motifs)):
            ax.axvline(0, color=P.AXIS, linewidth=1, zorder=1)
            ax.axhline(1, color=P.MUTED, linewidth=0.8, linestyle="--", zorder=1)
            for r in design[design["replicate"] == rep].itertuples():
                y = ratio[label][r.sample_id]
                ax.plot(y.index, y.to_numpy(), color=color[r.condition], linestyle=style[r.sample_id],
                        linewidth=1.1)
            ax.set_title(f"Replicate {rep}, {label}" if rep else label, loc="left")
            ax.set_xlim(*WINDOW)
            ax.set_xticks(range(WINDOW[0], WINDOW[1] + 1, 15))
            P.round_ticks(ax.yaxis, 4)
            ax.tick_params(labelsize=P.TICK_LABEL_SIZE)
            if j == 0:
                ax.set_ylabel("5' ends, observed\nover expected", fontsize=P.AXIS_LABEL_SIZE)
            if i == n - 1:
                ax.set_xlabel("5' end to the motif (nt)", fontsize=P.AXIS_LABEL_SIZE)
    for ax in axes.flat:  # after all lines: the shared limits include them all
        ax.set_ylim(bottom=0)
    P.condition_legend(fig, design, color)
    return fig


def main():
    if len(sys.argv) < 8:
        sys.exit(__doc__)
    design_tsv, transcripts_tsv, transcripts_fa, motifs, disome, prefix, *libraries = sys.argv[1:]
    motifs = json.loads(motifs)
    disome = tuple(map(int, disome.split("-")))
    design = P.read_design(design_tsv)
    files = {}
    for library in libraries:
        sample, *paths = library.split(":")
        files[sample] = paths
    if set(design["sample_id"]) != set(files):
        sys.exit("design.tsv and the BAM files have different libraries")
    transcripts = pd.read_csv(transcripts_tsv, sep="\t").set_index("transcript_id")
    with pysam.FastxFile(transcripts_fa) as fa:
        sequences = {r.name: r.sequence for r in fa}
    sites = F.motif_sites(sequences, transcripts, motifs, MARGINS)

    ratio = {label: {} for label in motifs}
    rows = []
    for sample in design["sample_id"]:
        bam, counts_tsv = files[sample]
        fp = F.read_bam(bam)
        P.check_footprints(fp, bam, counts_tsv)
        fp = F.annotate(fp[fp["length"].between(*disome)], transcripts)
        for label in motifs:
            observed, expected = F.motif_counts(fp, transcripts, sites[label], WINDOW)
            ratio[label][sample] = observed / expected if expected > 0 else observed * np.nan
            rows.append({"sample": sample, "motif": label, "pattern": motifs[label], "sites": len(sites[label]),
                         "footprints": int(observed.sum()), "expected": expected, **ratio[label][sample]})

    info = design.set_index("sample_id")[["replicate", "condition"]]
    data = info.join(pd.DataFrame(rows).set_index("sample")).rename_axis("sample")
    source = (f"disome footprints ({disome[0]}–{disome[1]} nt; molecules, so no PCR duplicates) on the main "
              f"transcripts (reads with one best alignment, to a main transcript)")
    listed = "; ".join(f"{label} ({pattern}): {len(sites[label]):,}" for label, pattern in motifs.items())
    P.save(motifs_figure(design, ratio, motifs), prefix, data,
           f"5' ends of the {source} around the sites of peptide motifs in the CDSs (0: the first nt of the first "
           f"codon of the motif), observed over expected: the 5' ends at each position summed over the sites, "
           f"divided by the mean 5' ends per nt of the CDS of each site's transcript summed over the sites (1, "
           f"dashed: the mean density of the CDSs). Sites: the matches of a regular expression on the amino acids "
           f"of the CDS, in the reading frame (overlapping matches are all sites), at least {MARGINS[0]} nt after "
           f"the start codon and {MARGINS[1]} nt before the stop codon. Motifs (pattern): sites: {listed}. No "
           f"P-site offsets: if the A site of the leading ribosome is 45 nt from the 5' end, a peak at -45 puts "
           f"the first codon of the motif in the A site, -42 the second codon, -39 the third codon from the start "
           f"of the motif.")


if __name__ == "__main__":
    main()
