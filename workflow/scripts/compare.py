"""Compare the libraries and their conditions (library preparations): metrics, saturation, gene counts, report.

Inputs: qc_summary.py's summary.tsv, the transcripts.tsv of extract_transcripts.py and, for each library,
the tables of qc_library.py. Usable footprints: disome length, approximate A-site in the CDS, collapsed.

metrics.tsv, one row for each library:
  yield                 usable / sequenced reads
  usable_at_<d>M        expected usable footprints at the same depth d (the smallest library, in M reads)
  usable_duplicates     duplicate fraction of the usable footprints (1 - molecules / reads)
  <fate>_fraction       the read fates of summary.tsv (all but transcript), as fractions of the sequenced reads
  gel_fraction          fraction of the CDS footprints (all lengths) in the gel lengths
  monosome_fraction     fraction of the CDS footprints in 26-34 nt (qc_summary.py MONOSOME_LENGTHS)
  length_iqr            interquartile range (nt) of the usable footprint lengths
  frame_score           usable footprints in the main frame of their length (5' end), as a fraction
  stop_peak             stop-codon metagene (5' ends, disome lengths): highest position in STOP_PEAK over the
                        mean of STOP_BACKGROUND
  cds_utr3_density      footprints per nt in CDSs over footprints per nt in 3' UTRs (disome lengths)
  replicate_correlation mean Pearson correlation of the gene counts (log2 CPM + 1, genes with a mean of 10
                        or more) with the other libraries of the same condition
saturation.tsv: expected usable footprints at the depths of DEPTHS (M reads), from the reads of each usable
footprint: subsampling to depth d keeps a footprint with r reads with probability 1 - (1 - d / sequenced)^r.
gene_counts.tsv, correlation.tsv: the usable footprints for each gene, and the correlations.
report.pdf: the plots.

Usage: compare.py <summary.tsv> <transcripts.tsv> <disome_min> <disome_max> <gel_min> <gel_max> <outdir>
       <sample>:<regions.tsv>:<metagene.tsv>:<genes.tsv>:<dups.tsv> ...
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

STOP_PEAK = (-70, -30)
STOP_BACKGROUND = (-150, -90)
DEPTHS = [1, 2, 5, 10, 15, 20]
# categorical slots (validated palette), in fixed order; libraries of one condition differ by line style
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
STYLES = ["-", "--", ":", "-."]
INK, MUTED = "#0b0b0b", "#52514e"


def expected_usable(dups, sequenced, depth):
    p = min(depth / sequenced, 1.0)
    return float((dups["molecules"] * (1 - (1 - p) ** dups["reads"])).sum())


def main():
    if len(sys.argv) < 9:
        sys.exit(__doc__)
    summary_tsv, transcripts_tsv, dmin, dmax, gmin, gmax, outdir, *libraries = sys.argv[1:]
    dmin, dmax, gmin, gmax = int(dmin), int(dmax), int(gmin), int(gmax)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(summary_tsv, sep="\t", dtype={"sample": str}).set_index("sample")
    tx = pd.read_csv(transcripts_tsv, sep="\t")
    cds_nt = (tx["cds_end"] - tx["cds_start"]).sum()
    utr3_nt = (tx["length"] - tx["cds_end"]).sum()
    fates = list(summary.columns[summary.columns.get_loc("sequenced") + 1:summary.columns.get_loc("transcript") + 1])

    tables = {}
    for library in libraries:
        sample, *paths = library.split(":")
        tables[sample] = [pd.read_csv(p, sep="\t", dtype={"gene_id": str}) for p in paths]
    samples = list(tables)
    conditions = list(dict.fromkeys(summary.loc[samples, "condition"]))
    if len(conditions) > len(COLORS):
        sys.exit(f"{len(conditions)} conditions: the report colours at most {len(COLORS)}")
    style = {}
    for c, condition in enumerate(conditions):
        for i, s in enumerate(s for s in samples if summary.loc[s, "condition"] == condition):
            style[s] = {"color": COLORS[c], "linestyle": STYLES[i % len(STYLES)]}

    depth = np.floor(summary.loc[samples, "sequenced"].min() / 1e5) / 10  # M reads, 1 decimal, rounded down
    genes = pd.concat({s: t[2].set_index("gene_id")["molecules"] for s, t in tables.items()}, axis=1).fillna(0)
    genes = genes.astype(int).sort_index()
    kept = genes[genes.mean(axis=1) >= 10]
    logcpm = np.log2(kept / kept.sum() * 1e6 + 1)
    corr = logcpm.corr()

    rows, saturation, metagenes = [], {}, {}
    for s in samples:
        regions, metagene, _, dups = tables[s]
        row = {"sample": s, "condition": summary.loc[s, "condition"], "sequenced": summary.loc[s, "sequenced"]}
        usable = summary.loc[s, "usable"]
        row["usable"] = usable
        row["yield"] = usable / row["sequenced"]
        row[f"usable_at_{depth:g}M"] = expected_usable(dups, row["sequenced"], depth * 1e6)
        row["usable_duplicates"] = 1 - dups["molecules"].sum() / max((dups["molecules"] * dups["reads"]).sum(), 1)
        for fate in fates[:-1]:
            row[f"{fate}_fraction"] = summary.loc[s, fate] / row["sequenced"]
        row["gel_fraction"] = summary.loc[s, "gel_cds"] / max(summary.loc[s, "cds_molecules"], 1)
        row["monosome_fraction"] = summary.loc[s, "monosome_cds"] / max(summary.loc[s, "cds_molecules"], 1)

        disome = regions[(regions["length"] >= dmin) & (regions["length"] <= dmax)]
        cds = disome[disome["region"] == "cds"]
        by_length = cds.groupby("length")["molecules"].sum()
        cum = by_length.cumsum() / max(by_length.sum(), 1)
        row["length_iqr"] = (cum[cum >= 0.75].index[0] - cum[cum >= 0.25].index[0]) if len(cum) else np.nan
        frames = cds.groupby(["length", "frame"])["molecules"].sum().unstack(fill_value=0)
        row["frame_score"] = frames.max(axis=1).sum() / max(frames.values.sum(), 1)

        stop = metagene[(metagene["anchor"] == "stop") & (metagene["length"] >= dmin) & (metagene["length"] <= dmax)]
        profile = stop.groupby("position")["molecules"].sum().reindex(range(-150, 151), fill_value=0)
        background = profile.loc[STOP_BACKGROUND[0]:STOP_BACKGROUND[1]].mean()
        row["stop_peak"] = profile.loc[STOP_PEAK[0]:STOP_PEAK[1]].max() / background if background else np.nan
        utr3 = disome.loc[disome["region"] == "utr3", "molecules"].sum()
        row["cds_utr3_density"] = (cds["molecules"].sum() / cds_nt) / (utr3 / utr3_nt) if utr3 else np.nan
        same = [o for o in samples if o != s and summary.loc[o, "condition"] == row["condition"]]
        row["replicate_correlation"] = corr.loc[s, same].mean() if same and s in corr else np.nan
        rows.append(row)

        depths = [d for d in DEPTHS if d * 1e6 < row["sequenced"]] + [row["sequenced"] / 1e6]
        saturation[s] = pd.Series({d: expected_usable(dups, row["sequenced"], d * 1e6) for d in depths})
        metagenes[s] = {anchor: metagene[(metagene["anchor"] == anchor) & (metagene["length"] >= dmin)
                                         & (metagene["length"] <= dmax)].groupby("position")["molecules"].sum()
                        .reindex(range(-150, 151), fill_value=0) for anchor in ("start", "stop")}
        metagenes[s]["by_length"] = stop.pivot_table(index="length", columns="position", values="molecules",
                                                     aggfunc="sum", fill_value=0)
        metagenes[s]["frames"] = frames
        metagenes[s]["lengths"] = regions[regions["region"] == "cds"].groupby("length")["molecules"].sum()

    metrics = pd.DataFrame(rows)
    metrics.to_csv(outdir / "metrics.tsv", sep="\t", index=False, float_format="%.4g")
    sat = pd.DataFrame({s: v for s, v in saturation.items()})
    sat.index.name = "depth_M"
    sat.to_csv(outdir / "saturation.tsv", sep="\t", float_format="%.1f")
    genes.to_csv(outdir / "gene_counts.tsv", sep="\t")
    corr.to_csv(outdir / "correlation.tsv", sep="\t", float_format="%.3f")

    plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.5, "font.size": 8,
                         "lines.linewidth": 1.5})
    with PdfPages(outdir / "report.pdf") as pdf:
        # read fates: at most 8 colours, so the 4 largest fates after transcript, multi and unaligned; the
        # rest are "other"
        fractions = summary.loc[samples, fates].div(summary.loc[samples, "sequenced"], axis=0)
        rest = fractions.drop(columns=["transcript", "multi", "unaligned"]).mean().sort_values(ascending=False)
        shown = ["transcript", "multi", "unaligned"] + list(rest.index[:4])
        fractions = fractions[shown].assign(other=fractions[rest.index[4:]].sum(axis=1))
        fig, ax = plt.subplots(figsize=(8, 4.5))
        left = np.zeros(len(samples))
        for i, fate in enumerate(fractions.columns):
            ax.barh(samples, fractions[fate], left=left, color=COLORS[i], edgecolor="white", linewidth=1,
                    label=fate)
            left += fractions[fate].values
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("fraction of the sequenced reads")
        ax.set_title(f"Read fates (other: {', '.join(rest.index[4:]) or '-'})", loc="left")
        ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        def lines(ax, series, title, xlabel, ylabel, span=None):
            for s in samples:
                y = series(s)
                ax.plot(y.index, y.values, label=f"{s} ({summary.loc[s, 'condition']})", **style[s])
            if span:
                ax.axvspan(*span, color="#f0efec", zorder=0)
            ax.set_title(title, loc="left")
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        fig, axes = plt.subplots(3, 1, figsize=(8, 10))
        lines(axes[0], lambda s: metagenes[s]["lengths"] / metagenes[s]["lengths"].sum(),
              f"CDS footprint lengths (shaded: gel {gmin}-{gmax} nt)", "footprint length (nt)",
              "fraction of CDS footprints", span=(gmin, gmax))
        for ax, anchor in zip(axes[1:], ("start", "stop")):
            lines(ax, lambda s: metagenes[s][anchor] / max(metagenes[s][anchor].mean(), 1e-9),
                  f"5' ends around the {anchor} codon ({dmin}-{dmax} nt)", f"5' end to the {anchor} codon (nt)",
                  "density / mean")
        axes[0].legend(ncol=3, frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        ncol = 3
        nrow = -(-len(samples) // ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(8, 2.6 * nrow), squeeze=False)
        for ax, s in zip(axes.flat, samples):
            m = metagenes[s]["by_length"]
            if m.size:
                m = m.div(m.sum(axis=1).replace(0, 1), axis=0)
                ax.imshow(m.values, aspect="auto", cmap="Blues", origin="lower",
                          extent=(m.columns.min() - 0.5, m.columns.max() + 0.5, m.index.min() - 0.5,
                                  m.index.max() + 0.5))
            ax.set_title(f"{s} ({summary.loc[s, 'condition']})", loc="left")
            ax.set_xlabel("5' end to the stop codon (nt)")
            ax.set_ylabel("length (nt)")
            ax.grid(False)
        for ax in axes.flat[len(samples):]:
            ax.axis("off")
        fig.suptitle("5' ends around the stop codon, by length (each row scaled to its sum)", x=0.01, ha="left")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(2, 1, figsize=(8, 8))
        lines(axes[0], lambda s: metagenes[s]["frames"].max(axis=1) / metagenes[s]["frames"].sum(axis=1)
              .replace(0, 1), "Main-frame fraction of the CDS footprints, by length", "footprint length (nt)",
              "fraction in the main frame")
        lines(axes[1], lambda s: saturation[s] / 1e3, "Saturation: expected usable footprints after subsampling",
              "sequenced reads (M)", "usable footprints (thousands)")
        axes[1].legend(ncol=3, frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(corr.values, cmap="Blues", vmin=min(0.5, np.nanmin(corr.values)), vmax=1)
        ax.set_xticks(range(len(corr)), corr.columns, rotation=90)
        ax.set_yticks(range(len(corr)), corr.index)
        for i in range(len(corr)):
            for j in range(len(corr)):
                ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if corr.values[i, j] > 0.85 else INK)
        ax.grid(False)
        ax.set_title(f"Gene count correlation (Pearson, log2 CPM + 1, {len(kept)} genes)", loc="left")
        fig.colorbar(im, ax=ax, shrink=0.7)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        table = metrics.set_index("sample").drop(columns="condition").T
        fig, ax = plt.subplots(figsize=(8.5, 0.22 * len(table) + 1))
        ax.axis("off")
        cells = [[f"{v:.3g}" for v in r] for r in table.values]
        t = ax.table(cellText=cells, rowLabels=table.index,
                     colLabels=[f"{s}\n{summary.loc[s, 'condition']}" for s in table.columns], loc="center")
        t.auto_set_font_size(False)
        t.set_fontsize(6)
        ax.set_title("Metrics (metrics.tsv)", loc="left")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


if __name__ == "__main__":
    main()
