"""Shared plot style and helpers for the plot_*.py scripts (plots.smk).

Colours: the categorical palette below is validated for colour-vision deficiency in this order, for marks that
touch (stacks, bars, lines); give series its slots in order, never cycled. Text uses the ink colours, never a
series colour.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
# label colour inside a fill: near-black on yellow, white on the other slots
LABEL_ON = {c: (INK if c == "#eda100" else "white") for c in CATEGORICAL}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10.5, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "axes.labelsize": 10.5, "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "axes.axisbelow": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "xtick.color": AXIS, "ytick.color": AXIS,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "legend.frameon": False, "legend.fontsize": 10, "lines.linewidth": 1.6, "figure.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42,
})


def read_design(path):
    """design.tsv (rule design_table): one row for each library, in plot order."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def round_ticks(axis, n=6):
    axis.set_major_locator(MaxNLocator(nbins=n, steps=[1, 2, 2.5, 5, 10]))


def percent_ticks(axis, n=6):
    round_ticks(axis, n)
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))


def row_major(items, ncol):
    """Reorder legend entries: a legend with ncol columns fills them column by column; this makes it read row by
    row."""
    rows = [items[i:i + ncol] for i in range(0, len(items), ncol)]
    return [r[c] for c in range(ncol) for r in rows if c < len(r)]


def library_rows(design):
    """y positions of one bar per library (design order, from the top), with a gap and a label before each
    replicate if there are several. Returns the positions, the (y, text) labels and the total height."""
    replicates = design["replicate"]
    grouped = replicates.ne("").all() and replicates.nunique() > 1
    y, labels, pos = [], [], 0.0
    for i, rep in enumerate(replicates):
        if grouped and (i == 0 or rep != replicates.iloc[i - 1]):
            pos += 0.6 if i == 0 else 1.1
            labels.append((pos - 0.75, f"Replicate {rep}"))
        y.append(pos)
        pos += 1
    return y, labels, pos


def library_axis(ax, design, y, labels, height):
    """Look of the one-bar-per-library charts: library and replicate labels, no horizontal grid."""
    for yl, text in labels:
        ax.text(0, yl, text, va="center", ha="left", fontsize=10, fontweight="bold", color=INK)
    ax.set_yticks(y, [f"{r.sample_id}  {r.condition}" for r in design.itertuples()])
    ax.tick_params(axis="y", length=0, labelcolor=INK)
    ax.set_ylim(height, -0.2)
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)


def save(fig, prefix, data, caption):
    """<prefix>.pdf, .png, .tsv (the numbers in the figure) and .txt (the caption)."""
    prefix = Path(prefix)
    fig.savefig(prefix.with_suffix(".pdf"))
    fig.savefig(prefix.with_suffix(".png"), dpi=200)
    plt.close(fig)
    data.to_csv(prefix.with_suffix(".tsv"), sep="\t", float_format="%.6g")
    prefix.with_suffix(".txt").write_text(caption.strip() + "\n")
