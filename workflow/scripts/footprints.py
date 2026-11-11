"""Footprints on transcripts: the computations of the figures.

This module has no workflow files and no plots: the plot_*.py scripts read the workflow files, call these
functions and draw the figures. Errors raise ValueError.

Footprint table (read_bam): one row for each molecule (a collapsed read), with the columns
  transcript  the transcript (categorical: the references of the BAM)
  start       the 5' end on the transcript (0-based)
  length      the footprint length (nt)
  reads       the reads of the molecule
Annotation table: one row for each transcript (index: transcript id), with the columns
  cds_start   the first nt of the start codon (0-based)
  cds_end     the nt after the stop codon (the stop codon is in the CDS)
"""
import re

import numpy as np
import pandas as pd

REGIONS = ["5' UTR", "CDS", "3' UTR"]
COUNT = re.compile(r"_x(\d+)_[^_]*$")  # <library>_<rank>_x<count>_<5' random nt>:<3' random nt>
NOT_PRIMARY_FORWARD = 0x4 | 0x10 | 0x100 | 0x800  # unmapped, reverse, secondary, supplementary
# standard genetic code: codons in TCAG order -> one-letter amino acids (*: stop)
CODE = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def read_bam(path):
    """The footprint table of a BAM of collapsed reads: one primary forward-strand alignment for each molecule,
    its read count in the name (x<count>, collapse_reads.py)."""
    import pysam  # only here: the other functions do not need it

    with pysam.AlignmentFile(path) as bam:
        names = list(bam.references)
        tid, start, length, reads = [], [], [], []
        for a in bam:
            if a.flag & NOT_PRIMARY_FORWARD:
                raise ValueError(f"{path}: {a.query_name} is not a primary forward-strand alignment")
            count = COUNT.search(a.query_name)
            if count is None:
                raise ValueError(f"{path}: {a.query_name} has no read count (_x<count>_)")
            tid.append(a.reference_id)
            start.append(a.reference_start)
            length.append(a.query_length)
            reads.append(int(count.group(1)))
    return pd.DataFrame({"transcript": pd.Categorical.from_codes(tid, names), "start": np.array(start, dtype=int),
                         "length": np.array(length, dtype=int), "reads": np.array(reads, dtype=int)})


def annotate(fp, annotation):
    """The footprint table with the positions of the 5' end to the start codon (to_start; 0: the A of the ATG)
    and to the stop codon (to_stop; 0: its first nt), the frame of the 5' end (frame: to_start modulo 3; 0: the
    first nt of a codon) and the region of the footprint midpoint (5' end + length // 2): 5' UTR, CDS (the stop
    codon included) or 3' UTR."""
    missing = set(fp["transcript"].unique()) - set(annotation.index)
    if missing:
        raise ValueError(f"{len(missing)} transcripts are not in the annotation, e.g. {sorted(missing)[0]}")
    cds = annotation.loc[fp["transcript"].astype(str), ["cds_start", "cds_end"]]
    cds_start, cds_end = cds["cds_start"].to_numpy(), cds["cds_end"].to_numpy()
    start = fp["start"].to_numpy()
    middle = start + fp["length"].to_numpy() // 2
    region = np.where(middle < cds_start, 0, np.where(middle < cds_end, 1, 2))
    return fp.assign(to_start=start - cds_start, to_stop=start - (cds_end - 3), frame=(start - cds_start) % 3,
                     region=pd.Categorical.from_codes(region, REGIONS))


def regions(fp, lengths):
    """The footprints of each length (rows: lengths[0] to lengths[1], both included) in each region (columns:
    REGIONS; region of annotate)."""
    counts = fp.groupby(["length", "region"], observed=False).size().unstack(fill_value=0)
    return counts.reindex(index=range(lengths[0], lengths[1] + 1), columns=REGIONS, fill_value=0)


def frames(fp, lengths):
    """The footprints of each length (rows: lengths[0] to lengths[1], both included) in each frame (columns 0, 1,
    2: frame of annotate)."""
    counts = fp.groupby(["length", "frame"]).size().unstack(fill_value=0)
    return counts.reindex(index=range(lengths[0], lengths[1] + 1), columns=range(3), fill_value=0)


def metagene(positions, window):
    """The footprints at each position of the window (first, last: both included), from their positions to a
    codon (to_start or to_stop of annotate)."""
    positions = np.asarray(positions)
    inside = positions[(positions >= window[0]) & (positions <= window[1])]
    return pd.Series(np.bincount(inside - window[0], minlength=window[1] - window[0] + 1),
                     index=range(window[0], window[1] + 1))


def translate(seq):
    """The one-letter amino acids of the codons of seq, from its first nt (X: a codon with a letter other than
    ACGT)."""
    seq = seq.upper()
    return "".join(CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def motif_sites(sequences, annotation, motifs, margins):
    """The sites of each peptide motif in the CDSs: motif label -> table with the columns transcript and anchor
    (the first nt of the first codon of the match; 0-based, on the transcript). motifs: label -> regular
    expression on the one-letter amino acids of the CDS (translate; the stop codon excluded), matched in the
    reading frame; overlapping matches are all sites. Only the sites at least margins[0] nt after the first nt
    of the start codon and margins[1] nt before the first nt of the stop codon are kept. sequences: transcript ->
    nucleotide sequence, for each transcript of the annotation."""
    missing = set(annotation.index) - set(sequences)
    if missing:
        raise ValueError(f"{len(missing)} transcripts have no sequence, e.g. {sorted(missing)[0]}")
    patterns = {label: re.compile(f"(?=(?:{pattern}))") for label, pattern in motifs.items()}
    rows = {label: [] for label in motifs}
    for tx, cds_start, cds_end in annotation[["cds_start", "cds_end"]].itertuples():
        protein = translate(sequences[tx][cds_start:cds_end - 3])
        for label, pattern in patterns.items():
            for m in pattern.finditer(protein):
                anchor = cds_start + 3 * m.start()
                if anchor - cds_start >= margins[0] and cds_end - 3 - anchor >= margins[1]:
                    rows[label].append((tx, anchor))
    return {label: pd.DataFrame(rows[label], columns=["transcript", "anchor"]) for label in motifs}


def motif_counts(fp, annotation, sites, window):
    """The 5' ends of the footprints (annotate) around the sites of a motif (motif_sites), as observed and
    expected. Observed: the 5' ends at each position of the window (first, last: both included; 0: the anchor),
    summed over the sites. Expected: the mean 5' ends per nt of the CDS (stop codon included) of the site's
    transcript, summed over the sites; the observed at each position if the 5' ends were spread evenly on each
    CDS. Returns the observed (a Series) and the expected."""
    categories = fp["transcript"].cat.categories
    site_codes = pd.Categorical(sites["transcript"], categories=categories).codes.astype(np.int64)
    if (site_codes < 0).any():
        raise ValueError(f"sites on transcripts that are not references of the footprints, e.g. "
                         f"{sites['transcript'][site_codes < 0].iloc[0]}")
    codes = fp["transcript"].cat.codes.to_numpy().astype(np.int64)
    in_cds = ((fp["to_start"] >= 0) & (fp["to_stop"] < 3)).to_numpy()
    cds = annotation.loc[sites["transcript"], ["cds_start", "cds_end"]]
    density = np.bincount(codes[in_cds], minlength=len(categories))[site_codes] / (cds["cds_end"] -
                                                                                    cds["cds_start"]).to_numpy()
    # 5' end counts by key transcript code * 2^32 + position; a last key above all others for searchsorted
    keys, counts = np.unique(codes * 2**32 + fp["start"].to_numpy(), return_counts=True)
    keys, counts = np.append(keys, np.iinfo(np.int64).max), np.append(counts, 0)
    anchors = sites["anchor"].to_numpy()
    observed = pd.Series(0, index=range(window[0], window[1] + 1))
    for position in observed.index:  # one position at a time: memory grows with the sites only
        query = site_codes * 2**32 + anchors + position
        i = np.searchsorted(keys, query)
        observed[position] = np.where((keys[i] == query) & (anchors + position >= 0), counts[i], 0).sum()
    return observed, float(density.sum())


def saturation(fp, fractions, min_molecules, rng):
    """Subsamples of the reads of the footprint table: each read is kept with the probability fraction (one
    random subsample for each fraction, from rng, a numpy Generator; fraction 1 keeps all reads). For each
    fraction: the reads, the molecules with one read or more, and the transcripts with min_molecules or more of
    these molecules."""
    codes, reads = fp["transcript"].cat.codes.to_numpy(), fp["reads"].to_numpy()
    rows = []
    for fraction in fractions:
        kept = reads if fraction == 1 else rng.binomial(reads, fraction)
        present = kept > 0
        rows.append({"fraction": fraction, "reads": int(kept.sum()), "molecules": int(present.sum()),
                     "transcripts": int((np.bincount(codes[present]) >= min_molecules).sum())})
    return pd.DataFrame(rows)
