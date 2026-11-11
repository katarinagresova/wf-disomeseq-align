"""Tests of workflow/scripts/footprints.py, with synthetic footprints.

Run from the repository: python -m pytest .test/unit (the packages of workflow/envs/python.yaml and pytest).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow" / "scripts"))
import footprints as F  # noqa: E402

TRANSCRIPTS = ["tx1", "tx2"]
# tx1: 5' UTR 0-99, CDS 100-399 (stop codon 397-399), 3' UTR 400-599
ANNOTATION = pd.DataFrame({"cds_start": [100, 50], "cds_end": [400, 200]}, index=TRANSCRIPTS)


def table(rows):
    """A footprint table from (transcript, start, length, reads) rows."""
    fp = pd.DataFrame(rows, columns=["transcript", "start", "length", "reads"])
    fp["transcript"] = pd.Categorical(fp["transcript"], categories=TRANSCRIPTS)
    return fp


def write_bam(path, records):
    """A BAM on tx1 (600 nt) and tx2 (300 nt) from (name, transcript, start, length, flag) records."""
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "tx1", "LN": 600}, {"SN": "tx2", "LN": 300}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for name, transcript, start, length, flag in records:
            a = pysam.AlignedSegment(bam.header)
            a.query_name, a.query_sequence, a.flag = name, "A" * length, flag
            a.reference_id, a.reference_start = TRANSCRIPTS.index(transcript), start
            a.cigarstring, a.mapping_quality = f"{length}M", 255
            bam.write(a)
    return path


def test_read_bam(tmp_path):
    bam = write_bam(tmp_path / "a.bam", [("lib_x1_0_x12_ACGTACGT:ACGT", "tx2", 7, 60, 0),
                                         ("lib_1_x3_AAAAAAAA:CCCC", "tx1", 352, 61, 0)])
    fp = F.read_bam(bam)
    assert list(fp.columns) == ["transcript", "start", "length", "reads"]
    assert list(fp["transcript"].cat.categories) == TRANSCRIPTS
    assert fp.astype({"transcript": str}).values.tolist() == [["tx2", 7, 60, 12], ["tx1", 352, 61, 3]]


@pytest.mark.parametrize("name, flag", [("lib_0_x2_A:C", 0x4), ("lib_0_x2_A:C", 0x10), ("lib_0_x2_A:C", 0x100),
                                        ("lib_0_x2_A:C", 0x800), ("lib_0_2_A:C", 0)])
def test_read_bam_rejects(tmp_path, name, flag):
    bam = write_bam(tmp_path / "a.bam", [("lib_0_x1_A:C", "tx1", 0, 30, 0), (name, "tx1", 0, 30, flag)])
    with pytest.raises(ValueError):
        F.read_bam(bam)


def test_annotate_positions():
    fp = F.annotate(table([("tx1", 100, 60, 1), ("tx1", 352, 60, 1), ("tx2", 0, 30, 1)]), ANNOTATION)
    assert fp["to_start"].tolist() == [0, 252, -50]
    assert fp["to_stop"].tolist() == [-297, -45, -197]
    assert fp["frame"].tolist() == [0, 0, 1]


def test_annotate_regions():
    # midpoint (5' end + 30) at 99, 100, 399, 400: the last nt of the 5' UTR, the first nt of the CDS, the last
    # nt of the stop codon, the first nt of the 3' UTR
    fp = F.annotate(table([("tx1", s, 60, 1) for s in (69, 70, 369, 370)]), ANNOTATION)
    assert fp["region"].tolist() == ["5' UTR", "CDS", "CDS", "3' UTR"]
    assert list(fp["region"].cat.categories) == F.REGIONS


def test_annotate_rejects_missing_transcript():
    with pytest.raises(ValueError):
        F.annotate(table([("tx2", 0, 30, 1)]), ANNOTATION.loc[["tx1"]])


def test_regions():
    # midpoints (5' end + 27) at 96 (5' UTR), 127 and 397 (CDS), 427 (3' UTR) of tx1
    fp = F.annotate(table([("tx1", 69, 54, 1), ("tx1", 100, 54, 4), ("tx1", 370, 54, 1), ("tx1", 400, 55, 1)]),
                    ANNOTATION)
    counts = F.regions(fp, (53, 55))
    assert list(counts.index) == [53, 54, 55] and list(counts.columns) == F.REGIONS
    assert counts.values.tolist() == [[0, 0, 0], [1, 2, 0], [0, 0, 1]]  # molecules, not reads
    assert F.regions(fp.iloc[:0], (54, 54)).values.tolist() == [[0, 0, 0]]


def test_frames():
    # 5' ends at codon nt 1, 2 and 3 (frames 0, 1, 2) of tx1 (CDS from 100)
    fp = F.annotate(table([("tx1", 100, 54, 1), ("tx1", 103, 54, 5), ("tx1", 104, 54, 1), ("tx1", 105, 56, 1)]),
                    ANNOTATION)
    counts = F.frames(fp, (53, 56))
    assert list(counts.index) == [53, 54, 55, 56] and list(counts.columns) == [0, 1, 2]
    assert counts.values.tolist() == [[0, 0, 0], [2, 1, 0], [0, 0, 0], [0, 0, 1]]  # molecules, not reads
    assert F.frames(fp.iloc[:0], (54, 55)).values.tolist() == [[0, 0, 0], [0, 0, 0]]


def test_metagene_window():
    profile = F.metagene(pd.Series([-61, -60, -60, 0, 150, 151]), (-60, 150))
    assert list(profile.index) == list(range(-60, 151))
    assert (profile[-60], profile[0], profile[150], profile.sum()) == (2, 1, 1, 4)


def test_saturation():
    fp = table([("tx1", 0, 60, 1), ("tx1", 3, 60, 5), ("tx2", 0, 60, 2)])
    curve = F.saturation(fp, [0, 0.5, 1], 2, np.random.default_rng(1))
    assert curve.loc[[0, 2], ["reads", "molecules", "transcripts"]].values.tolist() == [[0, 0, 0], [8, 3, 1]]
    assert 0 <= curve.loc[1, "molecules"] <= 3 and 0 <= curve.loc[1, "reads"] <= 8


def test_translate():
    assert F.translate("ATGccaNNNTAAG") == "MPX*"


def cds_sequence(length, codons, utr):
    """A transcript sequence: utr nt, then the CDS (ATG, GCT except the codons {index: codon}, TAA), then utr nt
    to the length."""
    cds = ["ATG"] + [codons.get(k, "GCT") for k in range(1, 99)] + ["TAA"]
    return (utr * 100)[:100] + "".join(cds) + (utr * 100)[:length - 400]


def test_motif_sites():
    # tx1 (CDS 100-400): Pro-Pro-Pro at codons 10-12 (nt 130), Arg-Ala-Lys at codons 20-22 (nt 160), Pro-Pro at
    # codons 90-91 (nt 370: 27 nt before the stop codon); 5' and 3' UTRs that code Pro-Pro in frame
    sequences = {"tx1": cds_sequence(600, {10: "CCA", 11: "CCG", 12: "CCT", 20: "CGT", 22: "AAA", 90: "CCC",
                                           91: "CCC"}, "CCA"),
                 "tx2": "A" * 50 + "ATG" + "GCT" * 48 + "TAG" + "A" * 100}
    sites = F.motif_sites(sequences, ANNOTATION, {"PP": "PP", "RXK": "R.K"}, (30, 30))
    assert list(sites) == ["PP", "RXK"]
    assert sites["PP"].values.tolist() == [["tx1", 130], ["tx1", 133]]  # overlapping matches; not 370
    assert sites["RXK"].values.tolist() == [["tx1", 160]]
    assert F.motif_sites(sequences, ANNOTATION, {"PP": "PP"}, (0, 0))["PP"]["anchor"].tolist() == [130, 133, 370]
    with pytest.raises(ValueError):
        F.motif_sites({"tx1": sequences["tx1"]}, ANNOTATION, {"PP": "PP"}, (30, 30))


def test_motif_counts():
    # CDS 5' ends: tx1 (300 nt) 195 (2 molecules), 250, 260; tx2 (150 nt) 120. Not in a CDS: tx1 50, tx2 200
    fp = F.annotate(table([("tx1", 195, 60, 1), ("tx1", 195, 60, 5), ("tx1", 250, 60, 1), ("tx1", 260, 60, 1),
                           ("tx1", 50, 60, 1), ("tx2", 120, 60, 1), ("tx2", 200, 60, 1)]), ANNOTATION)
    sites = pd.DataFrame({"transcript": ["tx1", "tx1", "tx2"], "anchor": [200, 250, 120]})
    observed, expected = F.motif_counts(fp, ANNOTATION, sites, (-5, 5))
    assert list(observed.index) == list(range(-5, 6))
    assert (observed[-5], observed[0], observed.sum()) == (2, 2, 4)  # molecules; tx2 200 is not at tx1 200
    assert expected == pytest.approx(4 / 300 + 4 / 300 + 1 / 150)
    observed, expected = F.motif_counts(fp, ANNOTATION, sites.iloc[:0], (-5, 5))
    assert (observed.sum(), expected) == (0, 0)
    with pytest.raises(ValueError):
        F.motif_counts(fp, ANNOTATION, pd.DataFrame({"transcript": ["tx3"], "anchor": [0]}), (-5, 5))
