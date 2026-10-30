"""Transcripts to align to and ncRNA contaminants, from a genome fasta and its NCBI RefSeq GTF.

Protein-coding (transcript_biotype mRNA): one transcript for each gene. A candidate has a complete
CDS: not partial, CDS + stop_codon contiguous on the transcript, a multiple of 3 long, ATG start,
TAA / TAG / TGA stop and no stop before it. The last excludes the models that NCBI corrected for
genome errors ("unclassified transcription discrepancy": the genome sequence has a frameshift) and
selenoproteins. Of the candidates: the longest CDS, then the longest transcript, then the
smallest transcript_id. Genes without a candidate are left out (counted in the log).

ncRNA contaminants: the transcripts with transcript_biotype rRNA, tRNA, snRNA, snoRNA, guide_RNA
(scaRNA) or miRNA, named <type>_<transcript_id>, type rRNA, tRNA or ncRNA (build_contaminants.py).
lnc_RNA and misc_RNA (non-coding isoforms of coding genes) are not contaminants.

Writes transcripts.fa (>transcript_id), transcripts.tsv (transcript_id, gene_id, length, cds_start,
cds_end: 0-based on the transcript, end exclusive, stop codon included) and ncrna.fa. pysam indexes
the genome (.fai next to it) if it has no index.

Usage: extract_transcripts.py <genome.fa> <annotation.gtf> <transcripts.fa> <transcripts.tsv> <ncrna.fa>
"""
import re
import sys
from collections import Counter, defaultdict

import pysam

NCRNA_TYPES = {"rRNA": "rRNA", "tRNA": "tRNA", "snRNA": "ncRNA", "snoRNA": "ncRNA", "guide_RNA": "ncRNA",
               "miRNA": "ncRNA"}
STOPS = {"TAA", "TAG", "TGA"}
ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_gtf(path):
    """transcript_id -> its transcript attributes, plus its exons and CDS + stop_codon intervals."""
    info, exons, cds = {}, defaultdict(list), defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = line.rstrip("\n").split("\t")
            if feature not in ("transcript", "exon", "CDS", "stop_codon"):
                continue
            attrs = dict(ATTRIBUTE.findall(attributes))
            tx = attrs["transcript_id"]
            interval = (int(start) - 1, int(end))  # 0-based, end exclusive
            if feature == "transcript":
                info[tx] = {"chrom": chrom, "strand": strand, "gene_id": attrs["gene_id"],
                            "biotype": attrs.get("transcript_biotype"), "partial": attrs.get("partial") == "true"}
            elif feature == "exon":
                exons[tx].append((chrom, strand, interval))
            else:
                cds[tx].append(interval)
    return info, exons, cds


def transcript_sequence(genome, chrom, strand, intervals):
    seq = "".join(genome.fetch(chrom, s, e) for s, e in sorted(intervals)).upper()
    return seq.translate(COMPLEMENT)[::-1] if strand == "-" else seq


def cds_span(strand, intervals, cds_intervals):
    """(cds_start, cds_end) on the transcript, or None if the CDS is not contiguous on it."""
    lo, hi = min(s for s, _ in cds_intervals), max(e for _, e in cds_intervals)
    first, last = (lo, hi - 1) if strand == "+" else (hi - 1, lo)  # 5'-most and 3'-most CDS nt
    offset, start, end = 0, None, None
    for s, e in sorted(intervals, reverse=(strand == "-")):
        if s <= first < e:
            start = offset + (first - s if strand == "+" else e - 1 - first)
        if s <= last < e:
            end = offset + (last - s if strand == "+" else e - 1 - last) + 1
        offset += e - s
    if start is None or end is None or end - start != sum(e - s for s, e in cds_intervals):
        return None
    return start, end


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    genome_fa, gtf, out_fa, out_tsv, out_ncrna = sys.argv[1:]

    info, exons, cds = read_gtf(gtf)
    genome = pysam.FastaFile(genome_fa)
    log = Counter()

    candidates = defaultdict(list)  # gene -> (sort key, transcript, seq, cds start, cds end)
    with open(out_ncrna, "w") as ncrna:
        for tx, t in info.items():
            places = {(chrom, strand) for chrom, strand, _ in exons[tx]}
            if places != {(t["chrom"], t["strand"])}:
                log["skipped: exons missing or on another chromosome / strand"] += 1
                continue
            intervals = [interval for _, _, interval in exons[tx]]
            if t["biotype"] in NCRNA_TYPES:
                seq = transcript_sequence(genome, t["chrom"], t["strand"], intervals)
                ncrna.write(f">{NCRNA_TYPES[t['biotype']]}_{tx}\n{seq}\n")
                log[f"ncRNA contaminant: {t['biotype']}"] += 1
                continue
            if t["biotype"] != "mRNA":
                continue
            log["mRNA"] += 1
            if t["partial"] or not cds[tx]:
                log["mRNA not a candidate: partial or no CDS"] += 1
                continue
            span = cds_span(t["strand"], intervals, cds[tx])
            if span is None:
                log["mRNA not a candidate: CDS not contiguous on the transcript"] += 1
                continue
            seq = transcript_sequence(genome, t["chrom"], t["strand"], intervals)
            start, end = span
            if (end - start) % 3 or seq[start:start + 3] != "ATG" or seq[end - 3:end] not in STOPS:
                log["mRNA not a candidate: CDS length, start or stop codon"] += 1
                continue
            if any(seq[i:i + 3] in STOPS for i in range(start, end - 3, 3)):
                log["mRNA not a candidate: stop codon inside the CDS"] += 1
                continue
            candidates[t["gene_id"]].append(((-(end - start), -len(seq), tx), tx, seq, start, end))

    genes = {t["gene_id"] for t in info.values() if t["biotype"] == "mRNA"}
    log["protein-coding genes"] = len(genes)
    log["protein-coding genes without a candidate (left out)"] = len(genes - candidates.keys())
    with open(out_fa, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("transcript_id\tgene_id\tlength\tcds_start\tcds_end\n")
        for gene in sorted(candidates):
            _, tx, seq, start, end = min(candidates[gene])
            fa.write(f">{tx}\n{seq}\n")
            tsv.write(f"{tx}\t{gene}\t{len(seq)}\t{start}\t{end}\n")
            log["transcripts written"] += 1
    for key, n in sorted(log.items()):
        print(f"{key}\t{n}", file=sys.stderr)


if __name__ == "__main__":
    main()
