"""Footprint positions of one library on the transcripts (the BAM of bowtie_align).

Each BAM record is one molecule with <count> reads (x<count> in its name, collapse_reads.py). Its
approximate A-site, ASITE_FROM_3P nt before its 3' end (of the leading ribosome of a disome; approximately
also of a monosome), gives its region: utr5, cds or utr3. The step that calibrates the offsets for each
length replaces this approximation. Writes:

  regions.tsv   length, region, frame of the 5' end to the CDS start (cds only, else -), molecules, reads
  metagene.tsv  anchor (start: first nt of the start codon; stop: first nt of the stop codon), length,
                position of the 5' end to the anchor (-WINDOW..WINDOW), molecules
  genes.tsv     gene_id, molecules: the usable footprints (disome length, A-site in the CDS)
  dups.tsv      reads, molecules: the usable footprints by their number of reads (saturation curve)

Usage: qc_library.py <in.bam> <transcripts.tsv> <disome_min> <disome_max> <regions.tsv> <metagene.tsv>
       <genes.tsv> <dups.tsv>
"""
import re
import sys
from collections import Counter

import pysam

ASITE_FROM_3P = 15
WINDOW = 150
COUNT = re.compile(r"_x(\d+)_[^_]*$")


def write(path, header, counter):
    with open(path, "w") as o:
        o.write("\t".join(header) + "\n")
        for key in sorted(counter):
            o.write("\t".join(map(str, (*key, counter[key]))) + "\n")


def main():
    if len(sys.argv) != 9:
        sys.exit(__doc__)
    bam, transcripts, dmin, dmax, out_regions, out_metagene, out_genes, out_dups = sys.argv[1:]
    dmin, dmax = int(dmin), int(dmax)

    cds = {}
    with open(transcripts) as f:
        next(f)
        for line in f:
            tx, gene, _, start, end = line.split()
            cds[tx] = (gene, int(start), int(end))

    molecules, reads, metagene, genes, dups = Counter(), Counter(), Counter(), Counter(), Counter()
    with pysam.AlignmentFile(bam) as f:
        for a in f:
            if a.is_secondary or a.is_unmapped:
                continue
            gene, start, end = cds[a.reference_name]
            five, length = a.reference_start, a.query_length
            asite = a.reference_end - ASITE_FROM_3P
            region = "utr5" if asite < start else "cds" if asite < end else "utr3"
            frame = (five - start) % 3 if region == "cds" else "-"
            count = int(COUNT.search(a.query_name).group(1))
            molecules[length, region, frame] += 1
            reads[length, region, frame] += count
            for anchor, pos in (("start", five - start), ("stop", five - (end - 3))):
                if -WINDOW <= pos <= WINDOW:
                    metagene[anchor, length, pos] += 1
            if region == "cds" and dmin <= length <= dmax:
                genes[gene,] += 1
                dups[count,] += 1

    with open(out_regions, "w") as o:
        o.write("length\tregion\tframe\tmolecules\treads\n")
        for key in sorted(molecules, key=lambda k: (k[0], k[1], str(k[2]))):
            o.write("\t".join(map(str, (*key, molecules[key], reads[key]))) + "\n")
    write(out_metagene, ["anchor", "length", "position", "molecules"], metagene)
    write(out_genes, ["gene_id", "molecules"], genes)
    write(out_dups, ["reads", "molecules"], dups)


if __name__ == "__main__":
    main()
