"""Assign each read from its bowtie alignments (bowtie_align): a contaminant type, multi, transcript or unaligned.

bowtie aligns each trimmed read once to the transcripts (one for each gene) and the contaminants in one
index, forward strand only, end to end with at most 2 mismatches, up to 256 alignments, best first. This
reads its SAM on stdin. The best stratum (the alignments with the fewest mismatches) gives the class:

  <type>      a contaminant in the best stratum, also if a transcript is there too: the type of the
              first one (from its name, <n>_<type>_<id>, build_contaminants.py)
  multi       more than one alignment in the best stratum, all to transcripts
  transcript  one alignment in the best stratum, to a transcript (a worse contaminant alignment is allowed)
  unaligned   no alignment

The transcript reads go to the BAM (uncompressed, in input order; header without the contaminants, no
base qualities). counts.tsv has one row for each class and read length: molecules (records) and reads
(the x<count> of the collapsed read name, collapse_reads.py). Fails unless each read's alignments are
best first.

Usage: assign_reads.py <contaminants.fa> <out.bam> <counts.tsv> < bowtie.sam
"""
import re
import sys
from collections import Counter

import pysam

COUNT = re.compile(r"_x(\d+)_[^_]*$")  # <library>_<rank>_x<count>_<5' random nt>:<3' random nt>


def fasta_names(path):
    with open(path) as f:
        return {line[1:].split()[0] for line in f if line.startswith(">")}


def reads(sam):
    """(name, its records), per read: bowtie writes a read's records together."""
    name, group = None, []
    for a in sam:
        if a.query_name != name:
            if group:
                yield name, group
            name, group = a.query_name, []
        group.append(a)
    if group:
        yield name, group


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    contaminants_fa, bam_out, counts_out = sys.argv[1:]

    contaminant_names = fasta_names(contaminants_fa)
    molecules, copies = Counter(), Counter()
    with pysam.AlignmentFile("-", "r") as sam:
        names = sam.references
        if not contaminant_names <= set(names):
            sys.exit(f"{len(contaminant_names - set(names))} records of {contaminants_fa} are not in the index")
        # None for a transcript
        ctype = [n.split("_", 2)[1] if n in contaminant_names else None for n in names]
        header = sam.header.to_dict()
        header["SQ"] = [sq for sq in header["SQ"] if sq["SN"] not in contaminant_names]
        out_tid = {tid: i for i, tid in enumerate(t for t, c in enumerate(ctype) if c is None)}

        with pysam.AlignmentFile(bam_out, "wbu", header=header) as out:
            for name, group in reads(sam):
                if group[0].is_unmapped:
                    cls = "unaligned"
                else:
                    nms = [a.get_tag("NM") for a in group]
                    if any(x > y for x, y in zip(nms, nms[1:])):
                        sys.exit(f"{name}: alignments not best first (NM {nms})")
                    best = group[:nms.count(nms[0])]
                    types = [ctype[a.reference_id] for a in best if ctype[a.reference_id] is not None]
                    cls = types[0] if types else "multi" if len(best) > 1 else "transcript"
                key = (cls, group[0].query_length)
                molecules[key] += 1
                copies[key] += int(COUNT.search(name).group(1))
                if cls == "transcript":
                    a = group[0]
                    a.reference_id = out_tid[a.reference_id]
                    a.mapping_quality = 255
                    a.set_tags([("NH", 1), ("NM", a.get_tag("NM")), ("MD", a.get_tag("MD"))])
                    a.query_qualities = None
                    out.write(a)

    with open(counts_out, "w") as o:
        o.write("class\tlength\tmolecules\treads\n")
        for cls, length in sorted(molecules):
            o.write(f"{cls}\t{length}\t{molecules[cls, length]}\t{copies[cls, length]}\n")
    print(f"molecules\t{molecules.total()}\nreads\t{copies.total()}", file=sys.stderr)


if __name__ == "__main__":
    main()
