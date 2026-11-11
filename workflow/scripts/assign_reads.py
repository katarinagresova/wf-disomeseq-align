"""Assign each read from its bowtie alignments (bowtie_align): a contaminant type, multi, transcript, spike_in or unaligned.

bowtie aligns each trimmed read once to the transcripts (one for each gene) of the main species and of the
optional spike-in, and the contaminants, in one index, forward strand only, end to end with at most 2
mismatches, up to 256 alignments, best first. This reads its SAM on stdin. The best stratum (the alignments
with the fewest mismatches) gives the class:

  <type>      a contaminant in the best stratum, also if a transcript is there too: the type of the
              one with the smallest n (from its name, <n>_<type>_<id>, build_contaminants.py), so a tie
              between types goes to the main annotation, then the configured fastas, then the spike-in
  multi       more than one alignment in the best stratum, all to transcripts (of one or both species)
  transcript  one alignment in the best stratum, to a main transcript (a worse alignment elsewhere is allowed)
  spike_in    one alignment in the best stratum, to a spike-in transcript
  unaligned   no alignment

The transcript and spike_in reads go to one BAM each (uncompressed, in input order; header with the
transcripts of that species only, no base qualities). Without a spike-in (<spike_in.fa> is -), its BAM has
no records. counts.tsv has one row for each class and read length: molecules (records) and reads (the
x<count> of the collapsed read name, collapse_reads.py). Fails unless each read's alignments are best first.

Usage: assign_reads.py <contaminants.fa> <spike_in.fa | -> <out.bam> <spike_in.bam> <counts.tsv> < bowtie.sam
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
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    contaminants_fa, spike_in_fa, bam_out, spike_in_bam_out, counts_out = sys.argv[1:]

    contaminant_names = fasta_names(contaminants_fa)
    spike_in_names = fasta_names(spike_in_fa) if spike_in_fa != "-" else set()
    molecules, copies = Counter(), Counter()
    with pysam.AlignmentFile("-", "r") as sam:
        names = sam.references
        if len(set(names)) != len(names):
            sys.exit("record names in the index are not unique")
        if not (contaminant_names | spike_in_names) <= set(names):
            sys.exit(f"records of {contaminants_fa} or {spike_in_fa} are not in the index")
        # None for a transcript
        ctype = [n.split("_", 2)[1] if n in contaminant_names else None for n in names]
        cnumber = [int(n.split("_", 1)[0]) if n in contaminant_names else None for n in names]
        species = ["spike_in" if n in spike_in_names else "transcript" for n in names]
        out, out_tid = {}, {}
        for cls, path in (("transcript", bam_out), ("spike_in", spike_in_bam_out)):
            tids = [t for t, c in enumerate(ctype) if c is None and species[t] == cls]
            out_tid.update({tid: i for i, tid in enumerate(tids)})
            header = sam.header.to_dict()
            header["SQ"] = [header["SQ"][t] for t in tids]
            out[cls] = pysam.AlignmentFile(path, "wbu", header=header)

        for name, group in reads(sam):
            if group[0].is_unmapped:
                cls = "unaligned"
            else:
                nms = [a.get_tag("NM") for a in group]
                if any(x > y for x, y in zip(nms, nms[1:])):
                    sys.exit(f"{name}: alignments not best first (NM {nms})")
                best = group[:nms.count(nms[0])]
                hits = [a.reference_id for a in best if ctype[a.reference_id] is not None]
                cls = (ctype[min(hits, key=cnumber.__getitem__)] if hits
                       else "multi" if len(best) > 1 else species[best[0].reference_id])
            key = (cls, group[0].query_length)
            molecules[key] += 1
            copies[key] += int(COUNT.search(name).group(1))
            if cls in out:
                a = group[0]
                a.reference_id = out_tid[a.reference_id]
                a.mapping_quality = 255
                a.set_tags([("NH", 1), ("NM", a.get_tag("NM")), ("MD", a.get_tag("MD"))])
                a.query_qualities = None
                out[cls].write(a)
        for f in out.values():
            f.close()

    with open(counts_out, "w") as o:
        o.write("class\tlength\tmolecules\treads\n")
        for cls, length in sorted(molecules):
            o.write(f"{cls}\t{length}\t{molecules[cls, length]}\t{copies[cls, length]}\n")
    print(f"molecules\t{molecules.total()}\nreads\t{copies.total()}", file=sys.stderr)


if __name__ == "__main__":
    main()
