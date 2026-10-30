"""The contaminant set for bowtie_index: the ncRNAs of the annotation and the configured fastas.

ncrna.fa (extract_transcripts.py) has records named <type>_<id>, sense strand only. Each record of a
<type>=<fasta> argument is added on both strands (<id>, <id>_rc): mtDNA is transcribed from both.
Each sequence is written once (the first record with it), as >{n}_{type}_{id}: assign_reads.py takes
the type from the name, and bowtie_index needs the names unique.

Usage: build_contaminants.py <out.fa> <ncrna.fa> [<type>=<fasta> ...]
"""
import sys

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    out_fa, ncrna_fa, *extra = sys.argv[1:]

    records = [tuple(name.split("_", 1)) + (seq,) for name, seq in read_fasta(ncrna_fa)]
    for arg in extra:
        rtype, path = arg.split("=", 1)
        for name, seq in read_fasta(path):
            records += [(rtype, name, seq), (rtype, f"{name}_rc", seq.translate(COMPLEMENT)[::-1])]

    seen = set()
    with open(out_fa, "w") as out:
        for rtype, name, seq in records:
            if seq in seen:
                continue
            seen.add(seq)
            out.write(f">{len(seen)}_{rtype}_{name}\n{seq}\n")
    print(f"records\t{len(records)}\nwritten\t{len(seen)}", file=sys.stderr)


if __name__ == "__main__":
    main()
