"""Remove retrocopies (processed pseudogenes annotated as protein-coding genes) from the transcripts to align to.

A retrocopy is a copy of a spliced mRNA in the genome: it has no introns, and its sequence is almost the same as
that of its parent gene. The annotation can give it a protein-coding model (e.g. "glyceraldehyde-3-phosphate
dehydrogenase-like"). A footprint of the parent then aligns equally well to both, is multi-mapping and is lost.

A gene is a retrocopy if its transcript has at most MAX_EXONS exons (a model of a retrocopy can have a small
second exon) and at least MIN_FRACTION of its CDS is also in a gene with at least MIN_PARENT_EXONS exons: the
CDS is cut into tiles of TILE nt, every STEP nt, and a tile is in that gene if it aligns to its transcript with
up to 2 mismatches (bowtie, rule retrocopies). Only genes with few exons can be removed, and only for a parent
with many exons. Intronless gene families (e.g. the replication-dependent histones) keep their genes, and
paralogs with introns (e.g. actins) are never removed.

  tiles  <transcripts.fa> <transcripts.tsv> <partners.fa> <tiles.fa>
         writes the transcripts with >= MIN_PARENT_EXONS exons (partners.fa) and the CDS tiles of the
         transcripts with <= MAX_EXONS exons (tiles.fa, >transcript_id:position)
  filter <transcripts.fa> <transcripts.tsv> <hits.txt> <out.fa> <out.tsv> <retrocopies.tsv>
         hits.txt: bowtie output of the tiles, columns tile and partner transcript. Writes the transcripts
         without the retrocopies, and retrocopies.tsv (transcript_id, gene_id, exons, tiles, tiles_in_parent,
         fraction, parent_gene: the gene with the most tiles)

The counts go to stderr.
"""
import sys
from collections import Counter, defaultdict

TILE, STEP = 30, 10
MAX_EXONS, MIN_PARENT_EXONS, MIN_FRACTION = 2, 3, 0.9


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
    return {name: "".join(chunks) for name, chunks in seqs.items()}


def read_tsv(path):
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        return header, [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]


def n_tiles(row):
    return max(0, (int(row["cds_end"]) - int(row["cds_start"]) - TILE) // STEP + 1)


def tiles(fa, tsv, partners_fa, tiles_fa):
    seqs = read_fasta(fa)
    _, rows = read_tsv(tsv)
    n = Counter()
    with open(partners_fa, "w") as partners, open(tiles_fa, "w") as out:
        for r in rows:
            tx, exons = r["transcript_id"], int(r["exons"])
            if exons >= MIN_PARENT_EXONS:
                partners.write(f">{tx}\n{seqs[tx]}\n")
                n["partner transcripts"] += 1
            if exons <= MAX_EXONS:
                cds = seqs[tx][int(r["cds_start"]):int(r["cds_end"])]
                for i in range(n_tiles(r)):
                    out.write(f">{tx}:{i * STEP}\n{cds[i * STEP:i * STEP + TILE]}\n")
                n["tiles"] += n_tiles(r)
                n["candidate transcripts"] += 1
    for key, value in n.items():
        print(f"{key}\t{value}", file=sys.stderr)


def filter_retrocopies(fa, tsv, hits_txt, out_fa, out_tsv, out_retrocopies):
    seqs = read_fasta(fa)
    header, rows = read_tsv(tsv)
    by_tx = {r["transcript_id"]: r for r in rows}
    in_parent = defaultdict(set)  # transcript -> its tiles that align to a partner
    parents = defaultdict(Counter)  # transcript -> partner gene -> tiles
    with open(hits_txt) as f:
        for line in f:
            tile, partner = line.rstrip("\n").split("\t")
            tx = tile.rsplit(":", 1)[0]
            if int(by_tx[tx]["exons"]) > MAX_EXONS or int(by_tx[partner]["exons"]) < MIN_PARENT_EXONS:
                sys.exit(f"{tile} -> {partner}: not a candidate tile and a partner transcript")
            in_parent[tx].add(tile)
            parents[tx][by_tx[partner]["gene_id"]] += 1

    retrocopies = []
    for tx, hit in in_parent.items():
        r = by_tx[tx]
        if len(hit) >= MIN_FRACTION * n_tiles(r):
            retrocopies.append((r, len(hit), parents[tx].most_common(1)[0][0]))
    removed = {r["transcript_id"] for r, _, _ in retrocopies}

    with open(out_fa, "w") as fa_out, open(out_tsv, "w") as tsv_out:
        tsv_out.write("\t".join(header) + "\n")
        for r in rows:
            if r["transcript_id"] not in removed:
                fa_out.write(f">{r['transcript_id']}\n{seqs[r['transcript_id']]}\n")
                tsv_out.write("\t".join(r[c] for c in header) + "\n")
    with open(out_retrocopies, "w") as out:
        out.write("transcript_id\tgene_id\texons\ttiles\ttiles_in_parent\tfraction\tparent_gene\n")
        for r, hit, parent in sorted(retrocopies, key=lambda x: x[0]["gene_id"]):
            out.write(f"{r['transcript_id']}\t{r['gene_id']}\t{r['exons']}\t{n_tiles(r)}\t{hit}\t"
                      f"{hit / n_tiles(r):.3f}\t{parent}\n")
    print(f"transcripts\t{len(rows)}\nretrocopies removed\t{len(removed)}\ntranscripts written\t"
          f"{len(rows) - len(removed)}", file=sys.stderr)


def main():
    if len(sys.argv) == 6 and sys.argv[1] == "tiles":
        tiles(*sys.argv[2:])
    elif len(sys.argv) == 8 and sys.argv[1] == "filter":
        filter_retrocopies(*sys.argv[2:])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
