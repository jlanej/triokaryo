#!/usr/bin/env python3
"""How a diploid assembly represents the sex chromosomes (or any chromosome): per haplotype, the reference bases its contigs
align to, twice over where two copies were assembled separately, and the bases covered once (the union), so that an extra
copy collapsed into one haplotype (1x aligned, 1x covered), assembled twice in one haplotype (2x aligned, 1x covered) or
split over the two haplotypes (1x in each) can be told apart.

  assembly_sexchrom.py --out table.tsv hap1=hap1.tsv hap2=hap2.paf [...] [--chroms chrX,chrY,chr21] [--min-mapq 20]

An input is a PAF from minimap2 (query, qlen, qstart, qend, strand, target, tlen, tstart, tend, matches, alnlen, mapq) or a
table with hap, chrom, contig, start, end, span, mapq, kind (kind primary/secondary/supplementary). Primary and supplementary
alignments count, the latter being the further pieces of a long contig; secondary ones, the same piece placed again, do not.
The pseudoautosomal regions of GRCh38 are excluded from the X and the Y."""
import argparse
import sys

PAR = {"chrX": [(10_001, 2_781_479), (155_701_383, 156_030_895)], "chrY": [(10_001, 2_781_479), (56_887_903, 57_217_415)]}
LENGTH = {"chrX": 156_040_895, "chrY": 57_227_415, "chr21": 46_709_983}


def clip(chrom, s, e):
    """The interval minus the PARs: a list of (start, end)."""
    pieces = [(s, e)]
    for a, b in PAR.get(chrom, []):
        out = []
        for x, y in pieces:
            if y <= a or x >= b:
                out.append((x, y))
            else:
                if x < a:
                    out.append((x, a))
                if y > b:
                    out.append((b, y))
        pieces = out
    return [(x, y) for x, y in pieces if y > x]


def read(label, path, min_mapq):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if not f or f[0].startswith("#"):
            continue
        if len(f) >= 12 and f[1].isdigit() and f[6].isdigit():                     # PAF
            contig, chrom, s, e, mapq = f[0], f[5], int(f[7]), int(f[8]), int(f[11])
            primary = not any(t.startswith("tp:A:S") for t in f[12:])
        else:                                                                      # the BAM-derived table
            contig, chrom, s, e, mapq, primary = f[2], f[1], int(f[3]) - 1, int(f[4]) - 1, int(f[6]), f[7] != "secondary"
        if primary and mapq >= min_mapq:
            yield label, chrom, contig, s, e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="label=path")
    ap.add_argument("--chroms", default="chrX,chrY")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    chroms = a.chroms.split(",")
    per = {}
    for item in a.inputs:
        label, path = item.split("=", 1)
        for lab, chrom, contig, s, e in read(label, path, a.min_mapq):
            if chrom not in chroms:
                continue
            for x, y in clip(chrom, s, e):
                d = per.setdefault((lab, chrom), dict(intervals=[], contigs={}))
                d["intervals"].append((x, y))
                d["contigs"][contig] = d["contigs"].get(contig, 0) + (y - x)
    with open(a.out, "w") as fh:
        fh.write("haplotype\tchrom\tlength_non_par\taligned_bases\tcovered_bases\taligned_over_covered\tcontigs\tlargest_contig_bases\tlargest_contig\n")
        for (lab, chrom), d in sorted(per.items()):
            iv = sorted(d["intervals"])
            covered, cur = 0, None
            for x, y in iv:
                if cur is None or x > cur[1]:
                    if cur:
                        covered += cur[1] - cur[0]
                    cur = [x, y]
                else:
                    cur[1] = max(cur[1], y)
            if cur:
                covered += cur[1] - cur[0]
            aligned = sum(y - x for x, y in iv)
            L = LENGTH.get(chrom, 0) - sum(b - a_ for a_, b in PAR.get(chrom, []))
            big = max(d["contigs"].items(), key=lambda kv: kv[1])
            fh.write("%s\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%s\n" % (lab, chrom, L, aligned, covered, aligned / covered if covered else float("nan"), len(d["contigs"]), big[1], big[0]))
    print(open(a.out).read(), end="", file=sys.stderr)


if __name__ == "__main__":
    main()
