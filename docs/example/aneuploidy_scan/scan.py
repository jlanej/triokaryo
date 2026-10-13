#!/usr/bin/env python3
"""Whole-chromosome aneuploidies across the 1000 Genomes DRAGEN cohort from the per-sample ploidy VCFs, with the trio and the
assembly each sample has: candidates for a trio analysis whose reading can be held against a long-read assembly.

  scan.py --ploidy-dir ploidy/ --ped 20130606_g1k_3202_samples_ped_population.txt --hprc-y1 Year1_assemblies_v2_genbank.index \
          --hprc-r2 assemblies_release2_v1.0.index.csv --hgsvc3 supp-table_accessions-assembly-data.2024v2.0.tsv --out .

A DRAGEN ploidy VCF holds one record per chromosome with DC (depth of coverage) and NDC (depth normalised to the expected
ploidy). Copies per chromosome are taken as 2 DC / median DC over the autosomes, so that the X and the Y are read on the same
scale whatever sex DRAGEN assumed. Writes aneuploidies.tsv (every sample with a whole-chromosome deviation), candidates.tsv
(those with both parents in the cohort and an assembly) and summary.md."""
import argparse
import csv
import glob
import gzip
import os
import re
from statistics import median

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]
GC_RICH = {"chr16", "chr17", "chr19", "chr20", "chr22"}   # a library's GC bias lifts these a few per cent: a mosaic gain there needs the trio analysis
GAIN, MOSAIC = 2.5, 2.2        # copies: a trisomy in at least half the cells; a mosaic gain in at least a fifth
LOSS, MOSAIC_LOSS = 1.5, 1.8   # a monosomy in at least half the cells; a mosaic loss in at least a fifth
SEX_OFF = 0.25                 # the X or the Y this far from an integer: a mosaic


def read_ploidy(path):
    """{chrom: (DC, NDC, FILTER, ALT)}"""
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            keys, vals = f[8].split(":"), f[9].split(":")
            d = dict(zip(keys, vals))
            out[f[0]] = (float(d["DC"]), float(d.get("NDC", "nan")), f[6], f[4])
    return out


def copies(p):
    auto = [p[c][0] for c in AUTOSOMES if c in p]
    med = median(auto) if auto else float("nan")
    return {c: 2.0 * v[0] / med for c, v in p.items()}, med


def sex_complement(cp):
    x, y = cp.get("chrX", float("nan")), cp.get("chrY", float("nan"))
    rx, ry = int(round(x)), int(round(y)) if y == y else 0
    comp = "X" * max(rx, 0) + "Y" * max(ry, 0)
    mosaic = abs(x - rx) >= SEX_OFF or (y == y and abs(y - ry) >= SEX_OFF)
    return comp, mosaic, x, y


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ploidy-dir", required=True)
    ap.add_argument("--ped", required=True, help="the 1000 Genomes 3,202-sample pedigree with populations")
    ap.add_argument("--hprc-y1", help="HPRC year-1 assembly index")
    ap.add_argument("--hprc-r2", help="HPRC release-2 assembly index (csv)")
    ap.add_argument("--hgsvc3", help="HGSVC3 assembly accessions table")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    ped = {}
    for line in open(a.ped):
        f = line.split()
        if len(f) < 7 or f[0] == "FamilyID":
            continue
        ped[f[1]] = dict(family=f[0], father=f[2], mother=f[3], sex={"1": "M", "2": "F"}.get(f[4], ""), pop=f[5], superpop=f[6])
    children = {s for s, r in ped.items() if r["father"] in ped and r["mother"] in ped}
    parent_of = {}
    for s in children:
        parent_of.setdefault(ped[s]["father"], []).append(s)
        parent_of.setdefault(ped[s]["mother"], []).append(s)
    asm = {}
    if a.hprc_y1:
        for line in open(a.hprc_y1):
            s = line.split("\t")[0]
            if re.fullmatch(r"(HG|NA)\d{5}", s):
                asm.setdefault(s, set()).add("HPRC-y1 hifiasm trio")
    if a.hprc_r2:
        for r in csv.DictReader(open(a.hprc_r2)):
            asm.setdefault(r["sample_id"], set()).add("HPRC-r2 %s %s" % (r["assembly_method"], r["phasing"]))
    if a.hgsvc3:
        for s in set(re.findall(r"(?:HG|NA)\d{5}", open(a.hgsvc3).read())):
            asm.setdefault(s, set()).add("HGSVC3 verkko")
    rows = []
    n_samples = 0
    for path in sorted(glob.glob(os.path.join(a.ploidy_dir, "*.ploidy.vcf.gz"))):
        sample = os.path.basename(path).split(".")[0]
        p = read_ploidy(path)
        if not p:
            continue
        n_samples += 1
        cp, med = copies(p)
        flags = []
        for c in AUTOSOMES:
            v = cp.get(c)
            if v is None:
                continue
            dragen = ", DRAGEN %s" % p[c][3].strip("<>") if p[c][3] not in (".", "") else ""
            if v >= GAIN:
                flags.append("%s trisomy (%.2f copies%s)" % (c, v, dragen))
            elif v >= MOSAIC:
                flags.append("%s mosaic gain (%.2f copies%s%s)" % (c, v, dragen, "; GC-rich, possibly the library's GC bias" if c in GC_RICH else ""))
            elif v <= LOSS:
                flags.append("%s monosomy (%.2f copies)" % (c, v))
            elif v <= MOSAIC_LOSS:
                flags.append("%s mosaic loss (%.2f copies)" % (c, v))
        comp, mosaic, x, y = sex_complement(cp)
        info = ped.get(sample, {})
        expected = {"M": "XY", "F": "XX"}.get(info.get("sex", ""), "")
        if expected and comp != expected:
            flags.append("sex chromosomes %s against the pedigree's %s (X %.2f, Y %.2f copies)" % (comp or "none", expected, x, y))
        elif mosaic:
            flags.append("mosaic sex chromosomes (X %.2f, Y %.2f copies)" % (x, y))
        if not flags:
            continue
        role = "child of %s x %s" % (info["father"], info["mother"]) if sample in children else ("parent of %s" % ",".join(parent_of[sample]) if sample in parent_of else "no trio")
        rows.append(dict(sample=sample, sex=info.get("sex", ""), population=info.get("pop", ""), superpopulation=info.get("superpop", ""), trio=role,
                         flags="; ".join(flags), x_copies="%.2f" % x, y_copies="%.2f" % y, autosomal_depth="%.1f" % med,
                         assembly="; ".join(sorted(asm.get(sample, ()))) or "none"))
    cols = ["sample", "sex", "population", "superpopulation", "trio", "flags", "x_copies", "y_copies", "autosomal_depth", "assembly"]
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "aneuploidies.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    cand = [r for r in rows if r["trio"].startswith("child") and r["assembly"] != "none"]
    with open(os.path.join(a.out, "candidates.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in cand:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    kinds = dict(trisomy=sum("trisomy" in r["flags"] for r in rows), mosaic_gain=sum("mosaic gain" in r["flags"] for r in rows),
                 loss=sum("monosomy" in r["flags"] or "mosaic loss" in r["flags"] for r in rows), sex=sum("sex chromosomes" in r["flags"] for r in rows))
    with open(os.path.join(a.out, "summary.md"), "w") as fh:
        fh.write("# Whole-chromosome aneuploidies in the 1000 Genomes DRAGEN cohort\n\n")
        fh.write("%d samples read; %d with a whole-chromosome deviation: %d with an autosomal trisomy (at least half the cells), %d with a mosaic "
                 "autosomal gain (a fifth to a half), %d with an autosomal loss, %d with a sex-chromosome complement other than the pedigree's.\n\n"
                 % (n_samples, len(rows), kinds["trisomy"], kinds["mosaic_gain"], kinds["loss"], kinds["sex"]))
        fh.write("%d of them are children of a trio in the cohort, %d have an assembly (HPRC year 1 or release 2, HGSVC3), %d both: the candidates.\n\n"
                 % (sum(r["trio"].startswith("child") for r in rows), sum(r["assembly"] != "none" for r in rows), len(cand)))
        fh.write("| sample | sex | population | trio | aneuploidy | assembly |\n| --- | --- | --- | --- | --- | --- |\n")
        for r in cand:
            fh.write("| %s | %s | %s | %s | %s | %s |\n" % (r["sample"], r["sex"], r["population"], r["trio"], r["flags"], r["assembly"]))
    print("%d samples, %d aneuploid, %d candidates -> %s" % (n_samples, len(rows), len(cand), a.out))


if __name__ == "__main__":
    main()
