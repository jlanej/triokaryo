#!/usr/bin/env python3
"""DRAGEN's own copy-number calls for a sample, in the generic events form triokaryo draws and matches (`--events`): the segmental
CNV calls (cnv.vcf.gz: PASS DRAGEN:GAIN and DRAGEN:LOSS records of at least --min-mb, with the integer copy number CN and the
segment mean SM) and the whole-chromosome ploidy calls (ploidy.vcf.gz: <DUP> and <DEL> records, with the chromosome's copies as
2 DC over the autosomal median DC).

  dragen_events.py SAMPLE.cnv.vcf.gz SAMPLE.ploidy.vcf.gz [...] [--min-mb=2] [--min-dev=0.3] > events.external.tsv

Several samples' files may be listed; the sample name is the VCF's. An autosomal ploidy call is kept when the chromosome's
copies are at least --min-dev from two (DRAGEN flags the GC-rich chromosomes 19 and 22 at 2.1 copies in many libraries);
the X and the Y are kept as DRAGEN flagged them. Writes sample, chrom, start, end, label, type."""
import gzip
import sys
from statistics import median

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}


def records(path):
    sample = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                sample = f[9]
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in f[7].split(";"))
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            yield sample, f[0], int(f[1]), f[2], f[4], f[6], info, fmt


def cnv_rows(path, min_bp):
    for sample, chrom, pos, vid, alt, filt, info, fmt in records(path):
        if filt != "PASS" or not (vid.startswith("DRAGEN:GAIN") or vid.startswith("DRAGEN:LOSS")):
            continue
        end = int(info["END"])
        if end - pos + 1 < min_bp:
            continue
        kind = "gain" if "GAIN" in vid else "loss"
        yield sample, chrom, pos - 1, end, "DRAGEN CNV: CN %s (segment mean %.2f)" % (fmt.get("CN", "?"), float(fmt.get("SM", "nan"))), kind


def ploidy_rows(path, min_dev=0.3):
    recs = list(records(path))
    dc = {chrom: float(fmt["DC"]) for _, chrom, _, _, _, _, _, fmt in recs if "DC" in fmt}
    auto = [v for c, v in dc.items() if c in AUTOSOMES]
    med = median(auto) if auto else float("nan")
    for sample, chrom, pos, vid, alt, filt, info, fmt in recs:
        if alt not in ("<DUP>", "<DEL>"):
            continue
        copies = 2.0 * dc[chrom] / med
        if chrom in AUTOSOMES and abs(copies - 2.0) < min_dev:
            continue
        yield sample, chrom, 0, int(info["END"]), "DRAGEN ploidy: %s, %.2f copies (%.2f of the expected depth%s)" % (
            alt.strip("<>"), copies, float(fmt.get("NDC", "nan")), "" if filt == "PASS" else "; " + filt), "gain" if alt == "<DUP>" else "loss"


def main(argv):
    min_mb, min_dev = 2.0, 0.3
    paths = [a for a in argv if not a.startswith("--")]
    for a in argv:
        if a.startswith("--min-mb="):
            min_mb = float(a.split("=", 1)[1])
        if a.startswith("--min-dev="):
            min_dev = float(a.split("=", 1)[1])
    print("sample\tchrom\tstart\tend\tlabel\ttype")
    for p in paths:
        rows = cnv_rows(p, int(min_mb * 1e6)) if ".cnv." in p else ploidy_rows(p, min_dev) if ".ploidy." in p else ()
        for r in rows:
            print("\t".join(str(x) for x in r))


if __name__ == "__main__":
    main(sys.argv[1:])
