"""A reference panel: the depth structure every genome shares along the reference - centromere flanks, segmental
duplications, the acrocentric short arms, the caller's own depth filter - as the median LRR per bin over other genomes,
with its spread. Subtracted from each member's LRR it leaves the member's own events; bins the panel cannot pin (a
wide spread, too few genomes) are left out of the calls. Built from any genomes counted the same way: per-sample or
multi-sample VCFs (`triokaryo panel --vcfs`), or earlier runs' bins.tsv (`--runs`), with the cohort itself the natural
panel when there are many trios."""
import csv
import os

import numpy as np

from .model import NA
from .vcfscan import scan_vcf


def sample_profile(sites_by_chrom, genome, bin_size, member=0, min_sites=20, min_dp=8, min_gq=20, min_het=5):
    """One sample's self-normalised bin LRR (log2 of the bin's median depth over the autosomal median), the X aligned to
    two copies by the sample's own X median, and its band deviation per bin (the median |BAF - 1/2| at heterozygous sites):
    and its heterozygosity rate per bin (heterozygous over confident calls):
    ({(chrom, start): lrr}, {(chrom, start): bdev}, {(chrom, start): het_rate})."""
    out, depth, bdev, hetr = {}, {}, {}, {}
    for chrom, s in sites_by_chrom.items():
        L = genome.length[chrom]
        nb = int(np.ceil(L / bin_size))
        idx = np.clip(s.pos // bin_size, 0, nb - 1)
        dp, alt, gt, gq = s.dp[member], s.alt[member], s.gt[member], s.gq[member]
        called = (dp >= min_dp) & (gq >= min_gq) & (gt >= 0)
        het = called & (gt == 1)
        for b in np.unique(idx):
            sel = (idx == b) & (dp > 0) & ~s.par
            if sel.sum() >= min_sites:
                depth[(chrom, int(b) * bin_size)] = float(np.median(dp[sel]))
            h = (idx == b) & het & ~s.par
            if h.sum() >= min_het:
                bdev[(chrom, int(b) * bin_size)] = float(np.median(np.abs(alt[h] / dp[h] - 0.5)))
            c_ = (idx == b) & called & ~s.par
            if c_.sum() >= min_sites:
                hetr[(chrom, int(b) * bin_size)] = float(h.sum() / c_.sum())
    auto = [v for (c, _), v in depth.items() if c in genome.autosomes]
    if len(auto) < 100:
        return out, bdev, hetr
    med = float(np.median(auto))
    for key, v in depth.items():
        out[key] = float(np.log2(v / med))
    xs = [v for (c, _), v in out.items() if c == "chrX"]
    if len(xs) >= 20:
        xmed = float(np.median(xs))
        shift = -xmed if xmed < -0.5 else 0.0            # a male's X read as two copies
        for key in list(out):
            if key[0] == "chrX":
                out[key] += shift
    return out, bdev, hetr


def build_panel(genome, bin_size, vcfs=(), samples=None, runs=(), log=None, thin=1):
    """Per bin: n genomes, the median LRR and its robust SD. vcfs: VCF paths (every sample of each unless `samples` names
    them); runs: directories of earlier runs (bins.tsv + summary.tsv)."""
    profiles = []
    for v in vcfs:
        import pysam
        names = list(pysam.VariantFile(v).header.samples)
        use = [n for n in names if (not samples or n in samples)]
        if not use:
            continue
        scan = scan_vcf(v, use, genome, thin=thin, log=None)
        for m, name in enumerate(use):
            lr, bd, hr = sample_profile(scan.chroms, genome, bin_size, m)
            profiles.append((name, lr, bd, hr))
            if log:
                log("panel: %s from %s (%d bins)" % (name, os.path.basename(v), len(lr)))
    for d in runs:
        bpath, spath = os.path.join(d, "bins.tsv"), os.path.join(d, "summary.tsv")
        if not (os.path.exists(bpath) and os.path.exists(spath)):
            continue
        summ = next(csv.DictReader(open(spath), delimiter="\t"))
        if int(summ.get("bin_size", bin_size)) != bin_size:
            continue
        rows = list(csv.DictReader(open(bpath), delimiter="\t"))
        for role in ("child", "father", "mother"):
            name = summ.get(role, role)
            try:
                xc = float(summ.get("%s_x_copies" % role, "NA"))
            except ValueError:
                xc = NA
            shift = float(np.log2(2.0 / xc)) if np.isfinite(xc) and xc > 0 else 0.0
            prof, bd, hr = {}, {}, {}
            for r in rows:
                key = (r["chrom"], int(r["start"]))
                v = r.get("%s_lrr" % role, "NA")
                if v != "NA":
                    prof[key] = float(v) + (shift if r["chrom"] == "chrX" else 0.0)
                b = r.get("%s_bdev" % role, "NA")
                if b != "NA":
                    bd[key] = float(b)
                h = r.get("%s_het_rate" % role, "NA")
                if h != "NA" and int(r.get("%s_n_called" % role, 0) or 0) >= 20:
                    hr[key] = float(h)
            profiles.append((name, prof, bd, hr))
            if log:
                log("panel: %s from %s (%d bins)" % (name, d, len(prof)))
    keys = sorted({k for _, p, _, _ in profiles for k in p} | {k for _, _, b, _ in profiles for k in b},
                  key=lambda k: (genome.chroms.index(k[0]) if k[0] in genome.chroms else 99, k[1]))
    out = []
    for key in keys:
        v = np.array([p[key] for _, p, _, _ in profiles if key in p])
        v = v[np.isfinite(v)]
        b = np.array([bd[key] for _, _, bd, _ in profiles if key in bd])
        b = b[np.isfinite(b)]
        h = np.array([hr[key] for _, _, _, hr in profiles if key in hr])
        h = h[np.isfinite(h)]
        if not len(v) and not len(b):
            continue
        med = float(np.median(v)) if len(v) else NA
        rsd = float(1.4826 * np.median(np.abs(v - med))) if len(v) > 1 else NA
        out.append(dict(chrom=key[0], start=key[1], end=min(key[1] + bin_size, genome.length[key[0]]), n=int(len(v)), lrr_median=med, lrr_rsd=rsd,
                        n_bdev=int(len(b)), bdev_median=float(np.median(b)) if len(b) else NA, n_het=int(len(h)), het_rate_median=float(np.median(h)) if len(h) else NA))
    return out, [n for n, _, _, _ in profiles]


def write_panel(path, rows, samples, bin_size):
    with open(path, "w") as fh:
        fh.write("#triokaryo panel: bin %d, %d genomes: %s\n" % (bin_size, len(samples), ",".join(samples)))
        fh.write("chrom\tstart\tend\tn\tlrr_median\tlrr_rsd\tn_bdev\tbdev_median\tn_het\thet_rate_median\n")
        f = lambda x: ("%.4f" % x) if np.isfinite(x) else "NA"
        for r in rows:
            fh.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%d\t%s\n" % (r["chrom"], r["start"], r["end"], r["n"], f(r["lrr_median"]), f(r["lrr_rsd"]), r["n_bdev"],
                                                                      f(r["bdev_median"]), r["n_het"], f(r["het_rate_median"])))


def load_panel(path):
    """{(chrom, start): (n, median, rsd, n_bdev, bdev_median, n_het, het_rate_median)}"""
    out = {}
    fl = lambda x: float(x) if x != "NA" else NA
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom\t") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            c, s, e, n, med, rsd = parts[:6]
            nb, bmed = (parts[6], parts[7]) if len(parts) >= 8 else ("0", "NA")
            nh, hmed = (parts[8], parts[9]) if len(parts) >= 10 else ("0", "NA")
            out[(c, int(s))] = (int(n), fl(med), fl(rsd), int(nb), fl(bmed), int(nh), fl(hmed))
    return out
