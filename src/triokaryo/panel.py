"""A reference panel: the depth structure shared by every genome along the reference (centromere flanks, segmental
duplications, acrocentric short arms, the caller's depth filter) as the per-bin median LRR over other genomes with its
spread, together with the per-bin median band deviation and heterozygosity rate. Subtracted from each member's tracks
it leaves the member's own events; bins the panel cannot characterise (wide spread, too few genomes) are excluded from
calling. Built from any genomes processed the same way: per-sample or multi-sample VCFs (`triokaryo panel --vcfs`), or
earlier runs' bins.tsv (`--runs`); with many trios the cohort itself is the natural panel."""
import csv
import os

import numpy as np

from .model import NA, homref_depth_ratio, reference_block_depth, trimmed_mean
from .vcfscan import scan_vcf


NO_Y_LEVEL = -1.5           # a genome whose Y bins' median LRR (diploid scale) lies under this carries no Y: what depth its Y bins hold
                            # is stray reads (X-homologous sequence in a female), not a chromosome, and gives the Y rows nothing


def align_sex_chromosomes(lrr, y_copies=None):
    """A genome's X and Y LRR on the diploid scale, in place: a one-copy X (median under -0.5: a male, or a 45,X) and a one-copy
    Y (median under -0.5) are each shifted up by exactly one unit. The rows then hold each chromosome's level with its
    mappability deficit whatever the genome's sex, so that the deficit cancels for a member: shifting a male's X by its own
    median instead would erase the deficit for males only, and a male-heavy panel would read every female's X as a loss.
    A genome without a Y contributes no Y rows: its Y bins, where a joint call gives a female depth from stray X-homologous
    reads, are dropped - by the run's own reading where given (y_copies 0), else by their level (median under NO_Y_LEVEL).
    Before this, every female with a few Y bins counted as a genome with a Y and her stray depth, shifted up a unit as if
    a one-copy Y, went into the Y rows: a cohort of 663 read as 663 with a Y."""
    ys = [k for k in lrr if k[0] == "chrY"]
    if ys:
        v = [lrr[k] for k in ys]
        if (y_copies is not None and y_copies == 0) or float(np.median(v)) < NO_Y_LEVEL:
            for k in ys:
                del lrr[k]
    for c, n_min in (("chrX", 20), ("chrY", 3)):
        v = [x for (cc, _), x in lrr.items() if cc == c]
        if len(v) >= n_min and float(np.median(v)) < -0.5:
            for key in list(lrr):
                if key[0] == c:
                    lrr[key] += 1.0
    return lrr


def sample_profile(sites_by_chrom, genome, bin_size, member=0, min_sites=20, min_dp=8, min_gq=20, min_het=5, depth_sites="all"):
    """One sample's self-normalised bin LRR (log2 of the bin's trimmed-mean depth over the autosomal median), the X and Y on the
    diploid scale (align_sex_chromosomes), its band deviation per bin (the median |BAF - 1/2| at heterozygous sites) and its
    heterozygosity rate per bin (heterozygous over confident calls):
    ({(chrom, start): lrr}, {(chrom, start): bdev}, {(chrom, start): het_rate})."""
    out, depth, bdev, hetr = {}, {}, {}, {}
    for chrom, s in sites_by_chrom.items():
        L = genome.length[chrom]
        nb = int(np.ceil(L / bin_size))
        idx = np.clip(s.pos // bin_size, 0, nb - 1)
        dp, alt, gt, gq = s.dp[member], s.alt[member], s.gt[member], s.gq[member]
        called = (dp >= min_dp) & (gq >= min_gq) & (gt >= 0)
        het = called & (gt == 1)
        with_depth = (dp > 0) & ((gt != 0) if depth_sites == "variant" else True)      # a joint caller's hom-ref DP may be a block's minimum
        for b in np.unique(idx):
            sel = (idx == b) & with_depth & ~s.par
            if sel.sum() >= min_sites:
                depth[(chrom, int(b) * bin_size)] = trimmed_mean(dp[sel])
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
    return align_sex_chromosomes(out), bdev, hetr


def build_panel(genome, bin_size, vcfs=(), samples=None, runs=(), log=None, thin=1, roles=None):
    """Per bin: n genomes, the median LRR and its robust SD, the median band deviation and heterozygosity rate. vcfs: VCF paths
    (every sample of each unless `samples` names them); runs: directories of earlier runs (bins.tsv + summary.tsv), taking the
    members named in `roles` (default all three). Returns (rows, sample names, genomes with Y rows). The Y rows come from the
    genomes with a Y (the males), so a cohort's fathers and sons are the Y reference for a father without a son."""
    profiles = []
    for v in vcfs:
        import pysam
        names = list(pysam.VariantFile(v).header.samples)
        use = [n for n in names if (not samples or n in samples)]
        if not use:
            continue
        scan = scan_vcf(v, use, genome, thin=thin, log=None)
        ds = "variant" if reference_block_depth(scan, homref_depth_ratio(scan, genome)) else "all"
        for m, name in enumerate(use):
            lr, bd, hr = sample_profile(scan.chroms, genome, bin_size, m, depth_sites=ds)
            profiles.append((name, lr, bd, hr))
            if log:
                log("panel: %s from %s (%d bins%s)" % (name, os.path.basename(v), len(lr), "; depth from the genome's own variant genotypes" if ds == "variant" else ""))
    for d in runs:
        bpath, spath = os.path.join(d, "bins.tsv"), os.path.join(d, "summary.tsv")
        if not (os.path.exists(bpath) and os.path.exists(spath)):
            continue
        summ = next(csv.DictReader(open(spath), delimiter="\t"))
        if int(summ.get("bin_size", bin_size)) != bin_size:
            continue
        rows = list(csv.DictReader(open(bpath), delimiter="\t"))
        for role in (roles or ("child", "father", "mother")):
            name = summ.get(role, role)
            yc = summ.get("%s_y_copies" % role, "")                   # the run's own reading of the member's Y copies, where it has one
            try:
                yc = int(float(yc)) if yc not in ("", "NA", "nan") else None
            except ValueError:
                yc = None
            prof, bd, hr = {}, {}, {}
            for r in rows:
                key = (r["chrom"], int(r["start"]))
                v = r.get("%s_lrr" % role, "NA")
                if v != "NA":
                    prof[key] = float(v)
                b = r.get("%s_bdev" % role, "NA")
                if b != "NA":
                    bd[key] = float(b)
                h = r.get("%s_het_rate" % role, "NA")
                if h != "NA" and int(r.get("%s_n_called" % role, 0) or 0) >= 20:
                    hr[key] = float(h)
            profiles.append((name, align_sex_chromosomes(prof, y_copies=yc), bd, hr))
            if log:
                log("panel: %s from %s (%d bins%s)" % (name, d, len(prof), "" if any(k[0] == "chrY" for k in prof) else "; no Y"))
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
    n_y = sum(1 for _, p, _, _ in profiles if sum(1 for k in p if k[0] == "chrY") >= 3)
    return out, [n for n, _, _, _ in profiles], n_y


def write_panel(path, rows, samples, bin_size):
    with open(path, "w") as fh:
        fh.write("#triokaryo panel: bin %d, %d genomes: %s\n" % (bin_size, len(samples), ",".join(samples)))
        fh.write("chrom\tstart\tend\tn\tlrr_median\tlrr_rsd\tn_bdev\tbdev_median\tn_het\thet_rate_median\n")
        f = lambda x: ("%.4f" % x) if np.isfinite(x) else "NA"
        for r in rows:
            fh.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%d\t%s\n" % (r["chrom"], r["start"], r["end"], r["n"], f(r["lrr_median"]), f(r["lrr_rsd"]), r["n_bdev"],
                                                                      f(r["bdev_median"]), r["n_het"], f(r["het_rate_median"])))


def load_panel(path):
    """{(chrom, start): (n, median, rsd, n_bdev, bdev_median, n_het, het_rate_median)}. A panel written before the Y was
    put on the diploid scale holds its Y rows at the male's one-copy level (median under -0.5): they are shifted up by one
    unit here, so that `triokaryo panel` output of either vintage reads the same."""
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
    ys = [v[1] for (c, _), v in out.items() if c == "chrY" and v[0] >= 3 and np.isfinite(v[1])]
    if len(ys) >= 3 and float(np.median(ys)) < -0.5:
        for key in list(out):
            if key[0] == "chrY":
                v = out[key]
                out[key] = (v[0], v[1] + 1.0 if np.isfinite(v[1]) else v[1]) + v[2:]
    return out
