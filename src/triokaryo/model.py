"""From sites to tracks: fixed-width bins along each chromosome; per member the depth (median over the bin's sites), the
LRR (log2 of the depth over the member's autosomal median; GC- and panel-corrected where given), the heterozygosity
rate and the B-allele band deviation at heterozygous sites; the within-trio depth tracks (child over the parents' mean,
father over mother, per site); and the band-deviation estimator used by the calls.

Expected values under an event in a cell fraction f:
  depth (LRR)   one-copy gain: log2(1 + f/2); one-copy loss: log2(1 - f/2)
  BAF bands     at a heterozygous site the alt fraction is 1/2 +- d: gain d = f / (2 (2 + f)) (1/3 and 2/3 at f = 1),
                loss d = f / (2 (2 - f)), copy-neutral LOH d = f / 2
  het rate      a constitutional loss of heterozygosity (isodisomy, deletion, run of homozygosity) leaves no heterozygous calls
The depth gives the copy number, the bands an independent estimate of f, the heterozygosity rate the constitutional case.
"""
from dataclasses import dataclass

import numpy as np

from .vcfscan import GT_HET, GT_HOMALT, GT_HOMREF

MEMBERS = ("child", "father", "mother")
NA = float("nan")


@dataclass
class Bins:
    chrom: list
    start: np.ndarray
    end: np.ndarray
    gc: np.ndarray
    n_dp: np.ndarray
    depth: np.ndarray
    n_called: np.ndarray
    n_het: np.ndarray
    het_rate: np.ndarray
    bdev: np.ndarray
    lrr: np.ndarray
    lrr_gc: np.ndarray
    child_vs_mid: np.ndarray
    father_vs_mother: np.ndarray
    index: dict                 # chrom -> (first bin, one past the last)
    bin_size: int
    autosomal: np.ndarray       # per bin
    par: np.ndarray             # per bin: mostly pseudoautosomal (X, Y)
    panel_median: np.ndarray = None   # the reference panel's LRR per bin (NaN without a panel)
    panel_rsd: np.ndarray = None
    masked: np.ndarray = None         # bins the panel cannot characterise: excluded from calling
    panel_bdev: np.ndarray = None     # the panel's band deviation per bin
    masked_bands: np.ndarray = None   # bins whose bands are split across the panel's genomes (paralogy): excluded from the LOH search
    bdev_adj: np.ndarray = None       # (3, B) the band deviation net of the panel's regional excess (the bands track the LOH search segments)
    het_rel: np.ndarray = None        # (3, B) the heterozygosity rate over the panel's (the rate track the LOH search segments)

    def of(self, chrom):
        a, b = self.index.get(chrom, (0, 0))
        return slice(a, b)

    @property
    def n(self):
        return len(self.start)


def load_gc_track(path):
    """chrom, start (0- or 1-based: the bin's left edge), end, gc per line -> {(chrom, bin_start_0based): gc}."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "chrom\t")):
                continue
            c, s, e, g = line.rstrip("\n").split("\t")[:4]
            try:
                out[(c, int(s))] = float(g)
            except ValueError:
                continue
    return out


def _median_or_nan(v):
    return float(np.median(v)) if len(v) else NA


PANEL_MAX_RSD = 0.25        # a bin whose LRR robust SD across the panel's genomes exceeds this is not called
PANEL_MIN_N = 5             # minimum genomes per bin: a median over fewer follows one genome's own event
PANEL_MIN_N_Y = 3           # the same for the Y, to which only the panel's males contribute
PANEL_BDEV_EXCESS = 0.03    # a bin whose panel band deviation exceeds the panel's genome-wide median by this much (paralogous sequence,
                            # where the alleles of two loci are counted as one) is excluded from the loss-of-heterozygosity search


def make_bins(scan, genome, bin_size=1_000_000, min_dp=8, min_gq=20, gc_track=None, min_het=5, panel=None):
    chroms = [c for c in genome.chroms if c in scan.chroms]
    rows = []
    index = {}
    for c in chroms:
        s = scan.sites(c)
        L = genome.length[c]
        nb = int(np.ceil(L / bin_size))
        a = len(rows)
        edges = np.arange(0, (nb + 1) * bin_size, bin_size)
        idx = np.clip(np.searchsorted(edges, s.pos, side="right") - 1, 0, nb - 1)
        order = np.argsort(idx, kind="stable")
        idx_s = idx[order]
        bounds = np.searchsorted(idx_s, np.arange(nb + 1))
        for b in range(nb):
            sel = order[bounds[b]:bounds[b + 1]]
            rows.append((c, b * bin_size, min((b + 1) * bin_size, L), sel))
        index[c] = (a, len(rows))
    B = len(rows)
    chrom = [r[0] for r in rows]
    start = np.array([r[1] for r in rows], dtype=np.int64)
    end = np.array([r[2] for r in rows], dtype=np.int64)
    gc = np.full(B, NA)
    if gc_track:
        for i, r in enumerate(rows):
            g = gc_track.get((r[0], int(r[1])))
            if g is None:
                g = gc_track.get((r[0], int(r[1]) + 1))
            if g is not None:
                gc[i] = g
    n_dp = np.zeros((3, B), dtype=np.int64)
    depth = np.full((3, B), NA)
    n_called = np.zeros((3, B), dtype=np.int64)
    n_het = np.zeros((3, B), dtype=np.int64)
    bdev = np.full((3, B), NA)
    child_vs_mid = np.full(B, NA)
    father_vs_mother = np.full(B, NA)
    par = np.zeros(B, dtype=bool)
    for i, (c, b0, b1, sel) in enumerate(rows):
        if not len(sel):
            continue
        s = scan.sites(c)
        dp = s.dp[:, sel]
        alt = s.alt[:, sel]
        gt = s.gt[:, sel]
        gq = s.gq[:, sel]
        par[i] = bool(s.par[sel].any())                            # a bin touching a pseudoautosomal region: a male is diploid there
        for m in range(3):
            ok = dp[m] > 0
            n_dp[m, i] = int(ok.sum())
            depth[m, i] = _median_or_nan(dp[m][ok])
            called = (dp[m] >= min_dp) & (gq[m] >= min_gq) & (gt[m] >= GT_HOMREF)
            n_called[m, i] = int(called.sum())
            het = called & (gt[m] == GT_HET)
            n_het[m, i] = int(het.sum())
            if n_het[m, i] >= min_het:
                baf = alt[m][het] / np.maximum(dp[m][het], 1)
                bdev[m, i] = float(np.median(np.abs(baf - 0.5)))
        all3 = (dp > 0).all(axis=0)
        if all3.sum() >= 5:
            d = dp[:, all3].astype(float)
            child_vs_mid[i] = float(np.median(np.log2(d[0] / ((d[1] + d[2]) / 2.0))))
            father_vs_mother[i] = float(np.median(np.log2(d[1] / d[2])))
    autosomal = np.array([c in genome.autosomes for c in chrom])
    with np.errstate(divide="ignore", invalid="ignore"):
        het_rate = np.where(n_called > 0, n_het / np.maximum(n_called, 1), NA)
    lrr = np.full((3, B), NA)
    lrr_gc = np.full((3, B), NA)
    pmed, prsd, masked = np.full(B, NA), np.full(B, NA), np.zeros(B, dtype=bool)
    pbdev, masked_bands, phet = np.full(B, NA), np.zeros(B, dtype=bool), np.full(B, NA)
    if panel:
        typical = np.nanmedian([e[4] for e in panel.values() if len(e) > 4 and np.isfinite(e[4])]) if any(len(e) > 4 and np.isfinite(e[4]) for e in panel.values()) else NA
        for i, r in enumerate(rows):
            e = panel.get((r[0], int(r[1])))
            if e is None:
                masked[i] = True                               # a bin the panel never saw
                continue
            n_, med_, rsd_ = e[:3]
            pmed[i], prsd[i] = med_, rsd_
            masked[i] = n_ < (PANEL_MIN_N_Y if r[0] == "chrY" else PANEL_MIN_N) or not np.isfinite(med_) or (np.isfinite(rsd_) and rsd_ > PANEL_MAX_RSD)
            if len(e) > 4 and np.isfinite(e[4]):
                pbdev[i] = e[4]
                masked_bands[i] = bool(np.isfinite(typical) and e[4] > typical + PANEL_BDEV_EXCESS)
            if len(e) > 6 and np.isfinite(e[6]) and e[5] >= PANEL_MIN_N:
                phet[i] = e[6]
    for m in range(3):
        ok = autosomal & np.isfinite(depth[m]) & (n_dp[m] >= 20)
        med = float(np.median(depth[m][ok])) if ok.any() else NA
        with np.errstate(divide="ignore", invalid="ignore"):
            lrr[m] = np.where(np.isfinite(depth[m]) & (depth[m] > 0), np.log2(depth[m] / med), NA)
        base = lrr[m].copy()
        if panel:
            base = np.where(np.isfinite(pmed) & ~masked, base - pmed, NA)       # the shared structure removed; masked bins excluded from calling
            fit = autosomal & (n_dp[m] >= 20) & np.isfinite(base)
            if fit.any():                                                      # re-centred: the panel's zero is not this genome's
                base[np.isfinite(base)] -= float(np.median(base[fit]))
        lrr_gc[m] = gc_correct(base, gc, autosomal & (n_dp[m] >= 20)) if np.isfinite(gc).any() else base
    bdev_adj, het_rel = bdev.copy(), het_rate.copy()
    if panel:                                                      # bands split in every genome: not a member's loss of heterozygosity;
        bdev[:, masked_bands] = NA                                 # the panel's regional band excess and regional heterozygosity removed
        with np.errstate(invalid="ignore", divide="ignore"):
            reg = np.where(np.isfinite(pbdev), pbdev - (typical if np.isfinite(typical) else 0.0), 0.0)
            bdev_adj = np.where(np.isfinite(bdev), bdev - reg[None, :], NA)
            het_rel = np.where(np.isfinite(het_rate) & np.isfinite(phet) & (phet > 0), het_rate / phet, NA)
    return Bins(chrom, start, end, gc, n_dp, depth, n_called, n_het, het_rate, bdev, lrr, lrr_gc, child_vs_mid, father_vs_mother, index, bin_size,
                autosomal, par, pmed, prsd, masked, pbdev, masked_bands, bdev_adj, het_rel)


def gc_correct(lrr, gc, fit_mask, window=41):
    """A running median of the LRR against GC over the autosomal bins (the window in bins, over GC-sorted bins), subtracted from
    every bin by its GC; the result re-centred so the autosomal median is zero. Robust to the events themselves (a median)."""
    out = lrr.copy()
    ok = fit_mask & np.isfinite(lrr) & np.isfinite(gc)
    if ok.sum() < 3 * window:
        return out
    g, y = gc[ok], lrr[ok]
    o = np.argsort(g)
    g, y = g[o], y[o]
    h = window // 2
    curve = np.array([np.median(y[max(0, i - h):i + h + 1]) for i in range(len(y))])
    # one value per distinct GC (the interpolation needs increasing x)
    gu, inv = np.unique(g, return_inverse=True)
    cu = np.array([np.median(curve[inv == k]) for k in range(len(gu))])
    have = np.isfinite(lrr) & np.isfinite(gc)
    out[have] = lrr[have] - np.interp(gc[have], gu, cu)
    med = np.median(out[fit_mask & np.isfinite(out)]) if (fit_mask & np.isfinite(out)).any() else 0.0
    out[np.isfinite(out)] -= med
    return out


def chrom_level(bins, m, chrom, min_sites=20, min_bins=5):
    """A chromosome's copy number from the member's corrected LRR over its usable bins (outside the pseudoautosomal regions,
    with at least min_sites sites): (raw copies 2 x 2^median LRR, rounded copies, usable bins); (NaN, NaN, bins) with fewer
    than min_bins usable bins."""
    sl = bins.of(chrom)
    if sl.start == sl.stop:
        return NA, NA, 0
    y = bins.lrr_gc[m][sl]
    ok = np.isfinite(y) & (bins.n_dp[m][sl] >= min_sites) & ~bins.par[sl]
    if ok.sum() < min_bins:
        return NA, NA, int(ok.sum())
    raw = 2.0 * 2.0 ** float(np.median(y[ok]))
    return raw, int(round(raw)), int(ok.sum())


def d_grid():
    return np.arange(0.0, 0.4751, 0.005)


def d_hat(alt, dp):
    """The band deviation d at heterozygous sites: the maximum-likelihood d of alt ~ Binomial(dp, 1/2 +- d), each side with
    probability 1/2; with the log-likelihood ratio against d = 0 (the evidence that the bands are split at all). Each site
    is weighted by its depth through the binomial, so the estimate does not inflate at low depth as the mean of |BAF - 1/2| does."""
    alt = np.asarray(alt, float)
    dp = np.asarray(dp, float)
    ok = dp > 0
    alt, dp = alt[ok], dp[ok]
    if len(dp) < 5:
        return NA, NA
    grid = d_grid()
    k, n = alt[:, None], dp[:, None]
    p1 = 0.5 + grid[None, :]
    p0 = 0.5 - grid[None, :]
    with np.errstate(divide="ignore"):
        l1 = k * np.log(p1) + (n - k) * np.log(np.maximum(p0, 1e-12))
        l0 = k * np.log(np.maximum(p0, 1e-12)) + (n - k) * np.log(p1)
    m = np.maximum(l1, l0)
    ll = (m + np.log(np.exp(l1 - m) + np.exp(l0 - m)) - np.log(2.0)).sum(axis=0)
    i = int(np.argmax(ll))
    return float(grid[i]), float(ll[i] - ll[0])


def f_from_lrr(lrr, kind):
    """The cell fraction from the depth: gain 2^lrr = 1 + f/2, loss 2^lrr = 1 - f/2."""
    r = 2.0 ** lrr
    return float(np.clip(2 * (r - 1), 0, 2)) if kind == "gain" else float(np.clip(2 * (1 - r), 0, 1))


def f_from_d(d, kind):
    """The cell fraction from the band deviation: gain d = f/(2(2+f)); loss d = f/(2(2-f)); LOH d = f/2."""
    if not np.isfinite(d):
        return NA
    if kind == "gain":
        return float(np.clip(4 * d / max(1 - 2 * d, 1e-9), 0, 2))
    if kind == "loss":
        return float(np.clip(4 * d / (1 + 2 * d), 0, 1))
    return float(np.clip(2 * d, 0, 1))


def site_stats(sites, m, lo, hi, min_dp, min_gq):
    """A member's called and heterozygous sites in [lo, hi], the band deviation at the heterozygous ones, the alt fractions."""
    sel = (sites.pos >= lo) & (sites.pos <= hi)
    dp, alt, gt, gq = sites.dp[m][sel], sites.alt[m][sel], sites.gt[m][sel], sites.gq[m][sel]
    called = (dp >= min_dp) & (gq >= min_gq) & (gt >= GT_HOMREF)
    het = called & (gt == GT_HET)
    d, llr = d_hat(alt[het], dp[het])
    return dict(n_called=int(called.sum()), n_het=int(het.sum()), het_rate=float(het.sum() / called.sum()) if called.sum() else NA, d_hat=d, llr_baf=llr)
