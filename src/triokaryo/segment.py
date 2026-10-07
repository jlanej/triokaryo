"""Segmentation of a member's tracks along each chromosome and the calls: gains and losses from the LRR, copy-neutral
loss of heterozygosity from the B-allele band deviation and the heterozygosity rate; each with the cell fraction
estimated from the depth and, independently, from the bands; classified as whole chromosome, arm or stretch."""
from dataclasses import dataclass, field, asdict

import numpy as np

from .model import NA, MEMBERS, f_from_d, f_from_lrr, site_stats

PARAMS = dict(min_abs=0.07,      # minimum |mean LRR| of a gain or loss: log2(1 + 0.10/2), a cell fraction of about 10%
              z=5.0,             # minimum split statistic for a segment boundary (binary segmentation)
              min_len=5,         # minimum bins per segment
              min_sites=20,      # minimum sites per bin for a depth value
              loh_d=0.04,        # minimum band deviation of a copy-neutral LOH (a cell fraction of about 8%) ...
              loh_llr=10.0,      # ... and its minimum log-likelihood ratio against d = 0
              het_rel=0.35,      # heterozygosity rate at or below this fraction of the member's own: a constitutional loss of heterozygosity
              bdev_excess=0.03,  # a segment's mean band deviation this far above the member's median is an LOH candidate
              min_f=0.10,        # minimum cell fraction reported
              span_frac=0.90)    # a segment covering this fraction of a chromosome (arm) is classified whole (p or q)


@dataclass
class Event:
    sample: str
    role: str
    chrom: str
    start: int
    end: int
    span: str
    type: str
    lrr: float = NA
    lrr_se: float = NA
    n_bins: int = 0
    f_lrr: float = NA
    d_hat: float = NA
    llr_baf: float = NA
    f_baf: float = NA
    het_rate: float = NA
    het_rate_rel: float = NA
    n_het: int = 0
    n_called: int = 0
    mie_rate: float = NA
    origin: str = ""
    origin_llr: float = NA
    origin_n: int = 0
    inheritance: str = ""
    external: str = ""
    note: str = ""
    source: str = "depth"       # depth (LRR), bands (the folded bands / the heterozygosity rate), phased (the phased scan)
    phase_shift: float = NA     # the phased fraction's shift from one half over the event (child: maternal allele; parent: transmitted allele)
    phase_se: float = NA
    n_phased: int = 0
    f_phase: float = NA         # the cell fraction from the phased shift
    origin_phase: str = ""      # what the sign says
    start_fine: float = NA      # the edges at site resolution, from the phased sites
    end_fine: float = NA
    edge_sites: int = 0
    homologues: str = ""        # a child's gain, LOH or heterodisomy: the parent's two copies one homologue or two (the auxiliary track)
    hetero_share: float = NA    # the share of the event's windows where they differ

    @property
    def f(self):
        """The cell fraction: from the depth for a depth-called gain or loss, else from the bands."""
        return self.f_lrr if np.isfinite(self.f_lrr) and self.type in ("gain", "loss") else self.f_baf

    def overlap(self, other):
        """The intersection as a fraction of the shorter segment (lenient: a caller's fragment inside this event scores 1)."""
        a = max(self.start, other.start)
        b = min(self.end, other.end)
        if b <= a:
            return 0.0
        return (b - a) / min(self.end - self.start, other.end - other.start)

    def reciprocal_overlap(self, other):
        """The intersection as a fraction of the longer segment: at least t means each segment covers at least t of the other."""
        a = max(self.start, other.start)
        b = min(self.end, other.end)
        if b <= a:
            return 0.0
        return (b - a) / max(self.end - self.start, other.end - other.start)

    def as_dict(self):
        d = asdict(self)
        d["f"] = self.f
        return d


def robust_sd(y):
    d = np.diff(y[np.isfinite(y)])
    if len(d) < 4:
        return NA
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2))


def binary_segmentation(y, min_len=5, z=5.0, sd=None, max_segments=20):
    """Boundaries in y (no NaN) by recursive binary splitting on the two-sample z of the means, with the noise sd given
    (robust, from the first differences) rather than estimated inside each piece. Returns [(i0, i1)] covering y."""
    n = len(y)
    if sd is None or not np.isfinite(sd) or sd <= 0:
        sd = robust_sd(y)
    if n < 2 * min_len or not np.isfinite(sd) or sd <= 0:
        return [(0, n)]
    out = []
    work = [(0, n)]
    while work and len(out) + len(work) < max_segments:
        a, b = work.pop()
        seg = y[a:b]
        m = b - a
        if m < 2 * min_len:
            out.append((a, b))
            continue
        cs = np.cumsum(seg)
        tot = cs[-1]
        k = np.arange(min_len, m - min_len + 1)
        left = cs[k - 1] / k
        right = (tot - cs[k - 1]) / (m - k)
        zs = np.abs(left - right) / (sd * np.sqrt(1.0 / k + 1.0 / (m - k)))
        j = int(np.argmax(zs))
        if zs[j] >= z:
            cut = a + int(k[j])
            work.append((cut, b))
            work.append((a, cut))
        else:
            out.append((a, b))
    out += work
    return sorted(out)


def refine_boundaries(y, segs, max_shift=4):
    """Each boundary moved within +-max_shift bins to the position minimising the two segments' squared deviations from their
    means: the recursive split lands a bin or two off an edge where a segment is short, and a stray bin left outside would
    otherwise form a weak event of its own."""
    segs = [list(s) for s in segs]
    for i in range(len(segs) - 1):
        a0, b1 = segs[i][0], segs[i + 1][1]
        best, best_cost = segs[i][1], None
        for b in range(max(a0 + 1, segs[i][1] - max_shift), min(b1 - 1, segs[i][1] + max_shift) + 1):
            left, right = y[a0:b], y[b:b1]
            cost = float(((left - left.mean()) ** 2).sum() + ((right - right.mean()) ** 2).sum())
            if best_cost is None or cost < best_cost - 1e-12:
                best, best_cost = b, cost
        segs[i][1] = segs[i + 1][0] = best
    return [tuple(s) for s in segs]


def merge_similar(y, segs, sd, z=3.0):
    """Adjacent segments whose means differ by less than z standard errors are one segment."""
    segs = list(segs)
    changed = True
    while changed and len(segs) > 1:
        changed = False
        for i in range(len(segs) - 1):
            a0, a1 = segs[i]
            b0, b1 = segs[i + 1]
            ma, mb = np.mean(y[a0:a1]), np.mean(y[b0:b1])
            se = sd * np.sqrt(1.0 / (a1 - a0) + 1.0 / (b1 - b0))
            if abs(ma - mb) < z * se:
                segs[i:i + 2] = [(a0, b1)]
                changed = True
                break
    return segs


def span_of(chrom, start, end, valid_starts, genome, frac):
    """whole / p / q / stretch, by the share of the chromosome's (arm's) usable bins the segment covers."""
    L = genome.length[chrom]
    pe = genome.p_end.get(chrom, 0)
    vs = np.asarray(valid_starts)
    inside = (vs >= start) & (vs < end)
    if inside.sum() >= frac * len(vs) and len(vs):
        return "whole"
    for arm, lo, hi in (("p", 0, pe), ("q", pe, L)):
        arm_bins = (vs >= lo) & (vs < hi)
        if arm_bins.sum() and inside[arm_bins].sum() >= frac * arm_bins.sum() and inside[~arm_bins].sum() < 0.5 * max(1, (~arm_bins).sum()):
            return arm
    return "stretch"


def join_pieces(events, bins, genome, max_gap_bins=3, max_dlrr=0.08, frac=0.90):
    """Adjacent events of one member, chromosome and type, separated by at most max_gap_bins bins (e.g. a masked centromere)
    and differing in level by at most max_dlrr, are joined into one event: the LRR re-averaged over the pieces, the span
    re-classified. Site counts are summed; d_hat and the bands' cell fraction are taken from the larger piece."""
    out = []
    by = {}
    for e in events:
        by.setdefault((e.sample, e.chrom, e.type), []).append(e)
    for key, evs in by.items():
        evs.sort(key=lambda e: e.start)
        merged = [evs[0]]
        for e in evs[1:]:
            last = merged[-1]
            gap = (e.start - last.end) / bins.bin_size
            if gap <= max_gap_bins and (abs(e.lrr - last.lrr) <= max_dlrr if e.type != "LOH" else abs(e.f_baf - last.f_baf) <= 0.15):
                n1, n2 = last.n_bins, e.n_bins
                last.lrr = (last.lrr * n1 + e.lrr * n2) / (n1 + n2)
                last.lrr_se = last.lrr_se * np.sqrt(n1) / np.sqrt(n1 + n2)
                last.n_bins = n1 + n2
                last.end = e.end
                last.n_het += e.n_het
                last.n_called += e.n_called
                if e.n_bins > n1:
                    last.d_hat, last.llr_baf, last.f_baf = e.d_hat, e.llr_baf, e.f_baf
                if last.type in ("gain", "loss"):
                    last.f_lrr = f_from_lrr(last.lrr, last.type)
                last.note = (last.note + "; " if last.note else "") + "joined across %d bin(s) the panel left out" % int(round(gap)) if gap > 0 else last.note
            else:
                merged.append(e)
        out += merged
    starts_by_chrom = {}
    for e in out:
        sl = bins.of(e.chrom)
        valid = np.isfinite(bins.lrr_gc[0][sl]) | np.isfinite(bins.lrr_gc[1][sl]) | np.isfinite(bins.lrr_gc[2][sl])
        starts_by_chrom.setdefault(e.chrom, bins.start[sl][valid])
        e.span = span_of(e.chrom, e.start, e.end, starts_by_chrom[e.chrom], genome, frac)
    return sorted(out, key=lambda e: (("child", "father", "mother").index(e.role), genome.chroms.index(e.chrom), e.start))


def call_member(bins, scan, m, sample, genome, params=None, sex=""):
    """Every event of one member, chromosome by chromosome. Returns (events, x_copies) - x_copies the copies of X the member's
    X depth reads (1 or 2; the X events are read against that baseline, so a 47,XXY's X is '2 copies', not a gain)."""
    P = dict(PARAMS, **(params or {}))
    role = MEMBERS[m]
    events = []
    # the member's own baselines: the band deviation and the heterozygosity rate over the autosomes. A bin carries a
    # heterozygosity rate from min_called confident sites: 50, or half the member's typical bin where the VCF is sparser (a
    # thinned scan, a small mock)
    auto = bins.autosomal
    typical = float(np.median(bins.n_called[m][auto & (bins.n_called[m] > 0)])) if (auto & (bins.n_called[m] > 0)).any() else 0.0
    min_called = int(min(50, max(10, typical // 2)))
    base_bdev = float(np.nanmedian(bins.bdev[m][auto])) if np.isfinite(bins.bdev[m][auto]).any() else NA
    base_het = float(np.nanmedian(bins.het_rate[m][auto & (bins.n_called[m] >= min_called)])) if (auto & (bins.n_called[m] >= min_called)).any() else NA
    x_copies = NA
    for chrom in bins.index:
        sl = bins.of(chrom)
        y_all = bins.lrr_gc[m][sl].copy()
        valid = np.isfinite(y_all) & (bins.n_dp[m][sl] >= P["min_sites"])
        if chrom == "chrY":
            continue                                           # too few sites in a VCF to read; the depth tool's job
        if chrom == "chrX":
            v = valid & ~bins.par[sl]
            if v.sum() >= P["min_len"]:
                med = float(np.median(y_all[v]))
                x_copies = int(round(2 * 2 ** med))
                y_all = y_all - med                            # the X read against the member's own X copy state
            valid = v
        if valid.sum() < P["min_len"]:
            continue
        y = y_all[valid]
        starts = bins.start[sl][valid]
        ends = bins.end[sl][valid]
        sd = robust_sd(y)
        if not np.isfinite(sd) or sd <= 0:
            continue
        segs = merge_similar(y, refine_boundaries(y, binary_segmentation(y, P["min_len"], P["z"], sd)), sd)
        sites = scan.sites(chrom)
        covered = np.zeros(len(y), dtype=bool)
        for a, b in segs:
            if b - a < P["min_len"]:
                continue                                       # a remnant shorter than a segment: not an event
            mean = float(np.mean(y[a:b]))
            se = sd / np.sqrt(b - a)
            if abs(mean) < max(P["min_abs"], 3 * se):
                continue
            kind = "gain" if mean > 0 else "loss"
            f_l = f_from_lrr(mean, kind)
            if f_l < P["min_f"]:
                continue
            lo, hi = int(starts[a]), int(ends[b - 1])
            st = site_stats(sites, m, lo + 1, hi, P["min_dp"], P["min_gq"]) if sites is not None else {}
            ev = Event(sample, role, chrom, lo, hi, span_of(chrom, lo, hi, starts, genome, P["span_frac"]), kind, lrr=mean, lrr_se=float(se),
                       n_bins=int(b - a), f_lrr=f_l, d_hat=st.get("d_hat", NA), llr_baf=st.get("llr_baf", NA), f_baf=f_from_d(st.get("d_hat", NA), kind),
                       het_rate=st.get("het_rate", NA), n_het=st.get("n_het", 0), n_called=st.get("n_called", 0))
            ev.het_rate_rel = ev.het_rate / base_het if np.isfinite(base_het) and base_het > 0 and np.isfinite(ev.het_rate) else NA
            if chrom == "chrX" and x_copies == 1 and kind == "loss":
                ev.note = "a loss on a single X: a mosaic loss of the one X"
            events.append(ev)
            covered[a:b] = True
        # copy-neutral loss of heterozygosity: bands split (bdev above the member's own) or no heterozygous calls,
        # over bins the depth did not call
        if chrom == "chrX" and x_copies == 1:
            continue                                           # a single X has no heterozygous sites to read
        banned = bins.masked_bands[sl][valid] if bins.masked_bands is not None else np.zeros(int(valid.sum()), dtype=bool)
        use_panel = bins.het_rel is not None and np.isfinite(bins.het_rel[m][auto]).sum() > 100
        if use_panel:                                              # the member's own baselines on the panel-relative tracks
            base_bdev_t = float(np.nanmedian(bins.bdev_adj[m][auto])) if np.isfinite(bins.bdev_adj[m][auto]).any() else NA
            base_het_t = float(np.nanmedian(bins.het_rel[m][auto & (bins.n_called[m] >= min_called)])) if (auto & (bins.n_called[m] >= min_called)).any() else NA
        else:
            base_bdev_t, base_het_t = base_bdev, base_het
        for track, kind_test in (("bdev", "bands"), ("het_rate", "rate")):
            src = (bins.bdev_adj if use_panel else bins.bdev) if track == "bdev" else (bins.het_rel if use_panel else bins.het_rate)
            t_all = src[m][sl][valid].copy()
            if track == "bdev":
                ok = np.isfinite(t_all) & (bins.n_het[m][sl][valid] >= 5) & ~banned
                if not np.isfinite(base_bdev_t):
                    continue
                t = np.where(ok, t_all, base_bdev_t)
                base_t = base_bdev_t
            else:
                ok = np.isfinite(t_all) & (bins.n_called[m][sl][valid] >= min_called)
                if not np.isfinite(base_het_t):
                    continue
                t = np.where(ok, t_all, base_het_t)
                base_t = base_het_t
            if ok.sum() < P["min_len"]:
                continue
            sdt = robust_sd(t)
            if not np.isfinite(sdt) or sdt <= 0:
                sdt = max(float(np.nanstd(t)), 1e-3)
            tsegs = merge_similar(t, refine_boundaries(t, binary_segmentation(t, P["min_len"], P["z"], sdt)), sdt)
            for a, b in tsegs:
                if b - a < P["min_len"] or covered[a:b].mean() > 0.5 or banned[a:b].mean() > 0.5:
                    continue
                if track == "bdev":
                    if np.mean(t[a:b]) - base_t < P["bdev_excess"]:
                        continue
                else:
                    if not (np.isfinite(base_t) and base_t > 0 and np.mean(t[a:b]) <= P["het_rel"] * base_t):
                        continue
                lo, hi = int(starts[a]), int(ends[b - 1])
                st = site_stats(sites, m, lo + 1, hi, P["min_dp"], P["min_gq"]) if sites is not None else {}
                rel = float(np.mean(t[a:b]) / base_t) if track == "het_rate" and np.isfinite(base_t) and base_t > 0 else (
                    st.get("het_rate", NA) / base_het if np.isfinite(base_het) and base_het > 0 and np.isfinite(st.get("het_rate", NA)) else NA)
                d = st.get("d_hat", NA)
                if use_panel and np.isfinite(d):                   # the bands' deviation net of what the panel's genomes show there
                    excess = float(np.nanmean(bins.panel_bdev[sl][valid][a:b])) - float(np.nanmedian(bins.panel_bdev)) if np.isfinite(bins.panel_bdev[sl][valid][a:b]).any() else 0.0
                    d = max(0.0, d - max(0.0, excess))
                banded = np.isfinite(d) and d >= P["loh_d"] and st.get("llr_baf", 0) >= P["loh_llr"]
                silent = np.isfinite(rel) and rel <= P["het_rel"] and st.get("n_called", 0) >= max(40, 2 * min_called)
                if not (banded or silent):
                    continue
                f_b = 1.0 if silent else f_from_d(d, "LOH")            # no heterozygous calls: in every cell, whatever the few bands left say
                if f_b < P["min_f"]:
                    continue
                ev = Event(sample, role, chrom, lo, hi, span_of(chrom, lo, hi, starts, genome, P["span_frac"]), "LOH", lrr=float(np.mean(y[a:b])),
                           lrr_se=float(sd / np.sqrt(b - a)), n_bins=int(b - a), d_hat=d, llr_baf=st.get("llr_baf", NA), f_baf=f_b,
                           het_rate=st.get("het_rate", NA), het_rate_rel=rel, n_het=st.get("n_het", 0), n_called=st.get("n_called", 0), source="bands",
                           note="no heterozygous calls: in every cell (a uniparental disomy, or a deletion the depth missed)" if silent else "")
                events.append(ev)
                covered[a:b] = True
    return join_pieces(events, bins, genome, frac=P["span_frac"]), x_copies
