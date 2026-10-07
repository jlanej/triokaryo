"""Transmission phasing of the B-allele fraction.

The alt-allele fraction at a heterozygous site is 1/2 +- d with the sign of d unknown, so BAF is conventionally read
folded (|BAF - 1/2|), an estimator biased upward by binomial noise (0.07 at 30x when d = 0). The trio supplies the
sign. Where the parents are opposite homozygotes, the parental origin of each of the child's alleles is known
regardless of the child's call, and the fraction of reads carrying the maternal allele (the maternal fraction) is 1/2
with equal homologue copy, 1/2 + d or 1/2 - d along an event according to the parent of origin, and 1 where no paternal
copy is present (deletion, isodisomy, heterodisomy). At a parent's heterozygous site where the child is homozygous,
the transmitted allele is known, and the fraction of reads carrying it shows whether an event in the parent lies on
the transmitted homologue. These are the main tracks: the pooled fraction over any set of sites is an unbiased
estimate of d with binomial error, and a phasing error (a genotype error elsewhere in the trio) flips one site's
contribution and so attenuates a signal rather than creating one. The reference bias (alt fraction slightly under
1/2 at heterozygous sites) is removed by the member's genome-wide median, applied with the sign of each site's tag.

Sites where only one parent is homozygous also phase the child (the homozygous parent's allele is known, the other
allele is the other parent's), but only while each parent contributed one homologue. Where the child carries two
different homologues of one parent (meiotic trisomy, heterodisomy), the sites tagged by the other parent's
homozygous genotype are misassigned. They are therefore pooled as auxiliary tracks: along a maternal meiotic trisomy
the track read at the father's homozygous sites departs from the main track (1/3 against 2/3) where the two maternal
copies differ and rejoins it where a crossover has made them identical, mapping the meiotic error and its crossovers;
a parent's auxiliary track (child heterozygous, other parent homozygous) does the same from the parent's side.

From the main tracks: pooled windows (a few hundred kb) and their step fit (total-variation denoising); per event the
phased shift, cell fraction, parent of origin and homologue count; boundaries at site resolution; per-homologue copy
number (the LRR-derived copy number split by the fraction: maternal and paternal along the child, transmitted and
untransmitted along a parent); and a scan of the track for events below the depth detection limit (a gain or loss in
a few per cent of cells, typed by the depth where it deviates) and for uniparental heterodisomy, which depth, folded
BAF and heterozygosity rate all miss (the child homozygous for one parent's allele at every informative site yet
heterozygous wherever that parent is)."""
from dataclasses import dataclass, field

import numpy as np

from .model import NA, MEMBERS, f_from_d, site_stats
from .segment import Event, binary_segmentation, merge_similar, refine_boundaries, robust_sd, span_of
from .smooth import tv_denoise_nan, tv_lambda, windows_by_count
from .trio import _confident
from .vcfscan import GT_HET, GT_HOMALT, GT_HOMREF

POOL_TARGET_BP = 400_000     # a window of pooled sites spans about this much ...
POOL_MIN_SITES = 20          # ... never fewer sites than this ...
POOL_MAX_SITES = 200         # ... nor more
TV_K = 2.5                   # the step fit's penalty, in units of the track's noise sd
EDGE_BINS = 1.5              # an event's edge is sought within this many bins of the bin edge (wider where the sites are sparse) ...
EDGE_MIN_SITES = 10          # ... with at least this many phased sites each side of the bin edge
EDGE_MIN_SIDE = 3            # sites each side of the split itself
ORIGIN_MIN_SITES = 40        # phased sites a parent-of-origin reading needs ...
ORIGIN_MIN_Z = 3.0           # ... and the shift's size in standard errors
HOMOLOGUE_MIN_WINDOWS = 4    # windows an event needs for the one-or-two-homologues reading ...
HOMOLOGUE_MIN_SHIFT = 0.05   # ... and the shift it needs (a gain in a fifth of the cells): below it the auxiliary tracks' signs are noise
SHARED_SHIFT = 0.05          # a window deviating this far from 1/2 in two or more members at once (the other's at least half of this member's) is shared (paralogy): no member's event
DISAGREE_SITES = 500         # minimum phased sites before an event's phased cell fraction is compared with the depth's
STAGE = dict(cen_windows=8,          # windows nearest the centromere whose state classifies the meiotic stage ...
             cen_max_bp=15_000_000,  # ... within this distance of the centromere
             smooth=5,               # running majority over this many windows before states and crossovers are read
             state_frac=0.7)         # share of the centromeric windows that must agree for a centromeric state
SCAN = dict(min_len=8,       # windows per segment of the phased scan
            min_bp=2_000_000,  # a phased find spans at least this (a dense cluster of sites makes many windows of a few hundred kb)
            z=5.0,           # the split statistic (binary segmentation, as the depth's)
            min_shift=0.015, # the smallest shift reported: a gain or loss in about 6% of cells ...
            min_shift_loh=0.025,  # ... and, where the depth is flat and the bands are the only evidence, a copy-neutral LOH in 5%
            min_z=5.0,       # the segment's shift in standard errors (the empirical, not the binomial, error)
            agree=0.75,      # minimum share of a segment's windows with the sign of its mean (a few deviating windows do not make an event)
            lrr_z=3.0,       # the depth's lean, in standard errors, that types a phased find as a gain or loss
            upd_shift=0.4,   # a flat-depth shift this large (both copies from one parent) ...
            upd_het=0.5)     # ... with the heterozygosity rate at least this share of the member's own: a heterodisomy, not an LOH


@dataclass
class PhasedTrack:
    chrom: str
    idx: np.ndarray          # indices into the chromosome's sites: the main track's, in order
    frac: np.ndarray         # per site: the tagged allele's fraction, reference bias out
    tag_alt: np.ndarray      # per site: the tagged allele is the alt one
    w_start: np.ndarray      # windows: first site's position ...
    w_end: np.ndarray        # ... last site's
    w_mid: np.ndarray
    w_n: np.ndarray
    w_dp: np.ndarray
    w_frac: np.ndarray       # the pooled fraction (depth-weighted)
    w_se: np.ndarray         # its binomial error at one half
    w_bin: np.ndarray        # the bin of the window's middle
    w_step: np.ndarray       # the step fit
    w_lrr: np.ndarray        # the bin's LRR (the member's GC- and panel-corrected one)
    w_copies_tag: np.ndarray    # the tagged homologue's copies: the LRR's step fit's copies (2 * 2^LRR) times the fraction's step fit
    w_copies_other: np.ndarray
    aux: dict = field(default_factory=dict)   # name -> (per-window pooled fraction, per-window sites): the auxiliary tracks
    aux_sites: dict = field(default_factory=dict)   # name -> (idx, frac) per site
    w_shared: np.ndarray = None   # windows deviating in two or more members (set across the trio): excluded from the scan and the fits
    w_rejected: np.ndarray = None # windows inside a segment the scan rejected (paralogy, a dense cluster): no copy number drawn there


def phase_classes(sites, m, min_dp, min_gq, x_hemizygous_child=False):
    """The member's phased site classes over one chromosome: {name: (indices, tag_alt)}.
    The child: 'main' - the parents opposite homozygotes, the child any confident call; tag the mother's allele.
               'mother_hom' - the mother homozygous, the father heterozygous, the child heterozygous; tag the mother's allele.
               'father_hom' - the father homozygous, the mother heterozygous, the child heterozygous; tag the allele the father did not give.
    A parent: 'main' - the parent heterozygous, the child homozygous; tag the child's allele.
              'child_het' - the parent and the child heterozygous, the other parent homozygous; tag the allele the other parent did not give.
    x_hemizygous_child: the X of a son, whose one X is maternal, so a father's allele differing from the son's is not an error
    and the mother's main class keeps those sites."""
    gt, dp, gq = sites.gt, sites.dp, sites.gq
    conf = [_confident(dp[i], gq[i], gt[i], min_dp, min_gq) & (gt[i] >= GT_HOMREF) for i in range(3)]
    hom = [conf[i] & (gt[i] != GT_HET) for i in range(3)]
    het = [conf[i] & (gt[i] == GT_HET) for i in range(3)]
    conf[m] &= dp[m] >= min_dp                                     # the member's own call needs reads to give a fraction (not a merge's absent call)
    hom[m] &= dp[m] >= min_dp
    out = {}
    if m == 0:
        opp = hom[1] & hom[2] & (gt[1] != gt[2])
        main = opp & conf[0]
        out["main"] = (np.flatnonzero(main), (gt[2] == GT_HOMALT)[main])
        c2 = hom[2] & het[1] & het[0]
        out["mother_hom"] = (np.flatnonzero(c2), (gt[2] == GT_HOMALT)[c2])
        c3 = hom[1] & het[2] & het[0]
        out["father_hom"] = (np.flatnonzero(c3), (gt[1] == GT_HOMREF)[c3])
    else:
        o = 3 - m
        main = het[m] & hom[0]
        if not (x_hemizygous_child and m == 2):
            main &= ~(hom[o] & (gt[o] != gt[0]))                      # a homozygous child opposite to the other parent: not this parent's to phase
        out["main"] = (np.flatnonzero(main), (gt[0] == GT_HOMALT)[main])
        cb = het[m] & het[0] & hom[o]
        out["child_het"] = (np.flatnonzero(cb), (gt[o] == GT_HOMREF)[cb])
    return out


def ref_bias(scan, genome, m, min_dp, min_gq, step=7):
    """The member's reference bias: the median of (alt fraction - 1/2) over its confident autosomal heterozygous sites."""
    v = []
    for c, s in scan.chroms.items():
        if c not in genome.autosomes:
            continue
        het = (s.gt[m] == GT_HET) & (s.dp[m] >= min_dp) & (s.gq[m] >= min_gq)
        i = np.flatnonzero(het)[::step]
        if len(i):
            v.append(s.alt[m][i] / s.dp[m][i] - 0.5)
    if not v:
        return 0.0
    v = np.concatenate(v)
    return float(np.median(v)) if len(v) >= 200 else 0.0


def _window_size(n_sites, span_bp):
    dens = n_sites / max(span_bp, 1.0)
    return int(np.clip(round(dens * POOL_TARGET_BP), POOL_MIN_SITES, POOL_MAX_SITES))


def _frac(sites, m, idx, tag, b):
    dp = sites.dp[m][idx].astype(float)
    alt = sites.alt[m][idx].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        f = np.where(tag, alt / dp, 1 - alt / dp) - np.where(tag, b, -b)
    return f, dp


def running_median(y, k=5):
    """The median of the finite values in a window of k around each point (NaN where the point is NaN)."""
    y = np.asarray(y, dtype=float)
    out = np.full(len(y), NA)
    h = k // 2
    for i in range(len(y)):
        if not np.isfinite(y[i]):
            continue
        v = y[max(0, i - h):i + h + 1]
        v = v[np.isfinite(v)]
        out[i] = float(np.median(v)) if len(v) else NA
    return out


def lrr_step(bins, m, median_bins=1):
    """The step fit of the member's LRR, chromosome by chromosome, one penalty for the genome; over a running median of
    median_bins bins first where asked (the copies: a dip of a bin or two is not a large event)."""
    y = bins.lrr_gc[m]
    lam = tv_lambda(y[bins.autosomal], TV_K)
    out = np.full(len(y), NA)
    for c in bins.index:
        sl = bins.of(c)
        v = running_median(y[sl], median_bins) if median_bins > 1 else y[sl]
        out[sl] = tv_denoise_nan(v, lam)
    return out


def _pool_by_windows(pos_all, frac_all, dp_all, w_start, w_end):
    """The depth-weighted pooled fraction of another site set over the same windows, and the sites in each."""
    n = len(w_start)
    out, cnt = np.full(n, NA), np.zeros(n, dtype=int)
    if not len(pos_all):
        return out, cnt
    a = np.searchsorted(pos_all, w_start, side="left")
    b = np.searchsorted(pos_all, w_end, side="right")
    for i in range(n):
        if b[i] > a[i]:
            d = dp_all[a[i]:b[i]]
            if d.sum() > 0:
                out[i] = float((frac_all[a[i]:b[i]] * d).sum() / d.sum())
                cnt[i] = b[i] - a[i]
    return out, cnt


def phased_tracks(scan, bins, genome, m, min_dp, min_gq, sex="", child_x_copies=NA):
    """The member's phased track per chromosome (autosomes and the X), with the pooled windows, the step fit, the auxiliary
    tracks and the parental copies. child_x_copies: the child's measured X copy number (1: a son, whose X is maternal).
    Returns ({chrom: PhasedTrack}, info)."""
    b = ref_bias(scan, genome, m, min_dp, min_gq)
    per = {}
    n_total, span = 0, 0
    hemi = np.isfinite(child_x_copies) and int(child_x_copies) == 1
    for c, s in scan.chroms.items():
        if c == "chrY" or c not in bins.index:
            continue
        cls = phase_classes(s, m, min_dp, min_gq, x_hemizygous_child=(c == "chrX" and hemi))
        keep = np.ones(s.n, dtype=bool)
        if c in genome.par or c == "chrX":
            keep &= ~s.par
        if bins.masked_bands is not None:                                   # bands split in every genome: paralogy, not phase
            a0 = bins.index[c][0]
            bi = a0 + np.minimum(s.pos // bins.bin_size, bins.index[c][1] - a0 - 1)
            keep &= ~bins.masked_bands[bi]
        cls = {k: (idx[keep[idx]], tag[keep[idx]]) for k, (idx, tag) in cls.items()}
        per[c] = cls
        if c in genome.autosomes:
            n_total += len(cls["main"][0])
            span += genome.length[c]
    w = _window_size(n_total, span)
    tracks = {}
    allw = []
    for c, cls in per.items():
        s = scan.sites(c)
        idx, tag = cls["main"]
        if len(idx) < POOL_MIN_SITES:
            continue
        frac, dp = _frac(s, m, idx, tag, b)
        pos = s.pos[idx]
        win = windows_by_count(pos, w)
        if not win:
            continue
        a0, a1 = bins.index[c]
        ws = np.array([pos[i0] for i0, i1 in win]); we = np.array([pos[i1 - 1] for i0, i1 in win])
        wm = np.array([int(np.median(pos[i0:i1])) for i0, i1 in win])
        wn = np.array([i1 - i0 for i0, i1 in win])
        wdp = np.array([dp[i0:i1].sum() for i0, i1 in win])
        wf = np.array([float((frac[i0:i1] * dp[i0:i1]).sum() / dp[i0:i1].sum()) for i0, i1 in win])
        wse = np.sqrt(0.25 / np.maximum(wdp, 1.0))
        wbin = a0 + np.minimum(wm // bins.bin_size, a1 - a0 - 1)
        lrr = bins.lrr_gc[m][wbin]
        t = PhasedTrack(c, idx, frac, tag, ws, we, wm, wn, wdp, wf, wse, wbin, np.full(len(win), NA), lrr, np.full(len(win), NA), np.full(len(win), NA))
        for k, (ai, at) in cls.items():
            if k == "main":
                continue
            af, adp = _frac(s, m, ai, at, b)
            t.aux[k] = _pool_by_windows(s.pos[ai], af, adp, ws, we)
            t.aux_sites[k] = (ai, af)
        tracks[c] = t
        allw.append(wf)
    fit_tracks(tracks, bins, m)
    return tracks, dict(ref_bias=b, window_sites=w, n_phased=int(sum(len(t.idx) for t in tracks.values())),
                        n_aux=int(sum(len(v[0]) for t in tracks.values() for v in t.aux_sites.values())))


def fit_tracks(tracks, bins, m):
    """The step fits of the member's windows (one penalty per genome; shared windows excluded) and the per-homologue copy
    number from them: the LRR step fit's copy number split by the fraction's step fit, each over a running median of five
    points first (a dip of a bin or two, or a deviating window or two, is not a large event); no copy number where the
    scan rejected a segment."""
    allw = [t.w_frac[np.isfinite(t.w_frac) & ~(t.w_shared if t.w_shared is not None else np.zeros(len(t.w_frac), dtype=bool))] for t in tracks.values()]
    allw = [v for v in allw if len(v)]
    lam = tv_lambda(np.concatenate(allw), TV_K) if allw else 0.0
    lstep = lrr_step(bins, m, median_bins=5)
    for t in tracks.values():
        y = t.w_frac.copy()
        if t.w_shared is not None:
            y[t.w_shared] = NA
        t.w_step = tv_denoise_nan(y, lam)
        smooth = tv_denoise_nan(running_median(y, 5), lam)
        if t.w_rejected is not None:
            smooth[t.w_rejected] = NA
        copies = 2.0 * 2.0 ** lstep[t.w_bin]
        t.w_copies_tag = copies * smooth
        t.w_copies_other = copies * (1 - smooth)


def mask_rejected(tracks_by_member, bins, rejected):
    """Windows inside segments the scan rejected carry no copy number; the fits are redone."""
    for m, tracks in enumerate(tracks_by_member):
        for c, t in tracks.items():
            rej = np.zeros(len(t.w_mid), dtype=bool)
            for r in rejected:
                if r["role"] == MEMBERS[m] and r["chrom"] == c:
                    rej |= (t.w_mid >= r["start"]) & (t.w_mid <= r["end"])
            t.w_rejected = rej
        fit_tracks(tracks, bins, m)


def mask_shared(tracks_by_member, bins, x_copies=None):
    """Windows deviating from 1/2 by SHARED_SHIFT or more in two or more members at once indicate paralogous sequence or an
    imbalance shared by the family, not a member's event: marked as shared on each member's track, and the fits redone.
    x_copies: {role: X copies}; a member with one X has a maternal fraction of 1 along the X by construction, which does
    not count against the others' X windows."""
    hemi = {i for i, role in enumerate(MEMBERS) if x_copies and np.isfinite(x_copies.get(role, NA)) and int(x_copies[role]) == 1}
    for m, tracks in enumerate(tracks_by_member):
        for c, t in tracks.items():
            mine = np.abs(t.w_frac - 0.5)
            n_parted = (mine >= SHARED_SHIFT).astype(int)
            for o, others in enumerate(tracks_by_member):
                if o == m or c not in others or (c == "chrX" and o in hemi):
                    continue
                u = others[c]
                j = np.searchsorted(u.w_start, t.w_mid, side="right") - 1          # the other member's window holding this one's middle
                ok = (j >= 0) & (j < len(u.w_mid))
                jj = np.clip(j, 0, max(len(u.w_mid) - 1, 0))
                inside = ok & (t.w_mid <= u.w_end[jj])
                parted = np.zeros(len(t.w_mid), dtype=int)
                theirs = np.abs(u.w_frac[jj[inside]] - 0.5)
                parted[inside] = ((theirs >= SHARED_SHIFT) & (theirs >= 0.5 * mine[inside])).astype(int)
                n_parted += parted
            t.w_shared = n_parted >= 2
    for m, tracks in enumerate(tracks_by_member):
        fit_tracks(tracks, bins, m)


def phased_reading(ev, track, sites, m):
    """The event's phased statistics: the pooled fraction's shift from 1/2, its binomial error, and the site count."""
    if track is None:
        return NA, NA, 0
    pos = sites.pos[track.idx]
    sel = (pos > ev.start) & (pos <= ev.end)
    n = int(sel.sum())
    if n == 0:
        return NA, NA, 0
    dp = sites.dp[m][track.idx[sel]].astype(float)
    shift = float((track.frac[sel] * dp).sum() / dp.sum()) - 0.5
    return shift, float(np.sqrt(0.25 / dp.sum())), n


def origin_from_shift(role, kind, shift):
    """What the sign of the phased shift says."""
    up = shift > 0
    if role == "child":
        if kind == "gain":
            return "extra copy maternal" if up else "extra copy paternal"
        if kind == "loss":
            return "paternal copy lost" if up else "maternal copy lost"
        if kind == "UPD":
            return "both copies maternal (heterodisomy)" if up else "both copies paternal (heterodisomy)"
        return "maternal copy retained (paternal replaced)" if up else "paternal copy retained (maternal replaced)"
    if kind == "gain":
        return "the duplicated homologue is the one passed to the child" if up else "the duplicated homologue is the one not passed to the child"
    if kind == "loss":
        return "the lost homologue is the one not passed to the child" if up else "the lost homologue is the one passed to the child"
    return "the retained homologue is the one passed to the child" if up else "the retained homologue is the one not passed to the child"


def _homologue_states(ev, track, absent_is_iso=False, upd_like=None):
    """Per window inside a child's event: True where the named parent's two copies in the child are different homologues, False
    where they are one, from the auxiliary track read at the other parent's homozygous sites (it departs from the main track,
    1/3 against 2/3 along a trisomy, where the copies differ; along a heterodisomy it sits at 1/2). Windows without that track
    are dropped, or counted as one homologue when absent_is_iso: a heterodisomy's isodisomic segments have no heterozygous
    child sites, so the track is absent there while the main track is present. Returns (window positions, states) or None."""
    if ev.role != "child" or ev.type not in ("gain", "LOH", "UPD") or not np.isfinite(ev.phase_shift) or track is None or abs(ev.phase_shift) < HOMOLOGUE_MIN_SHIFT:
        return None
    name = "father_hom" if ev.phase_shift > 0 else "mother_hom"          # a maternal event: the sites tagged by the father's homozygosity
    if name not in track.aux:
        return None
    af, an = track.aux[name]
    inside = (track.w_mid > ev.start) & (track.w_mid <= ev.end) & np.isfinite(track.w_frac)
    have = inside & np.isfinite(af) & (an >= 5)
    main = track.w_frac - 0.5
    other = np.where(np.isfinite(af), af, 0.5) - 0.5
    if (ev.type == "UPD") if upd_like is None else upd_like:
        # both copies from one parent in a cell fraction f: where they differ, the auxiliary track reads 1/(1+f) or f/(1+f),
        # i.e. 1/(1+f) - 1/2 from one half (0 at f = 1); where they are one homologue the child has no heterozygous site
        f = float(np.clip(ev.f, 0.0, 1.0)) if np.isfinite(ev.f) else 1.0
        differ = np.abs(other) < min(1.0 / (1.0 + f) - 0.5 + 0.12, 0.4)
    else:
        differ = np.sign(other) != np.sign(main)
    sel = inside if absent_is_iso else have
    if sel.sum() < HOMOLOGUE_MIN_WINDOWS:
        return None
    states = np.where(have, differ, False)[sel] if absent_is_iso else differ[sel]
    order = np.argsort(track.w_mid[sel])
    return track.w_mid[sel][order], states[order]


def homologues_reading(ev, track, upd_like=None):
    """For a child's gain, LOH or heterodisomy with a parent named: are the parent's two copies in the child one homologue or two?
    Returns (text, share of the event's windows in which the copies differ)."""
    st = _homologue_states(ev, track, absent_is_iso=bool(upd_like), upd_like=upd_like)
    if st is None:
        return "", NA
    share = float(st[1].mean())
    parent = "maternal" if ev.phase_shift > 0 else "paternal"
    if share >= 0.9:
        return "the two %s copies are different homologues throughout (meiotic)" % parent, share
    if share <= 0.1:
        return "the two %s copies are one homologue throughout (mitotic, or a meiosis II error without a crossover)" % parent, share
    return "the two %s copies differ over %.0f%% of the event and are one homologue over the rest (a meiotic error with crossovers)" % (parent, 100 * share), share


def running_majority(b, k=5):
    """The majority of a boolean sequence over a window of k around each point; a tie keeps the point's own value."""
    b = np.asarray(b, dtype=bool)
    out = b.copy()
    h = k // 2
    for i in range(len(b)):
        v = b[max(0, i - h):i + h + 1]
        s = int(v.sum())
        if 2 * s != len(v):
            out[i] = 2 * s > len(v)
    return out


def stage_reading(ev, track, genome, absent_is_iso=None, upd_like=None):
    """The meiotic stage of a child's whole-chromosome gain or heterodisomy, anchored at the centromere. The per-window states
    (two different homologues, or one) are smoothed by a running majority; the state of the windows nearest the centromere
    classifies the event: heterodisomic there, a meiosis I nondisjunction; isodisomic there with a heterodisomic segment
    elsewhere, meiosis II; isodisomic throughout, a mitotic duplication or a meiosis II error without a crossover. Each change
    of state along the chromosome is a crossover, placed midway between the two windows. Returns a dict (stage, centromere,
    n_crossovers, crossovers in Mb, n_windows) or None."""
    if ev.span != "whole":
        return None
    if absent_is_iso is None:
        absent_is_iso = (ev.type == "UPD") if upd_like is None else upd_like
    st = _homologue_states(ev, track, absent_is_iso=absent_is_iso, upd_like=upd_like)
    if st is None:
        return None
    mid, differ = st
    state = running_majority(differ, STAGE["smooth"])
    cen = genome.p_end.get(ev.chrom, 0)
    dist = np.abs(mid - cen)
    near = np.argsort(dist)[:STAGE["cen_windows"]]
    near = near[dist[near] <= STAGE["cen_max_bp"]]
    centro = ""
    if len(near) >= 3:
        h = float(state[near].mean())
        centro = "heterodisomic" if h >= STAGE["state_frac"] else "isodisomic" if h <= 1 - STAGE["state_frac"] else "mixed"
    change = np.flatnonzero(state[1:] != state[:-1])
    xo = [(mid[i] + mid[i + 1]) / 2.0 for i in change]
    het_any, iso_any = bool(state.any()), bool((~state).any())
    if centro == "heterodisomic":
        stage = "meiosis I"
    elif centro == "isodisomic":
        stage = "meiosis II" if het_any else "mitotic, or meiosis II without a crossover"
    elif not het_any:
        stage = "mitotic, or meiosis II without a crossover"                # isodisomic throughout, the centromere unresolved
    elif not iso_any:
        stage = "meiosis I"                                                 # heterodisomic throughout
    else:
        stage = ""                                                          # mixed, the centromere unresolved
    return dict(stage=stage, centromere=centro, n_crossovers=len(xo), crossovers=";".join("%.2f" % (x / 1e6) for x in xo), n_windows=int(len(mid)))


def refine_edges(ev, track, sites, bins, genome):
    """The event's start and end at site resolution: within EDGE_BINS bins of each bin edge, the split of the phased sites
    (signed by the event's shift) that best separates an outside at 1/2 from an inside away from it. (NA, NA, 0) where the
    event runs to the chromosome end, or the sites are too few."""
    if track is None or not np.isfinite(ev.phase_shift) or ev.phase_shift == 0:
        return NA, NA, 0
    pos = sites.pos[track.idx]
    y = np.sign(ev.phase_shift) * (track.frac - 0.5)
    L = genome.length[ev.chrom]
    half = EDGE_BINS * bins.bin_size
    out = [NA, NA]
    n_min = 0
    for k, edge in enumerate((ev.start, ev.end)):
        if (k == 0 and edge <= 0) or (k == 1 and edge >= L):
            continue
        h = half
        sel = np.flatnonzero((pos >= edge - h) & (pos <= edge + h))
        while min(int((pos[sel] < edge).sum()), int((pos[sel] >= edge).sum())) < EDGE_MIN_SITES and h < 6 * bins.bin_size:
            h *= 1.5                                                        # sparse sites: a wider look, until both sides of the bin edge have enough
            sel = np.flatnonzero((pos >= edge - h) & (pos <= edge + h))
        if min(int((pos[sel] < edge).sum()), int((pos[sel] >= edge).sum())) < EDGE_MIN_SITES:
            continue
        v = y[sel]
        n = len(v)
        cs = np.cumsum(v)
        tot = cs[-1]
        j = np.arange(EDGE_MIN_SIDE, n - EDGE_MIN_SIDE + 1)               # the split: left = v[:j], right = v[j:]
        ml = cs[j - 1] / j
        mr = (tot - cs[j - 1]) / (n - j)
        gain = j * ml ** 2 + (n - j) * mr ** 2                              # the squared-error reduction of a two-mean fit
        ok = (mr > ml) if k == 0 else (ml > mr)
        if not ok.any():
            continue
        gain = np.where(ok, gain, -np.inf)
        jj = int(j[int(np.argmax(gain))])
        out[k] = int((pos[sel[jj - 1]] + pos[sel[jj]]) // 2)
        n_min = n if not n_min else min(n_min, n)
    return out[0], out[1], n_min


def _hemizygous_x_reading(e, track):
    """The phased reading of a reported male's whole-X or segmental X event, whose baseline is one maternal X. A gain's extra
    X is maternal when the maternal fraction stays at 1 (shift about 1/2; the cell fraction is then not readable from the
    bands) and paternal when it falls to 1/(1+f) (shift (1-f)/(2(1+f)); f = (1-2s)/(1+2s)). A loss leaves the fraction at 1
    and has no parent to resolve. Sets f_phase, origin_phase, stage, homologues; returns True when handled."""
    if e.type == "gain":
        s = e.phase_shift
        if e.n_phased >= ORIGIN_MIN_SITES and (0.5 - s) >= ORIGIN_MIN_Z * e.phase_se and s < 0.4:
            e.origin_phase = "extra copy paternal"
            e.f_phase = float(np.clip((1 - 2 * s) / max(1 + 2 * s, 1e-9), 0, 2))
            if e.span == "whole":
                e.stage = "meiosis I (paternal: X and Y transmitted together)"
                e.centromere = ""
        elif e.n_phased >= ORIGIN_MIN_SITES and s >= 0.4:
            e.origin_phase = "extra copy maternal"
            e.f_phase = NA                                                  # the maternal fraction is 1 whatever the cell fraction
            e.homologues, e.hetero_share = homologues_reading(e, track, upd_like=True)
        return True
    if e.type == "loss":
        e.f_phase = NA
        e.origin_phase = "maternal copy lost (the single X)" if e.n_phased >= ORIGIN_MIN_SITES else ""
        return True
    return False


def annotate_events(events, tracks_by_member, scan, bins, genome, child_x_baseline=2):
    """Every event's phased statistics, parent of origin by the sign of the shift, homologue count, meiotic stage, and
    boundaries at site resolution; a note where the phased sign and the informative-site likelihood name different parents.
    An event annotated as a run of homozygosity receives no parent of origin or homologue count. child_x_baseline 1 (a
    reported male) reads the child's X events against one maternal X (_hemizygous_x_reading)."""
    for e in events:
        m = MEMBERS.index(e.role)
        track = tracks_by_member[m].get(e.chrom)
        sites = scan.sites(e.chrom)
        if track is None or sites is None:
            continue
        e.phase_shift, e.phase_se, e.n_phased = phased_reading(e, track, sites, m)
        if not np.isfinite(e.phase_shift):
            continue
        if e.role == "child" and e.chrom == "chrX" and child_x_baseline == 1:
            if _hemizygous_x_reading(e, track):
                if e.origin_phase == "extra copy maternal":
                    st = stage_reading(e, track, genome, upd_like=True)
                    if st:
                        e.stage, e.centromere, e.n_crossovers, e.crossovers = st["stage"], st["centromere"], st["n_crossovers"], st["crossovers"]
                e.start_fine, e.end_fine, e.edge_sites = refine_edges(e, track, sites, bins, genome) if e.span != "whole" else (NA, NA, 0)
                continue
        if e.role == "child" and e.chrom == "chrX" and e.type == "gain" and e.phase_shift < 0 and e.n_phased >= ORIGIN_MIN_SITES and -e.phase_shift >= ORIGIN_MIN_Z * e.phase_se:
            # a paternal extra X in a daughter (47,XXX, or a mosaic): the father has one X, so its two copies are one homologue
            e.f_phase = f_from_d(abs(e.phase_shift), e.type)
            e.origin_phase = "extra copy paternal"
            e.homologues, e.hetero_share = "the two paternal copies are one homologue (the father's single X)", 0.0
            if e.span == "whole":
                e.stage, e.centromere = "meiosis II, or post-zygotic (paternal: the father's single X duplicated)", "isodisomic"
            if e.origin and e.origin != e.origin_phase:
                e.note = (e.note + "; " if e.note else "") + "the phased bands and the opposite-homozygote sites name different parents"
            continue
        e.f_phase = f_from_d(abs(e.phase_shift), e.type)
        if e.source == "depth" and e.n_phased >= DISAGREE_SITES and np.isfinite(e.f_lrr) and np.isfinite(e.f_phase):
            if e.f_phase < 0.5 * e.f_lrr and e.f_lrr - e.f_phase > 0.05:
                e.note = (e.note + "; " if e.note else "") + ("the phased bands read a share of %.0f%% against the depth's %.0f%%: the depth's call may be an artefact "
                                                              "(GC, the panel's edge) rather than a copy change" % (100 * e.f_phase, 100 * e.f_lrr))
            elif e.f_phase > 1.5 * e.f_lrr and e.f_phase - e.f_lrr > 0.05:
                e.note = (e.note + "; " if e.note else "") + "the phased bands read a larger share (%.0f%%) than the depth (%.0f%%)" % (100 * e.f_phase, 100 * e.f_lrr)
        roh = "run of homozygosity" in e.note                               # a run of homozygosity has no parent of origin
        if e.n_phased >= ORIGIN_MIN_SITES and abs(e.phase_shift) >= ORIGIN_MIN_Z * e.phase_se and not roh:
            e.origin_phase = origin_from_shift(e.role, e.type, e.phase_shift)
            if e.role == "child" and e.origin and e.type != "UPD" and e.origin != e.origin_phase:
                e.note = (e.note + "; " if e.note else "") + "the phased bands and the opposite-homozygote sites name different parents"
            e.homologues, e.hetero_share = homologues_reading(e, track)
            st = stage_reading(e, track, genome)
            if st:
                e.stage, e.centromere, e.n_crossovers, e.crossovers = st["stage"], st["centromere"], st["n_crossovers"], st["crossovers"]
        e.start_fine, e.end_fine, e.edge_sites = refine_edges(e, track, sites, bins, genome)


def phased_scan(tracks, bins, scan, m, sample, genome, events, params, min_dp, min_gq, sex="", base_het=NA, rejected=None, x_copies=NA):
    """Events in the phased track not called by the depth: segments of the pooled windows shifted from 1/2, typed by the
    depth's deviation over the same bins (gain or loss); with flat depth, copy-neutral LOH at a cell fraction of 2d or,
    with a shift near 1/2 and heterozygosity retained, uniparental heterodisomy. Windows shared by two or more members
    are excluded first. A segment under min_bp, or whose windows are inconsistent with its mean, is rejected (appended
    to `rejected`, a list of dicts)."""
    S = dict(SCAN, **(params or {}))
    role = MEMBERS[m]
    out = []
    rejected = rejected if rejected is not None else []
    two_x = (np.isfinite(x_copies) and int(x_copies) == 2) if np.isfinite(x_copies) else (sex == "F")
    for chrom, t in tracks.items():
        if chrom not in genome.autosomes and not (chrom == "chrX" and two_x):
            continue                                                       # the X is scanned where the member has two copies
        mine = [e for e in events if e.sample == sample and e.chrom == chrom]
        free = np.ones(len(t.w_mid), dtype=bool)                            # the windows no event of the member's covers yet
        for e in mine:
            free &= (t.w_end < e.start) | (t.w_start > e.end)
        if t.w_shared is not None:
            free &= ~t.w_shared
        keep = np.flatnonzero(free & np.isfinite(t.w_frac))
        y = t.w_frac[keep] - 0.5
        if len(y) < max(4, S["min_len"] // 2):
            continue
        sd = robust_sd(y)
        if not np.isfinite(sd) or sd <= 0:
            continue
        if len(y) < 2 * S["min_len"]:
            segs = [(0, len(y))]                                            # too short to split: the chromosome as one segment
        else:
            segs = merge_similar(y, refine_boundaries(y, binary_segmentation(y, S["min_len"], S["z"], sd)), sd)
        sl = bins.of(chrom)
        lrr_all = bins.lrr_gc[m][sl]
        lsd = robust_sd(lrr_all[np.isfinite(lrr_all)])
        valid = np.isfinite(lrr_all)
        starts = bins.start[sl][valid]
        sites = scan.sites(chrom)
        for a, b in segs:
            ka, kb = keep[a], keep[b - 1]
            if kb - ka + 1 < S["min_len"] and b - a < len(y):                 # the span in windows, the left-out ones counted
                continue
            mean = float(np.mean(y[a:b]))
            se = sd / np.sqrt(b - a)
            if not np.isfinite(mean) or abs(mean) < max(S["min_shift"], S["min_z"] * se):
                continue
            lo, hi = int(t.w_start[ka]), int(t.w_end[kb])
            med = float(np.median(y[a:b]))
            agree = float(np.mean(np.sign(y[a:b]) == np.sign(mean)))
            why = ""
            if hi - lo < S["min_bp"]:
                why = "spans %.2f Mb, under the floor of %.0f Mb (a dense cluster of sites)" % ((hi - lo) / 1e6, S["min_bp"] / 1e6)
            elif abs(med) < S["min_shift"] or agree < S["agree"]:
                why = "its windows do not agree with its mean (median shift %+.3f, %.0f%% the mean's way): a few parted windows, not an event" % (med, 100 * agree)
            if why:
                rejected.append(dict(sample=sample, role=role, chrom=chrom, start=lo, end=hi, shift=mean, windows=int(b - a), reason=why))
                continue
            wb = np.unique(t.w_bin[keep[a:b]])
            v = bins.lrr_gc[m][wb]
            v = v[np.isfinite(v)]
            lz = float(np.mean(v) / (lsd / np.sqrt(len(v)))) if len(v) and np.isfinite(lsd) and lsd > 0 else 0.0
            kind = "gain" if lz >= S["lrr_z"] else "loss" if lz <= -S["lrr_z"] else "LOH"
            d = abs(mean)
            st = site_stats(sites, m, lo + 1, hi, min_dp, min_gq) if sites is not None else {}
            rel = st.get("het_rate", NA) / base_het if np.isfinite(base_het) and base_het > 0 and np.isfinite(st.get("het_rate", NA)) else NA
            if kind == "LOH" and d >= S["upd_shift"] and np.isfinite(rel) and rel >= S["upd_het"]:
                kind = "UPD"
            if kind == "LOH" and d < S["min_shift_loh"]:
                rejected.append(dict(sample=sample, role=role, chrom=chrom, start=lo, end=hi, shift=mean, windows=int(b - a),
                                     reason="the depth is flat and the shift (%.3f) is under the floor for the bands alone (%.3f)" % (d, S["min_shift_loh"])))
                continue
            if kind == "UPD":
                note = ("a uniparental heterodisomy: the depth is flat (%.1f standard errors), the heterozygosity kept (%.2f of the member's own), and the "
                        "child homozygous for one parent's allele wherever the parents are opposite homozygotes" % (lz, rel))
            elif kind == "LOH":
                note = ("from the phased bands: the depth is flat (%.1f standard errors); a copy-neutral LOH at this share, or a gain in %.0f%% / a loss in %.0f%% of cells" %
                        (lz, 100 * f_from_d(d, "gain"), 100 * f_from_d(d, "loss")))
            else:
                note = "from the phased bands: the depth leans the same way (%.1f standard errors) but under its own threshold" % lz
            ev = Event(sample, role, chrom, lo, hi, span_of(chrom, lo, hi, starts, genome, 0.90), kind, lrr=float(np.mean(v)) if len(v) else NA,
                       lrr_se=float(lsd / np.sqrt(len(v))) if len(v) and np.isfinite(lsd) else NA, n_bins=int(len(v)), d_hat=d, llr_baf=NA,
                       f_baf=f_from_d(d, kind), het_rate=st.get("het_rate", NA), het_rate_rel=rel, n_het=st.get("n_het", 0), n_called=st.get("n_called", 0),
                       source="phased", note=note)
            out.append(ev)
    return out
