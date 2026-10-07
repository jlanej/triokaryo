"""Trio analysis of each event: inheritance (an event of the same type in a parent), the parent of origin of the child's
events from the informative sites, and the Mendelian-error rate within the event.

At a site where the father is 0/0 and the mother 1/1 (or the reverse), the parental origin of each of the child's
alleles is known, and the child's alt-allele fraction indicates which parental copy is in excess, missing or doubled:
  gain  (cell fraction f)   alt fraction (1 + f) / (2 + f) if the duplicated copy carries the alt allele, else 1 / (2 + f)
  loss  (cell fraction f)   alt fraction (1 - f) / (2 - f) if the lost copy carries the alt allele, else 1 / (2 - f)
  LOH   (cell fraction f)   alt fraction (1 + f) / 2 if the retained copy carries the alt allele, else (1 - f) / 2
The log-likelihood ratio of the two parental assignments over the event's informative sites names the parent and
quantifies the support."""
import numpy as np

from .model import NA
from .vcfscan import GT_HET, GT_HOMALT, GT_HOMREF

EPS = 0.01          # sequencing noise at a homozygous site, so that a constitutional model never gives log(0)
GQ_ABSENT = -2      # a homozygous-reference call written by `bcftools merge -0` (vcfscan): no depth, taken as confident for a parent


def _confident(dp, gq, gt, min_dp, min_gq):
    """A call deep and sure enough, or a merge's homozygous-reference convention (no record in a per-sample VCF)."""
    return ((dp >= min_dp) & (gq >= min_gq)) | ((gt == GT_HOMREF) & (gq == GQ_ABSENT))


def _binom_ll(k, n, p):
    p = np.clip(p, EPS, 1 - EPS)
    return k * np.log(p) + (n - k) * np.log(1 - p)


def informative(sites, lo, hi, min_dp, min_gq):
    """Sites in [lo, hi] where the parents are opposite homozygotes and all three calls are deep and confident:
    (alt counts of the child, depths of the child, alt_from_mother bool)."""
    sel = (sites.pos >= lo) & (sites.pos <= hi)
    dp, alt, gt, gq = sites.dp[:, sel], sites.alt[:, sel], sites.gt[:, sel], sites.gq[:, sel]
    ok = (dp[0] >= min_dp) & (gq[0] >= min_gq) & _confident(dp[1], gq[1], gt[1], min_dp, min_gq) & _confident(dp[2], gq[2], gt[2], min_dp, min_gq)
    mat_alt = ok & (gt[1] == GT_HOMREF) & (gt[2] == GT_HOMALT)
    pat_alt = ok & (gt[1] == GT_HOMALT) & (gt[2] == GT_HOMREF)
    use = mat_alt | pat_alt
    return alt[0][use].astype(float), dp[0][use].astype(float), mat_alt[use]


def parent_of_origin(sites, ev, min_dp, min_gq):
    """(origin label, log-likelihood ratio, sites) for a child's gain, loss or LOH; f from the event (the bands' where the
    depth has none). The ratio is maternal-over-paternal for gains and LOH (which copy is extra / retained) and
    paternal-lost-over-maternal-lost for losses."""
    k, n, mat = informative(sites, ev.start + 1, ev.end, min_dp, min_gq)
    if len(n) < 10:
        return "", NA, int(len(n))
    f = ev.f if np.isfinite(ev.f) else 1.0
    f = float(np.clip(f, 0.02, 1.0))
    if ev.type == "gain":
        p_m = np.where(mat, (1 + f) / (2 + f), 1 / (2 + f))          # the extra copy maternal
        p_p = np.where(mat, 1 / (2 + f), (1 + f) / (2 + f))          # paternal
        llr = float((_binom_ll(k, n, p_m) - _binom_ll(k, n, p_p)).sum())
        return ("extra copy maternal" if llr > 0 else "extra copy paternal"), llr, int(len(n))
    if ev.type == "loss":
        p_pl = np.where(mat, 1 / (2 - f), (1 - f) / (2 - f))         # the paternal copy lost: the alt (maternal) stays whole
        p_ml = np.where(mat, (1 - f) / (2 - f), 1 / (2 - f))         # the maternal copy lost
        llr = float((_binom_ll(k, n, p_pl) - _binom_ll(k, n, p_ml)).sum())
        return ("paternal copy lost" if llr > 0 else "maternal copy lost"), llr, int(len(n))
    p_mr = np.where(mat, (1 + f) / 2, (1 - f) / 2)                   # the maternal copy retained (doubled)
    p_pr = np.where(mat, (1 - f) / 2, (1 + f) / 2)
    llr = float((_binom_ll(k, n, p_mr) - _binom_ll(k, n, p_pr)).sum())
    return ("maternal copy retained (paternal replaced)" if llr > 0 else "paternal copy retained (maternal replaced)"), llr, int(len(n))


def mie_rate(sites, lo, hi, min_dp, min_gq):
    """Mendelian errors among confident sites in [lo, hi]: a child homozygous for an allele a parent cannot give, or
    heterozygous where both parents are the same homozygote."""
    sel = (sites.pos >= lo) & (sites.pos <= hi)
    dp, gt, gq = sites.dp[:, sel], sites.gt[:, sel], sites.gq[:, sel]
    ok = (gt >= GT_HOMREF).all(axis=0)
    for m in range(3):
        ok &= _confident(dp[m], gq[m], gt[m], min_dp, min_gq)
    if ok.sum() < 20:
        return NA
    c, f, m = gt[0][ok], gt[1][ok], gt[2][ok]
    err = ((c == GT_HOMREF) & ((f == GT_HOMALT) | (m == GT_HOMALT))) | ((c == GT_HOMALT) & ((f == GT_HOMREF) | (m == GT_HOMREF))) | \
          ((c == GT_HET) & (((f == GT_HOMREF) & (m == GT_HOMREF)) | ((f == GT_HOMALT) & (m == GT_HOMALT))))
    return float(err.mean())


def read_trio(events, scan, genome, min_dp, min_gq, min_overlap=0.5):
    """Annotates the events in place: the child's with inheritance and the parent of origin; the parents' with whether the
    child carries the same event; every event with its Mendelian-error rate (the genome's baseline in the note)."""
    kid = [e for e in events if e.role == "child"]
    dad = [e for e in events if e.role == "father"]
    mom = [e for e in events if e.role == "mother"]
    base = {}
    for chrom, sites in scan.chroms.items():
        if chrom in genome.autosomes:
            r = mie_rate(sites, 1, genome.length[chrom], min_dp, min_gq)
            if np.isfinite(r):
                base[chrom] = r
    base_rate = float(np.median(list(base.values()))) if base else NA
    for e in events:
        sites = scan.sites(e.chrom)
        if sites is not None:
            e.mie_rate = mie_rate(sites, e.start + 1, e.end, min_dp, min_gq)
    for e in kid:
        same = [(p, src) for lst, src in ((dad, "father"), (mom, "mother")) for p in lst if p.chrom == e.chrom and p.type == e.type and e.overlap(p) >= min_overlap]
        if same:
            parts = []
            for p, src in same:
                parts.append("%s%s" % (src, " (in a share of the %s's cells)" % src if np.isfinite(p.f) and p.f < 0.8 else ""))
            e.inheritance = "inherited from the " + " and the ".join(dict.fromkeys(parts))
        else:
            e.inheritance = "new (neither parent carries it)"
        sites = scan.sites(e.chrom)
        if sites is not None and e.chrom in genome.autosomes:
            e.origin, e.origin_llr, e.origin_n = parent_of_origin(sites, e, min_dp, min_gq)
        if e.type == "LOH" and np.isfinite(e.f) and e.f >= 0.8:
            # in every cell: a uniparental isodisomy breaks Mendel at every site where the parents are opposite homozygotes (the child
            # is homozygous for one parent's allele); a run of homozygosity by descent has no such sites and no errors
            if not np.isfinite(e.mie_rate) or e.mie_rate < 0.05:
                e.note = (e.note + "; " if e.note else "") + "a run of homozygosity (both copies identical by descent), not a uniparental disomy: no Mendelian errors"
                e.origin, e.origin_llr = "", NA
            else:
                e.note = (e.note + "; " if e.note else "") + "Mendelian errors at the informative sites: a uniparental isodisomy"
    for lst, src in ((dad, "father"), (mom, "mother")):
        for p in lst:
            if p.type == "LOH" and np.isfinite(p.f) and p.f >= 0.8:
                p.note = (p.note + "; " if p.note else "") + "in every cell: a run of homozygosity (identical by descent) unless the depth says otherwise"
    for lst, src in ((dad, "father"), (mom, "mother")):
        for p in lst:
            hit = [e for e in kid if e.chrom == p.chrom and e.type == p.type and p.overlap(e) >= min_overlap]
            p.inheritance = "passed to the child" if hit else "not passed to the child"
    return base_rate
