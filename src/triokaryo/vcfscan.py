"""One pass over the trio's VCF: for every PASS biallelic SNV, per member, the depth, the alt-allele depth, the genotype
class and the genotype quality, stored per chromosome as arrays. No other record field is read, so an annotated VCF
costs no more than a bare one beyond its size. AD is required; DP falls back to the sum of AD; without GQ in the
header a called genotype with reads is taken as confident."""
from dataclasses import dataclass, field

import numpy as np
import pysam

GT_MISSING, GT_HOMREF, GT_HET, GT_HOMALT = -1, 0, 1, 2


@dataclass
class Sites:
    chrom: str
    pos: np.ndarray            # 1-based
    dp: np.ndarray             # (members, n) depth per member (DP, else the sum of AD)
    alt: np.ndarray            # (members, n) alt-allele depth (AD[1])
    gt: np.ndarray             # (members, n) genotype class
    gq: np.ndarray             # (members, n) genotype quality (-1 where absent; -2 a merge's homozygous reference)
    par: np.ndarray            # (n,) in a pseudoautosomal region

    @property
    def n(self):
        return len(self.pos)


@dataclass
class Scan:
    samples: tuple
    chroms: dict = field(default_factory=dict)      # chrom -> Sites, in the VCF's order
    n_records: int = 0
    n_used: int = 0
    skipped: dict = field(default_factory=dict)     # why records were left out
    gq_in_header: bool = True                       # False: the caller writes no GQ; a called genotype with reads is then taken as confident

    def sites(self, chrom):
        return self.chroms.get(chrom)


def _gt_class(gt):
    if gt is None:
        return GT_MISSING
    alleles = [a for a in gt if a is not None]
    if not alleles:
        return GT_MISSING
    if any(a > 1 for a in alleles):
        return GT_MISSING
    s = sum(alleles)
    if len(alleles) == 1:                       # haploid call (a male X): the class of its one allele
        return GT_HOMALT if s else GT_HOMREF
    return GT_HET if s == 1 else (GT_HOMALT if s == 2 else GT_HOMREF)


def scan_vcf(path, samples, genome, thin=1, pass_only=True, log=None, contigs=None):
    """samples: (child, father, mother) as named in the VCF. thin: keep every thin-th usable record (1: all).
    contigs: restrict to these (normalised names), else every chromosome of the genome."""
    vf = pysam.VariantFile(path)
    have = set(vf.header.samples)
    missing = [s for s in samples if s not in have]
    if missing:
        raise SystemExit("samples not in %s: %s (it holds: %s)" % (path, ", ".join(missing), ", ".join(list(have)[:6]) + (" ..." if len(have) > 6 else "")))
    vf.subset_samples(list(samples))
    formats = set(vf.header.formats.keys())
    if "AD" not in formats:
        raise SystemExit("%s: no AD (allelic depths) in FORMAT; triokaryo needs per-sample GT, AD and DP (GQ optional). Callers writing other "
                         "allelic-depth tags (freebayes AO/RO) need conversion" % path)
    if "DP" not in formats and log:
        log("WARNING: no DP in FORMAT; the depth is the sum of AD")
    has_gq = "GQ" in formats
    if not has_gq and log:
        log("WARNING: no GQ in FORMAT; a called genotype with reads is taken as confident on depth alone (--min-gq has no effect)")
    scan = Scan(samples=tuple(samples), gq_in_header=has_gq)
    skipped = scan.skipped
    buf = None
    cur = None
    cols = None

    def flush():
        if cur is None or not buf["pos"]:
            return
        n = len(buf["pos"])
        k = len(samples)
        pos = np.array(buf["pos"], dtype=np.int64)
        scan.chroms[cur] = Sites(cur, pos, np.array(buf["dp"], dtype=np.int32).reshape(n, k).T.copy(),
                                 np.array(buf["alt"], dtype=np.int32).reshape(n, k).T.copy(), np.array(buf["gt"], dtype=np.int8).reshape(n, k).T.copy(),
                                 np.array(buf["gq"], dtype=np.int16).reshape(n, k).T.copy(),
                                 np.array([genome.is_par(cur, p) for p in pos], dtype=bool) if cur in genome.par else np.zeros(n, dtype=bool))
        if log:
            log("  %s: %d sites" % (cur, n))

    k = 0
    for rec in vf:
        scan.n_records += 1
        chrom = genome.normalize(rec.contig)
        if chrom is None or (contigs and chrom not in contigs):
            skipped["contig"] = skipped.get("contig", 0) + 1
            continue
        alts = rec.alts
        if not alts or len(alts) != 1 or len(rec.ref) != 1 or len(alts[0]) != 1 or alts[0] == "*":
            skipped["not a biallelic SNV"] = skipped.get("not a biallelic SNV", 0) + 1
            continue
        if pass_only:
            f = rec.filter.keys()
            if f and "PASS" not in f:
                skipped["filtered"] = skipped.get("filtered", 0) + 1
                continue
        k += 1
        if thin > 1 and k % thin:
            skipped["thinned"] = skipped.get("thinned", 0) + 1
            continue
        if chrom != cur:
            flush()
            cur = chrom
            buf = dict(pos=[], dp=[], alt=[], gt=[], gq=[])
            if chrom in scan.chroms:
                raise SystemExit("%s: records of %s are not contiguous (sort the VCF)" % (path, chrom))
        buf["pos"].append(rec.pos)
        for s in rec.samples.values():
            ad = s.get("AD")
            a = int(ad[1]) if ad is not None and len(ad) > 1 and ad[1] is not None else 0
            d = s.get("DP")
            if d is None:
                d = sum(x for x in ad if x is not None) if ad is not None else 0
            q = s.get("GQ") if has_gq else None
            g = _gt_class(s.get("GT"))
            if not has_gq and g != GT_MISSING and (d or 0) > 0:
                q = 99                                     # no GQ anywhere: a called genotype with reads passes the confidence test
            # a per-sample VCF merged with `bcftools merge -0` writes 0/0 with no AD, DP or GQ where the sample had no record: a
            # homozygous reference call by convention, marked GQ -2 so that the trio reading can take it as confident
            if g == GT_HOMREF and d == 0 and q is None and (ad is None or all(x is None for x in ad)):
                q = -2
            buf["dp"].append(int(d))
            buf["alt"].append(a)
            buf["gt"].append(g)
            buf["gq"].append(int(q) if q is not None else -1)
        scan.n_used += 1
    flush()
    vf.close()
    return scan
