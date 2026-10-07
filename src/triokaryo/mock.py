"""A mock trio VCF with planted events, for the tests and for a demonstration: no real data anywhere in this package.

The genome is GRCh38 (the real lengths, a sparse site density). Each parent has two homologues per chromosome drawn from a
population allele frequency; the child inherits one from each. Events are written as changes to the copies of named
homologues in a share f of the cells, so that the depth, the B-allele bands, the heterozygosity rate, the parent of origin
and the inheritance all follow from one model. Depth is Poisson with a GC bias per sample (a synthetic GC track is written
beside the VCF); genotypes are called as a caller would (het where the alt fraction lies between 0.15 and 0.85)."""
import json
import os

import numpy as np
import pysam

from .genome import genome as load_genome

KID, DAD, MOM = "KID", "DAD", "MOM"
# the planted events: (member, chrom, start, end (None: chromosome end), f, kind, origin/inheritance note)
# kind: the homologue copy changes in the event cells - a dict of homologue -> delta copies, where homologues are
#   "pat" (the child's inherited paternal homologue), "pat_other" (the father's other, untransmitted one), "mat", "mat_other";
#   for a parent's own genome "h0"/"h1"
DEFAULT_EVENTS = [
    dict(member=KID, chrom="chr21", start=0, end=None, f=1.0, delta={"mat_other": +1}, label="trisomy 21, maternal (meiotic: both maternal homologues)", type="gain", origin="extra copy maternal", inherited=False),
    dict(member=KID, chrom="chr12", start=0, end=None, f=0.30, delta={"pat": +1}, label="mosaic +12 (30% of cells), paternal homologue duplicated", type="gain", origin="extra copy paternal", inherited=False),
    dict(member=DAD, chrom="chr10", start=39_800_000, end=None, f=1.0, delta={"h1": +1}, label="+10q in the father, constitutional", type="gain", origin="", inherited=False),
    dict(member=KID, chrom="chr10", start=39_800_000, end=None, f=1.0, delta={"pat": +1}, label="+10q inherited from the father", type="gain", origin="extra copy paternal", inherited=True),
    dict(member=KID, chrom="chr18", start=55_000_000, end=None, f=1.0, delta={"pat": -1}, label="del 18q, de novo, paternal copy lost", type="loss", origin="paternal copy lost", inherited=False),
    dict(member=KID, chrom="chr7", start=0, end=None, f=1.0, delta={"pat": -1, "mat": +1}, label="maternal isodisomy 7 (UPD, constitutional)", type="LOH", origin="maternal copy retained (paternal replaced)", inherited=False),
    dict(member=KID, chrom="chr6", start=0, end=59_800_000, f=0.40, delta={"pat": -1, "mat": +1}, label="mosaic CN-LOH 6p (40%), maternal copy retained", type="LOH", origin="maternal copy retained (paternal replaced)", inherited=False),
    dict(member=MOM, chrom="chr8", start=0, end=None, f=0.15, delta={"h0": +1}, label="mosaic +8 in the mother (15% of cells)", type="gain", origin="", inherited=False),
    dict(member=DAD, chrom="chr2", start=100_000_000, end=106_000_000, f=1.0, delta={"h1": -1}, label="6-Mb deletion in the father, constitutional", type="loss", origin="", inherited=False),
    dict(member=KID, chrom="chr2", start=100_000_000, end=106_000_000, f=1.0, delta={"pat": -1}, label="the father's 6-Mb deletion, inherited", type="loss", origin="paternal copy lost", inherited=True),
    dict(member=KID, chrom="chr15", start=0, end=None, f=1.0, delta={"pat": -1, "mat_other": +1}, label="maternal heterodisomy 15 (UPD, both maternal homologues, no paternal copy): the depth, the folded bands and the heterozygosity rate all flat", type="UPD", origin="both copies maternal (heterodisomy)", inherited=False),
]
# a share under the depth's threshold, for a dense mock (the phased scan needs the sites a real genome has: triokaryo mock --sites-per-mb 1000 --contigs chr15,chr16,chr17 --low-share)
LOW_SHARE_EVENTS = [
    dict(member=KID, chrom="chr16", start=0, end=None, f=0.08, delta={"pat": +1}, label="mosaic +16 (8% of cells), paternal homologue duplicated: under the depth's threshold, the phased bands' find", type="gain", origin="extra copy paternal", inherited=False),
    dict(member=MOM, chrom="chr16", start=0, end=36_800_000, f=0.08, delta={"h0": -1}, label="mosaic loss of 16p in the mother (8% of cells)", type="loss", origin="", inherited=False),
]
GC_CHROM = {"chr1": 0.0, "chr4": -0.03, "chr13": -0.025, "chr16": 0.03, "chr17": 0.04, "chr19": 0.07, "chr20": 0.02, "chr22": 0.06, "chrX": -0.02}


def _gt(alt, dp):
    if dp == 0:
        return "./."
    frac = alt / dp
    if alt >= 2 and 0.15 <= frac <= 0.85:
        return "0/1"
    return "1/1" if frac > 0.85 else "0/0"


def write_mock(out_dir, seed=1, sites_per_mb=60, no_events=False, xxy=False, depth=(30.0, 32.0, 28.0), gc_beta=(-0.8, -0.5, -1.0), bin_size=1_000_000,
               events=None, prefix="", contigs=None):
    """prefix: a tag before the sample names (KID, DAD, MOM), so that several mock trios can sit in one cohort.
    contigs: only these chromosomes (a dense small mock), else all."""
    os.makedirs(out_dir, exist_ok=True)
    nm = {KID: prefix + KID, DAD: prefix + DAD, MOM: prefix + MOM}
    rng = np.random.default_rng(seed)
    G = load_genome("grch38")
    events = [] if no_events else list(DEFAULT_EVENTS if events is None else events)
    # transmitted homologues per chromosome: the inherited paternal and maternal ones (chr10 and chr2 fixed so the father's events pass)
    tp = {c: int(rng.integers(2)) for c in G.chroms}
    tm = {c: int(rng.integers(2)) for c in G.chroms}
    tp["chr10"] = 1
    tp["chr2"] = 1
    # the GC track
    gc = {}
    with open(os.path.join(out_dir, "gc.tsv"), "w") as fh:
        fh.write("chrom\tstart\tend\tgc\n")
        for c in G.chroms:
            L = G.length[c]
            for s in range(0, L, bin_size):
                g = float(np.clip(0.41 + GC_CHROM.get(c, 0.0) + rng.normal(0, 0.035), 0.3, 0.62))
                gc[(c, s)] = g
                fh.write("%s\t%d\t%d\t%.4f\n" % (c, s, min(s + bin_size, L), g))
    truth = dict(seed=seed, samples=[nm[KID], nm[DAD], nm[MOM]], sexes=["M", "M", "F"], xxy=bool(xxy), events=[], transmitted_paternal=tp, transmitted_maternal=tm)
    vcf_txt = os.path.join(out_dir, "mock.vcf")
    with open(vcf_txt, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=triokaryo mock (no real data)\n")
        for c in G.chroms:
            fh.write("##contig=<ID=%s,length=%d>\n" % (c, G.length[c]))
        fh.write('##FILTER=<ID=PASS,Description="All filters passed">\n##FILTER=<ID=VQSRTrancheSNP99.90to100.00,Description="mock tranche">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
                 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s\t%s\n" % (nm[KID], nm[DAD], nm[MOM]))
        for c in G.chroms:
            if c == "chrY" or (contigs and c not in contigs):
                continue
            L = G.length[c]
            n = int(round(L / 1e6 * sites_per_mb))
            pos = np.sort(rng.integers(1, L, n))
            if c == "chrX":
                pos = pos[(pos > 2_781_479) & (pos < 155_701_383)]
                n = len(pos)
            p = np.clip(rng.beta(0.6, 0.6, n), 0.02, 0.98)
            # the four parental homologues
            d0, d1 = (rng.random(n) < p).astype(int), (rng.random(n) < p).astype(int)
            m0, m1 = (rng.random(n) < p).astype(int), (rng.random(n) < p).astype(int)
            hom = {"pat": [d0, d1][tp[c]], "pat_other": [d0, d1][1 - tp[c]], "mat": [m0, m1][tm[c]], "mat_other": [m0, m1][1 - tm[c]]}
            # base copies per member: the child's inherited pair; each parent's own two; the X: child (male) one maternal X, father one X
            base = {KID: {"pat": 1, "mat": 1}, DAD: {"h0": 1, "h1": 1}, MOM: {"h0": 1, "h1": 1}}
            if c == "chrX":
                base[KID] = {"mat": 1, "mat_other": 1} if xxy else {"mat": 1}
                base[DAD] = {"h0": 1}
            alleles = {KID: dict(hom), DAD: {"h0": d0, "h1": d1}, MOM: {"h0": m0, "h1": m1}}
            gcs = np.array([gc[(c, int(x // bin_size) * bin_size)] for x in pos])
            cols = {}
            for mi, member in enumerate((KID, DAD, MOM)):
                cp_base = base[member]
                tot_base = sum(cp_base.values())
                alt_base = sum(v * alleles[member][h] for h, v in cp_base.items())
                tot_ev = np.full(n, float(tot_base))
                alt_ev = alt_base.astype(float) if isinstance(alt_base, np.ndarray) else np.full(n, float(alt_base))
                f_site = np.zeros(n)
                for ev in events:
                    if ev["member"] != member or ev["chrom"] != c:
                        continue
                    e0, e1 = ev["start"], ev["end"] if ev["end"] is not None else L
                    sel = (pos >= e0) & (pos <= e1)
                    cp = dict(cp_base)
                    for h, dv in ev["delta"].items():
                        cp[h] = cp.get(h, 0) + dv
                    tot_ev[sel] = sum(cp.values())
                    alt_ev[sel] = sum(v * alleles[member][h][sel] for h, v in cp.items())
                    f_site[sel] = ev["f"]
                tot = f_site * tot_ev + (1 - f_site) * tot_base
                alt_mean = f_site * alt_ev + (1 - f_site) * (alt_base if isinstance(alt_base, np.ndarray) else float(alt_base))
                p_alt = np.where(tot > 0, alt_mean / np.maximum(tot, 1e-9), 0.0)
                p_alt = p_alt * (1 - 0.003) + (1 - p_alt) * 0.003
                mean_dp = depth[mi] * np.exp(gc_beta[mi] * (gcs - 0.41)) * tot / 2.0
                dp = rng.poisson(mean_dp)
                alt = rng.binomial(dp, p_alt)
                gq = np.where(rng.random(n) < 0.08, rng.integers(2, 20, n), 99)
                cols[member] = (dp, alt, gq)
            for i in range(n):
                filt = "PASS"
                u = rng.random()
                ref, alt_allele = "A", "G"
                if u < 0.01:
                    filt = "VQSRTrancheSNP99.90to100.00"
                elif u < 0.015:
                    ref, alt_allele = "AT", "A"                     # an indel: skipped by the scan
                elif u < 0.018:
                    alt_allele = "G,T"                              # multiallelic: skipped
                fields = []
                for member in (KID, DAD, MOM):
                    dp, al, gq = cols[member]
                    d, a, q = int(dp[i]), int(al[i]), int(gq[i])
                    fields.append("%s:%d,%d:%d:%s" % (_gt(a, d), d - a, a, d, q if d else "."))
                fh.write("%s\t%d\t.\t%s\t%s\t100\t%s\t.\tGT:AD:DP:GQ\t%s\n" % (c, int(pos[i]), ref, alt_allele, filt, "\t".join(fields)))
    vcf_gz = vcf_txt + ".gz"
    pysam.tabix_compress(vcf_txt, vcf_gz, force=True)
    pysam.tabix_index(vcf_gz, preset="vcf", force=True)
    os.remove(vcf_txt)
    for ev in events:
        truth["events"].append(dict(sample=nm[ev["member"]], chrom=ev["chrom"], start=ev["start"], end=ev["end"] if ev["end"] is not None else G.length[ev["chrom"]],
                                    f=ev["f"], type=ev["type"], label=ev["label"], origin=ev["origin"], inherited=ev["inherited"]))
    truth["x_copies"] = {nm[KID]: 2 if xxy else 1, nm[DAD]: 1, nm[MOM]: 2}
    with open(os.path.join(out_dir, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=1)
    with open(os.path.join(out_dir, "mock.trios.tsv"), "w") as fh:
        fh.write("#kid\tdad\tmom\tkid_sex\tdad_sex\tmom_sex\n%s\t%s\t%s\t1\t1\t2\n" % (nm[KID], nm[DAD], nm[MOM]))
    # the planted events as another caller would list them (NGS-DOSE's karyotype events columns), for the concordance check
    with open(os.path.join(out_dir, "events.external.tsv"), "w") as fh:
        fh.write("sample\tchrom\tspan\tstart_mb\tend_mb\tlabel\tkind\n")
        for ev in truth["events"]:
            if ev["type"] in ("LOH", "UPD") or ev["f"] < 0.10:
                continue                                           # a depth tool cannot see a copy-neutral event, nor a share under 10%
            L = G.length[ev["chrom"]]
            pe = G.p_end[ev["chrom"]]
            span = "whole" if ev["start"] == 0 and ev["end"] == L else "q" if ev["start"] == pe and ev["end"] == L else "p" if ev["start"] == 0 and ev["end"] == pe else "stretch"
            fh.write("%s\t%s\t%s\t%.3f\t%.3f\t%s%s%s\t%s\n" % (ev["sample"], ev["chrom"], span, ev["start"] / 1e6, ev["end"] / 1e6, "+" if ev["type"] == "gain" else "-",
                                                             ev["chrom"][3:], ("[%.2f]" % ev["f"]) if ev["f"] < 1 else "", "gain" if ev["type"] == "gain" else "loss"))
    with open(os.path.join(out_dir, "MOCK_DATA.txt"), "w") as fh:
        fh.write("MOCK DATA - triokaryo mock (seed %d): every record simulated, no real genome.\n" % seed)
    return dict(vcf=vcf_gz, trios=os.path.join(out_dir, "mock.trios.tsv"), truth=os.path.join(out_dir, "truth.json"), gc=os.path.join(out_dir, "gc.tsv"),
                events=os.path.join(out_dir, "events.external.tsv"))
