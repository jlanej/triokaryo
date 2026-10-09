"""triokaryo: large chromosomal events in a trio from its VCF.

  triokaryo run --vcf trio.vcf.gz --pedigree trios.tsv --child KID --out out/KID [--panel 1kg-dragen] [--gc-track gc.tsv|none] [--events other_calls.tsv]
  triokaryo run --vcf trio.vcf.gz --child KID --father DAD --mother MOM --sex M,M,F --out out/KID
  triokaryo merge --child kid.vcf.gz --father dad.vcf.gz --mother mom.vcf.gz --out trio.vcf.gz   # per-sample VCFs, with bcftools
  triokaryo panel --vcfs a.vcf.gz b.vcf.gz ... --out panel.tsv        # or --runs 'out/*' from earlier runs
  triokaryo gc-track --fasta ref.fa --out gc.tsv [--bin 1000000]
  triokaryo mock --out mock_dir [--seed 1] [--no-events] [--xxy] [--contigs chr15,chr16,chr17 --sites-per-mb 1000 --low-share | --meiosis]
  triokaryo batch --pedigree trios.tsv --vcf joint.vcf.gz --out out --jobs 8 --cohort cohort [--panel ...]   # every trio, in parallel
  triokaryo cohort --runs 'out/*' --out cohort [--events other_calls.tsv]   # cohort report and guide
  triokaryo calibrate --out calib [--cell-fractions 0.05,0.1,0.2,0.3,0.5,1 --depths 15,30,60 --sizes 5,10,20,50 --replicates 2]
  triokaryo report --runs 'out/*'                   # rebuild a run's page from its tables (no VCF needed)
  triokaryo guide --out guide.html                  # the meaning of every figure row, colour and column
"""
import argparse
import glob
import html
import json
import os
import sys

from . import __version__


def _log(s):
    print("[triokaryo] " + s, file=sys.stderr)


BUILT_IN_PANELS = {"1kg-dragen": "panel.1kg_dragen_3.7.6.1mb.tsv"}
BUILT_IN_GC_TRACKS = {"hg38": "gc.hg38.1mb.tsv", "grch38": "gc.hg38.1mb.tsv"}


def _data_path(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name)


def resolve_panel(panel):
    """A panel file, or a built-in name."""
    if panel and not os.path.exists(panel):
        if panel in BUILT_IN_PANELS:
            return _data_path(BUILT_IN_PANELS[panel])
        sys.exit("no panel %s (a file, or one of: %s)" % (panel, ", ".join(BUILT_IN_PANELS)))
    return panel


def resolve_gc_track(gc_track, genome_name, bin_size):
    """A GC track file, a built-in name ('hg38'), 'none', or, when not given, the shipped hg38 track where the genome and the
    1-Mb bin match it (otherwise no correction: `triokaryo gc-track` writes one for another reference or bin)."""
    if gc_track is None:
        return _data_path(BUILT_IN_GC_TRACKS["hg38"]) if genome_name.lower() in BUILT_IN_GC_TRACKS and bin_size == 1_000_000 else None
    if gc_track.lower() == "none":
        return None
    if os.path.exists(gc_track):
        return gc_track
    if gc_track.lower() in BUILT_IN_GC_TRACKS:
        return _data_path(BUILT_IN_GC_TRACKS[gc_track.lower()])
    sys.exit("no GC track %s (a file, 'hg38' for the shipped 1-Mb track, or 'none')" % gc_track)


def cmd_run(a):
    from .pedigree import read_trios, trio_from_args
    from .pipeline import run_trio
    if a.pedigree:
        trios = read_trios(a.pedigree)
        if a.child:
            trios = [t for t in trios if t.kid == a.child]
        if not trios:
            sys.exit("no trio %s in %s" % (a.child or "", a.pedigree))
        if len(trios) > 1 and not a.child:
            sys.exit("%d trios in %s: name one with --child" % (len(trios), a.pedigree))
        trio = trios[0]
    else:
        if not (a.child and a.father and a.mother):
            sys.exit("give --pedigree, or --child --father --mother")
        trio = trio_from_args(a.child, a.father, a.mother, a.sex)
    params = dict(min_abs=a.min_abs, z=a.z, min_len=a.min_len, min_f=a.min_f)
    a.panel = resolve_panel(a.panel)
    a.gc_track = resolve_gc_track(a.gc_track, a.genome, a.bin)
    note = ""
    for p in (os.path.join(os.path.dirname(os.path.abspath(a.vcf)), "MOCK_DATA.txt"),):
        if os.path.exists(p):
            note = open(p).read().strip()
    run_trio(a.vcf, trio, a.out, gc_track=a.gc_track, events_path=a.events, bin_size=a.bin, min_dp=a.min_dp, min_gq=a.min_gq, thin=a.thin,
             genome_name=a.genome, figures=not a.no_figures, params=params, log=_log, mock_note=note, panel=a.panel, genotypes=a.genotypes, depth_sites=a.depth_sites)
    return 0


def cmd_panel(a):
    from .genome import genome
    from .panel import build_panel, write_panel
    runs = sorted(d for pat in (a.runs or []) for d in glob.glob(pat))
    roles = tuple(r.strip() for r in a.roles.split(",") if r.strip()) if a.roles else None
    if roles and any(r not in ("child", "father", "mother") for r in roles):
        sys.exit("--roles takes child, father, mother (comma-separated); got %s" % a.roles)
    rows, samples, n_y = build_panel(genome(a.genome), a.bin, vcfs=a.vcfs or (), samples=set(a.samples.split(",")) if a.samples else None, runs=runs, log=_log,
                                     thin=a.thin, roles=roles)
    if not rows:
        sys.exit("no genome read for the panel (--vcfs and/or --runs)")
    write_panel(a.out, rows, samples, a.bin)
    import numpy as np
    from .model import PANEL_MAX_RSD, PANEL_MIN_N, PANEL_MIN_N_Y
    unp = sum(1 for r in rows if r["n"] < (PANEL_MIN_N_Y if r["chrom"] == "chrY" else PANEL_MIN_N) or (np.isfinite(r["lrr_rsd"]) and r["lrr_rsd"] > PANEL_MAX_RSD))
    if len(samples) < PANEL_MIN_N:
        _log("WARNING: %d genomes: a panel needs %d or more, or every bin is left out of the calls" % (len(samples), PANEL_MIN_N))
    if n_y < PANEL_MIN_N_Y:
        _log("WARNING: %d genome(s) with a Y: the Y rows need %d or more (males), or the Y is not corrected and mosaic loss of Y is read from the "
             "father/son ratio only" % (n_y, PANEL_MIN_N_Y))
    _log("panel of %d genomes (%d with a Y), %d bins (%d masked) -> %s" % (len(samples), n_y, len(rows), unp, a.out))
    return 0


def cmd_gc_track(a):
    from .genome import genome
    from .gctrack import write_gc_track
    n = write_gc_track(a.fasta, a.out, genome(a.genome), a.bin)
    _log("%d bins -> %s" % (n, a.out))
    return 0


def cmd_mock(a):
    from .mock import write_mock
    events = None
    if a.low_share:
        from .mock import LOW_SHARE_EVENTS
        events = LOW_SHARE_EVENTS
    if a.meiosis:
        from .mock import MEIOSIS_EVENTS
        events = (events or []) + MEIOSIS_EVENTS
    if a.sex_chromosomes:
        from .mock import SEX_EVENTS
        events = (events or []) + SEX_EVENTS
    deficit = tuple(float(x) for x in a.sex_deficit.split(",")) if a.sex_deficit else (1.0, 1.0)
    switches = {"mat": {}, "pat": {}}
    for item in (x for x in a.switch.split(",") if x.strip()):
        parts = item.split(":")
        if len(parts) != 3 or parts[0] not in switches:
            sys.exit("--switch: give mat:CHROM:POS or pat:CHROM:POS, comma-separated (got %r)" % item)
        switches[parts[0]][parts[1]] = int(float(parts[2]))
    paths = write_mock(a.out, seed=a.seed, sites_per_mb=a.sites_per_mb, no_events=a.no_events, xxy=a.xxy, prefix=a.prefix, events=events,
                       contigs=a.contigs.split(",") if a.contigs else None, child_sex=a.child_sex, xxx=a.xxx, sex_deficit=deficit,
                       switch_maternal=switches["mat"] or None, switch_paternal=switches["pat"] or None, refined=a.refined, ref_blocks=a.ref_blocks)
    _log("mock trio -> %s" % paths["vcf"])
    return 0


def _batch_one(args):
    """One trio of a batch, in a worker process: returns (kid, out dir, error or None)."""
    vcf, trio, out, kw = args
    try:
        from .pipeline import run_trio
        run_trio(vcf, trio, out, log=lambda s: None, **kw)
        return trio.kid, out, None
    except BaseException as e:  # noqa: BLE001  (a failed trio is reported, the batch goes on)
        return trio.kid, out, "%s: %s" % (type(e).__name__, e)


def cmd_batch(a):
    """Every trio of a pedigree, in parallel, from one joint VCF or per-trio VCFs named by a pattern; then the cohort report."""
    import multiprocessing
    from .pedigree import read_trios
    trios = read_trios(a.pedigree)
    if a.child:
        keep = set(a.child.split(","))
        trios = [t for t in trios if t.kid in keep]
    if not trios:
        sys.exit("no trio in %s" % a.pedigree)
    if not (a.vcf or a.vcf_pattern):
        sys.exit("give --vcf (one joint VCF) or --vcf-pattern (per-trio VCFs, e.g. 'vcfs/{kid}.vcf.gz')")
    panel = resolve_panel(a.panel)
    kw = dict(gc_track=resolve_gc_track(a.gc_track, a.genome, a.bin), events_path=a.events, bin_size=a.bin, min_dp=a.min_dp, min_gq=a.min_gq, thin=a.thin, genome_name=a.genome,
              figures=not a.no_figures, panel=panel, genotypes=a.genotypes, depth_sites=a.depth_sites)
    jobs = []
    for t in trios:
        vcf = a.vcf or a.vcf_pattern.format(kid=t.kid, dad=t.dad, mom=t.mom, family=t.family)
        if not os.path.exists(vcf):
            _log("WARNING: %s: no VCF at %s, skipped" % (t.kid, vcf))
            continue
        jobs.append((vcf, t, os.path.join(a.out, t.kid), kw))
    _log("%d trios, %d worker(s)" % (len(jobs), a.jobs))
    done, failed = [], []
    if a.jobs > 1:
        with multiprocessing.get_context("spawn").Pool(a.jobs) as pool:
            results = pool.imap_unordered(_batch_one, jobs)
            for kid, out, err in results:
                (failed if err else done).append((kid, out, err))
                _log("%s: %s" % (kid, err or ("-> " + out)))
    else:
        for job in jobs:
            kid, out, err = _batch_one(job)
            (failed if err else done).append((kid, out, err))
            _log("%s: %s" % (kid, err or ("-> " + out)))
    if failed:
        _log("%d trio(s) failed: %s" % (len(failed), "; ".join("%s (%s)" % (k, e) for k, _, e in failed)))
    if a.cohort and done:
        from .report import write_cohort
        write_cohort(a.cohort, [out for _, out, _ in done], events_path=a.events, genome_name=a.genome, log=_log)
    return 1 if failed and not done else 0


def cmd_merge(a):
    from .merge import merge_trio
    out = merge_trio([a.child, a.father, a.mother], a.out, log=_log, keep_intermediate=a.keep_intermediate, threads=a.threads)
    _log("trio VCF -> %s" % out)
    return 0


def cmd_calibrate(a):
    from .calibrate import run_calibration
    fl = lambda s: tuple(float(x) for x in s.split(",") if x.strip())  # noqa: E731
    rows = run_calibration(a.out, fractions=fl(a.cell_fractions), depths=fl(a.depths), sizes=fl(a.sizes), replicates=a.replicates, seed=a.seed,
                           sites_per_mb=a.sites_per_mb, log=_log, figures=not a.no_figures, contigs=a.contigs.split(",") if a.contigs else None, whole=a.whole)
    _log("%d planted events, %d detected -> %s" % (len(rows), sum(1 for r in rows if r["detected"]), a.out))
    return 0


def cmd_cohort(a):
    from .report import write_cohort
    runs = sorted(d for pat in a.runs for d in glob.glob(pat) if os.path.exists(os.path.join(d, "summary.json")))
    if not runs:
        sys.exit("no run with a summary.json under: " + " ".join(a.runs))
    write_cohort(a.out, runs, events_path=a.events, genome_name=a.genome, log=_log)
    return 0


def cmd_digest(a):
    from .digest import write_digest
    runs = sorted(d for pat in a.runs for d in glob.glob(pat) if os.path.exists(os.path.join(d, "summary.json")))
    if not runs:
        sys.exit("no run with a summary.json under: " + " ".join(a.runs))
    write_digest(a.out, runs, events_path=a.events, genome_name=a.genome, top=a.top, min_mb=a.min_mb, min_f=a.min_f, min_mb_loh=a.min_mb_loh,
                 parent_sex_mosaics=a.parent_sex_mosaics, allow_first_pass=a.allow_first_pass, max_width=a.max_width, title=a.title, log=_log)
    return 0


def cmd_guide(a):
    from .guide import write_guide
    write_guide(a.out, figures_dir=a.figures)
    _log("the guide -> %s" % a.out)
    return 0


def cmd_report(a):
    from .report import rebuild_run
    runs = sorted(d for pat in a.runs for d in glob.glob(pat) if os.path.exists(os.path.join(d, "summary.json")))
    if not runs:
        sys.exit("no run with a summary.json under: " + " ".join(a.runs))
    for d in runs:
        rebuild_run(d, log=_log)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="triokaryo", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version="triokaryo " + __version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="analyse one trio: scan the VCF, bin, segment, call, phase, write tables, figures and the page")
    r.add_argument("--vcf", required=True)
    r.add_argument("--pedigree", help="a trios file (#kid dad mom [kid_sex dad_sex mom_sex]) or a PED file")
    r.add_argument("--child")
    r.add_argument("--father")
    r.add_argument("--mother")
    r.add_argument("--sex", default="", help="sexes of child,father,mother as M/F when the members are named directly, e.g. M,M,F")
    r.add_argument("--out", required=True)
    r.add_argument("--gc-track", help="per-bin GC fraction the LRR is corrected against: a file (triokaryo gc-track), 'hg38' (the shipped 1-Mb track, "
                                      "the default where the genome is GRCh38 and the bin 1 Mb), or 'none'")
    r.add_argument("--panel", help="reference panel (triokaryo panel): per-bin median LRR, band deviation and heterozygosity rate of other genomes, "
                                   "subtracted from each member's tracks; recommended on real data, where centromere flanks and segmental duplications "
                                   "are otherwise called as events. A file, or '1kg-dragen' (twelve public 1000 Genomes genomes called by DRAGEN 3.7.6, 1-Mb bins)")
    r.add_argument("--events", help="events from another method to draw and match (sample chrom start end label [type]; or NGS-DOSE's karyotype/events.tsv)")
    r.add_argument("--genotypes", default="auto", choices=["auto", "pl", "vcf"], help="genotypes: as written in the VCF ('vcf'), re-derived from PL ('pl'), or re-derived where the header declares PP ('auto', the default): posterior genotypes (GATK CalculateGenotypePosteriors) were refined under a pedigree prior that hides the Mendelian errors the trio analysis reads and treats a son's X as diploid")
    r.add_argument("--depth-sites", default="auto", choices=["auto", "all", "variant"], help="the genotypes the bin depth is taken over: every one with reads ('all'), each member's own heterozygous and homozygous-alternate ones ('variant'), or the latter where homozygous-reference genotypes carry a reference block's depth ('auto', the default: MIN_DP declared, or hom-ref depth under 0.95 of the heterozygous depth, as GATK GenotypeGVCFs writes)")
    r.add_argument("--bin", type=int, default=1_000_000, help="bin width in bp")
    r.add_argument("--min-dp", type=int, default=8, help="minimum depth of a confident call")
    r.add_argument("--min-gq", type=int, default=20, help="minimum GQ of a confident call")
    r.add_argument("--thin", type=int, default=1, help="use every n-th usable site (faster; 1 = all)")
    r.add_argument("--genome", default="grch38")
    r.add_argument("--min-abs", type=float, default=0.07, help="minimum |mean LRR| of a gain or loss (0.07: a cell fraction of about 10%%)")
    r.add_argument("--z", type=float, default=5.0, help="minimum split statistic for a segment boundary (binary segmentation)")
    r.add_argument("--min-len", type=int, default=5, help="minimum bins per segment")
    r.add_argument("--min-f", type=float, default=0.10, help="minimum cell fraction reported")
    r.add_argument("--no-figures", action="store_true")
    r.set_defaults(fn=cmd_run)
    pn = sub.add_parser("panel", help="build a reference panel: per-bin median LRR, band deviation and heterozygosity rate over other genomes (VCFs or earlier runs)")
    pn.add_argument("--vcfs", nargs="*", help="VCFs to read (every sample of each unless --samples)")
    pn.add_argument("--samples", help="comma-separated sample names to take from the VCFs")
    pn.add_argument("--runs", nargs="*", help="earlier run directories (bins.tsv + summary.tsv), e.g. 'out/*'")
    pn.add_argument("--roles", default="", help="with --runs: the members to take from each run, comma-separated (child, father, mother; default all three); "
                                               "e.g. father,child for a panel whose Y rows come from every male")
    pn.add_argument("--out", required=True)
    pn.add_argument("--bin", type=int, default=1_000_000)
    pn.add_argument("--thin", type=int, default=1)
    pn.add_argument("--genome", default="grch38")
    pn.set_defaults(fn=cmd_panel)
    g = sub.add_parser("gc-track", help="per-bin GC fraction from a reference FASTA")
    g.add_argument("--fasta", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--bin", type=int, default=1_000_000)
    g.add_argument("--genome", default="grch38")
    g.set_defaults(fn=cmd_gc_track)
    m = sub.add_parser("mock", help="a simulated trio VCF with planted events (no real data)")
    m.add_argument("--out", required=True)
    m.add_argument("--seed", type=int, default=1)
    m.add_argument("--sites-per-mb", type=int, default=60)
    m.add_argument("--no-events", action="store_true")
    m.add_argument("--child-sex", default="M", choices=["M", "F"], help="the child's sex (default M)")
    m.add_argument("--xxy", action="store_true", help="a 47,XXY son with both maternal X homologues (a maternal meiosis I error)")
    m.add_argument("--xxx", action="store_true", help="a 47,XXX daughter with both maternal X homologues and the paternal X (a maternal meiosis I error); with --child-sex F")
    m.add_argument("--prefix", default="", help="prefix for the sample names KID, DAD, MOM (several simulated trios in one cohort)")
    m.add_argument("--contigs", default="", help="restrict to these chromosomes, comma-separated (a small dense simulation)")
    m.add_argument("--low-share", action="store_true", help="plant the low-cell-fraction events (below the depth threshold) instead of the default set")
    m.add_argument("--meiosis", action="store_true", help="plant the meiotic-stage trisomies (meiosis I and II with a crossover, mitotic) on chr13, chr16, chr17; use with --contigs chr10,chr11,chr12,chr13,chr16,chr17")
    m.add_argument("--sex-chromosomes", action="store_true", help="plant the sex-chromosome mosaics: loss of Y in the father, 45,X/46,XX in the mother, 46,XY/47,XXY in the son")
    m.add_argument("--sex-deficit", default="", help="depth factors on the X and the Y imitating real data's mappability deficit, e.g. 0.93,0.90")
    m.add_argument("--switch", default="", help="crossovers in the transmitted haplotypes, e.g. mat:chr8:70000000,pat:chr2:50000000: a parent's event changes sign there along the phased track")
    m.add_argument("--refined", action="store_true", help="write PL and PP, with GT and GQ refined under a pedigree prior as GATK CalculateGenotypePosteriors does (Mendelian violations penalised by 80 phred)")
    m.add_argument("--ref-blocks", action="store_true", help="write homozygous-reference genotypes with a reference block's depth (about 18%% under the site's), as GATK GenotypeGVCFs does")
    m.set_defaults(fn=cmd_mock)
    dg = sub.add_parser("digest", help="one self-contained page of the cohort's largest events, each with the trio's genome-wide and chromosome figures: the hand-off")
    dg.add_argument("--runs", nargs="+", required=True, help="run directories (globs)")
    dg.add_argument("--out", required=True, help="output directory: digest.html, digest.tsv, digest_left_out.tsv")
    dg.add_argument("--events", help="events from another method to match (NGS-DOSE's karyotype/events.tsv), where the runs were not made with them")
    dg.add_argument("--top", type=int, default=25, help="how many events to show (default 25)")
    dg.add_argument("--min-mb", type=float, default=10.0, help="a stretch under this many Mb is left out (default 10)")
    dg.add_argument("--min-f", type=float, default=0.2, help="an arm or stretch in under this share of the cells is left out (default 0.2); a whole chromosome never is")
    dg.add_argument("--min-mb-loh", type=float, default=20.0, help="a run of homozygosity (LOH stretch) under this many Mb is left out (default 20)")
    dg.add_argument("--parent-sex-mosaics", action="store_true", help="list a parent's X or Y lost in part of the cells (loss with age) too")
    dg.add_argument("--allow-first-pass", action="store_true", help="keep runs made without a panel (the first pass, whose calls are not to be read)")
    dg.add_argument("--max-width", type=int, default=1500, help="embedded figures downscaled to this width in pixels (default 1500)")
    dg.add_argument("--title", default="", help="the page's title")
    dg.add_argument("--genome", default="grch38")
    dg.set_defaults(fn=cmd_digest)
    c = sub.add_parser("cohort", help="cohort report over many runs: counts, landscape figure, every event, per-trio metrics, concordance, guide")
    c.add_argument("--runs", nargs="+", required=True, help="run directories (globs)")
    c.add_argument("--out", required=True)
    c.add_argument("--events", help="events from another method to match (as for run)")
    c.add_argument("--genome", default="grch38")
    c.set_defaults(fn=cmd_cohort)
    b = sub.add_parser("batch", help="every trio of a pedigree in parallel, from one joint VCF or per-trio VCFs, then the cohort report")
    b.add_argument("--pedigree", required=True, help="a trios file or PED with every trio")
    b.add_argument("--vcf", help="one joint VCF holding every member")
    b.add_argument("--vcf-pattern", help="per-trio VCF paths with {kid}, {dad}, {mom} or {family}, e.g. 'vcfs/{kid}.trio.vcf.gz'")
    b.add_argument("--out", required=True, help="output directory; each trio under <out>/<child>")
    b.add_argument("--child", help="comma-separated children to run (default all)")
    b.add_argument("--jobs", type=int, default=1, help="parallel workers")
    b.add_argument("--cohort", help="write the cohort report over the finished trios into this directory")
    b.add_argument("--gc-track", help="as for run: a file, 'hg38' (the default at 1 Mb on GRCh38) or 'none'")
    b.add_argument("--panel", help="as for run: a file or '1kg-dragen'")
    b.add_argument("--events")
    b.add_argument("--genotypes", default="auto", choices=["auto", "pl", "vcf"], help="as for run")
    b.add_argument("--depth-sites", default="auto", choices=["auto", "all", "variant"], help="as for run")
    b.add_argument("--bin", type=int, default=1_000_000)
    b.add_argument("--min-dp", type=int, default=8)
    b.add_argument("--min-gq", type=int, default=20)
    b.add_argument("--thin", type=int, default=1)
    b.add_argument("--genome", default="grch38")
    b.add_argument("--no-figures", action="store_true")
    b.set_defaults(fn=cmd_batch)
    mg = sub.add_parser("merge", help="per-sample VCFs into one trio VCF with bcftools: PASS biallelic SNVs, merged with -0 (absent = homozygous reference), indexed")
    mg.add_argument("--child", required=True, help="the child's VCF (bgzipped and indexed)")
    mg.add_argument("--father", required=True)
    mg.add_argument("--mother", required=True)
    mg.add_argument("--out", required=True, help="the trio VCF to write (.vcf.gz)")
    mg.add_argument("--threads", type=int, default=1)
    mg.add_argument("--keep-intermediate", action="store_true", help="keep the per-sample SNV VCFs")
    mg.set_defaults(fn=cmd_merge)
    cb = sub.add_parser("calibrate", help="detection and cell-fraction accuracy on the simulator over a grid of cell fractions, depths and event sizes")
    cb.add_argument("--out", required=True)
    cb.add_argument("--cell-fractions", default="0.05,0.1,0.2,0.3,0.5,1", help="comma-separated cell fractions (default 0.05,0.1,0.2,0.3,0.5,1)")
    cb.add_argument("--depths", default="30", help="comma-separated mean depths (default 30)")
    cb.add_argument("--sizes", default="5,10,20,50", help="comma-separated event sizes in Mb (default 5,10,20,50)")
    cb.add_argument("--replicates", type=int, default=1)
    cb.add_argument("--seed", type=int, default=1)
    cb.add_argument("--sites-per-mb", type=int, default=60, help="simulated site density (real WGS: about 1000 per Mb; 60 is quick but understates the phased scan)")
    cb.add_argument("--contigs", default="", help="restrict the simulation to these chromosomes, e.g. chr1,...,chr16 with a dense --sites-per-mb (chr1-18 with --whole)")
    cb.add_argument("--whole", action="store_true", help="add a whole-chromosome maternal meiosis I trisomy (chr17) and a maternal heterodisomy (chr18) at each cell fraction")
    cb.add_argument("--no-figures", action="store_true")
    cb.set_defaults(fn=cmd_calibrate)
    g = sub.add_parser("guide", help="the guide: the meaning of every figure row, colour, call and column, with pattern cards (one self-contained page)")
    g.add_argument("--out", required=True, help="the HTML file to write")
    g.add_argument("--figures", help="keep the pattern figures (PNG, SVG, PDF, sidecars, legends) in this directory")
    g.set_defaults(fn=cmd_guide)
    rp = sub.add_parser("report", help="rebuild a run's page, sidecars, legends and guide from its tables and figures (no VCF needed)")
    rp.add_argument("--runs", nargs="+", required=True, help="run directories (globs)")
    rp.set_defaults(fn=cmd_report)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
