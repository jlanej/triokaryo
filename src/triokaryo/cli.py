"""triokaryo: large chromosomal events in a trio from its VCF.

  triokaryo run --vcf trio.vcf.gz --pedigree trios.tsv --child KID --out out/KID [--gc-track gc.tsv] [--panel panel.tsv] [--events other_calls.tsv]
  triokaryo run --vcf trio.vcf.gz --child KID --father DAD --mother MOM --sex M,M,F --out out/KID
  triokaryo panel --vcfs a.vcf.gz b.vcf.gz ... --out panel.tsv        # or --runs 'out/*' from earlier runs
  triokaryo gc-track --fasta ref.fa --out gc.tsv [--bin 1000000]
  triokaryo mock --out mock_dir [--seed 1] [--no-events] [--xxy] [--contigs chr15,chr16,chr17 --sites-per-mb 1000 --low-share | --meiosis]
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
    if a.panel and not os.path.exists(a.panel):
        built_in = {"1kg-dragen": "panel.1kg_dragen_3.7.6.1mb.tsv"}
        if a.panel in built_in:
            a.panel = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", built_in[a.panel])
        else:
            sys.exit("no panel %s (a file, or one of: %s)" % (a.panel, ", ".join(built_in)))
    note = ""
    for p in (os.path.join(os.path.dirname(os.path.abspath(a.vcf)), "MOCK_DATA.txt"),):
        if os.path.exists(p):
            note = open(p).read().strip()
    run_trio(a.vcf, trio, a.out, gc_track=a.gc_track, events_path=a.events, bin_size=a.bin, min_dp=a.min_dp, min_gq=a.min_gq, thin=a.thin,
             genome_name=a.genome, figures=not a.no_figures, params=params, log=_log, mock_note=note, panel=a.panel)
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
    paths = write_mock(a.out, seed=a.seed, sites_per_mb=a.sites_per_mb, no_events=a.no_events, xxy=a.xxy, prefix=a.prefix, events=events,
                       contigs=a.contigs.split(",") if a.contigs else None, child_sex=a.child_sex, xxx=a.xxx, sex_deficit=deficit)
    _log("mock trio -> %s" % paths["vcf"])
    return 0


def cmd_calibrate(a):
    from .calibrate import run_calibration
    fl = lambda s: tuple(float(x) for x in s.split(",") if x.strip())  # noqa: E731
    rows = run_calibration(a.out, fractions=fl(a.cell_fractions), depths=fl(a.depths), sizes=fl(a.sizes), replicates=a.replicates, seed=a.seed,
                           sites_per_mb=a.sites_per_mb, log=_log, figures=not a.no_figures, contigs=a.contigs.split(",") if a.contigs else None)
    _log("%d planted events, %d detected -> %s" % (len(rows), sum(1 for r in rows if r["detected"]), a.out))
    return 0


def cmd_cohort(a):
    from .report import write_cohort
    runs = sorted(d for pat in a.runs for d in glob.glob(pat) if os.path.exists(os.path.join(d, "summary.json")))
    if not runs:
        sys.exit("no run with a summary.json under: " + " ".join(a.runs))
    write_cohort(a.out, runs, events_path=a.events, genome_name=a.genome, log=_log)
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
    r.add_argument("--gc-track", help="per-bin GC fraction (triokaryo gc-track); the LRR is GC-corrected against it")
    r.add_argument("--panel", help="reference panel (triokaryo panel): per-bin median LRR, band deviation and heterozygosity rate of other genomes, "
                                   "subtracted from each member's tracks; recommended on real data, where centromere flanks and segmental duplications "
                                   "are otherwise called as events. A file, or '1kg-dragen' (twelve public 1000 Genomes genomes called by DRAGEN 3.7.6, 1-Mb bins)")
    r.add_argument("--events", help="events from another method to draw and match (sample chrom start end label [type]; or NGS-DOSE's karyotype/events.tsv)")
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
    m.set_defaults(fn=cmd_mock)
    c = sub.add_parser("cohort", help="cohort report over many runs: counts, landscape figure, every event, per-trio metrics, concordance, guide")
    c.add_argument("--runs", nargs="+", required=True, help="run directories (globs)")
    c.add_argument("--out", required=True)
    c.add_argument("--events", help="events from another method to match (as for run)")
    c.add_argument("--genome", default="grch38")
    c.set_defaults(fn=cmd_cohort)
    cb = sub.add_parser("calibrate", help="detection and cell-fraction accuracy on the simulator over a grid of cell fractions, depths and event sizes")
    cb.add_argument("--out", required=True)
    cb.add_argument("--cell-fractions", default="0.05,0.1,0.2,0.3,0.5,1", help="comma-separated cell fractions (default 0.05,0.1,0.2,0.3,0.5,1)")
    cb.add_argument("--depths", default="30", help="comma-separated mean depths (default 30)")
    cb.add_argument("--sizes", default="5,10,20,50", help="comma-separated event sizes in Mb (default 5,10,20,50)")
    cb.add_argument("--replicates", type=int, default=1)
    cb.add_argument("--seed", type=int, default=1)
    cb.add_argument("--sites-per-mb", type=int, default=60, help="simulated site density (real WGS: about 1000 per Mb; 60 is quick but understates the phased scan)")
    cb.add_argument("--contigs", default="", help="restrict the simulation to these chromosomes, e.g. chr1,...,chr16 with a dense --sites-per-mb")
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
