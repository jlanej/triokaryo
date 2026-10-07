"""triokaryo: large chromosomal events in a trio, from its VCF.

  triokaryo run --vcf trio.vcf.gz --pedigree trios.tsv --child KID --out out/KID [--gc-track gc.tsv] [--events other_calls.tsv]
  triokaryo run --vcf trio.vcf.gz --child KID --father DAD --mother MOM --sex M,M,F --out out/KID
  triokaryo panel --vcfs a.vcf.gz b.vcf.gz ... --out panel.tsv        # or --runs 'out/*' from earlier runs
  triokaryo run ... --panel panel.tsv
  triokaryo gc-track --fasta ref.fa --out gc.tsv [--bin 1000000]
  triokaryo mock --out mock_dir [--seed 1] [--no-events] [--xxy] [--contigs chr15,chr16,chr17 --sites-per-mb 1000 --low-share]
  triokaryo cohort --runs 'out/*' --out cohort [--events other_calls.tsv]   # the cohort report + guide
  triokaryo report --runs 'out/*'                   # a run's page again from its tables (after a change to the page)
  triokaryo guide --out guide.html                  # how to read every figure, colour and column
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
    rows, samples = build_panel(genome(a.genome), a.bin, vcfs=a.vcfs or (), samples=set(a.samples.split(",")) if a.samples else None, runs=runs, log=_log, thin=a.thin)
    if not rows:
        sys.exit("no genome read for the panel (--vcfs and/or --runs)")
    write_panel(a.out, rows, samples, a.bin)
    import numpy as np
    from .model import PANEL_MAX_RSD, PANEL_MIN_N
    unp = sum(1 for r in rows if r["n"] < PANEL_MIN_N or (np.isfinite(r["lrr_rsd"]) and r["lrr_rsd"] > PANEL_MAX_RSD))
    if len(samples) < PANEL_MIN_N:
        _log("WARNING: %d genomes: a panel needs %d or more, or every bin is left out of the calls" % (len(samples), PANEL_MIN_N))
    _log("panel of %d genomes, %d bins (%d unpinned) -> %s" % (len(samples), len(rows), unp, a.out))
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
    paths = write_mock(a.out, seed=a.seed, sites_per_mb=a.sites_per_mb, no_events=a.no_events, xxy=a.xxy, prefix=a.prefix, events=events,
                       contigs=a.contigs.split(",") if a.contigs else None)
    _log("mock trio -> %s" % paths["vcf"])
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
    r = sub.add_parser("run", help="one trio: scan, bin, call, read, figures, page")
    r.add_argument("--vcf", required=True)
    r.add_argument("--pedigree", help="a trios file (#kid dad mom [kid_sex dad_sex mom_sex]) or a PED")
    r.add_argument("--child")
    r.add_argument("--father")
    r.add_argument("--mother")
    r.add_argument("--sex", default="", help="child,father,mother as M/F when given outright, e.g. M,M,F")
    r.add_argument("--out", required=True)
    r.add_argument("--gc-track", help="bin GC (triokaryo gc-track): the LRR is corrected with it")
    r.add_argument("--panel", help="a reference panel (triokaryo panel): the depth, band and heterozygosity structure every genome shares, taken "
                                   "out of each member's tracks; needed on real data, where centromere flanks and segmental duplications read as "
                                   "events otherwise. A file, or '1kg-dragen' (twelve public 1000 Genomes genomes called by DRAGEN 3.7.6, 1-Mb bins)")
    r.add_argument("--events", help="events from elsewhere to draw and compare (sample chrom start end label; or NGS-DOSE's karyotype/events.tsv)")
    r.add_argument("--bin", type=int, default=1_000_000)
    r.add_argument("--min-dp", type=int, default=8)
    r.add_argument("--min-gq", type=int, default=20)
    r.add_argument("--thin", type=int, default=1, help="keep every n-th usable site (speed; 1 = all)")
    r.add_argument("--genome", default="grch38")
    r.add_argument("--min-abs", type=float, default=0.07, help="|LRR| a gain or loss needs (0.07: a share of cells of about 10%%)")
    r.add_argument("--z", type=float, default=5.0, help="the split statistic a segment boundary needs")
    r.add_argument("--min-len", type=int, default=5, help="bins per segment")
    r.add_argument("--min-f", type=float, default=0.10, help="the smallest share of cells reported")
    r.add_argument("--no-figures", action="store_true")
    r.set_defaults(fn=cmd_run)
    pn = sub.add_parser("panel", help="a reference panel: the median LRR per bin over other genomes (per-sample or multi-sample VCFs, or earlier runs)")
    pn.add_argument("--vcfs", nargs="*", help="VCFs to read (every sample of each unless --samples)")
    pn.add_argument("--samples", help="comma-separated sample names to take from the VCFs")
    pn.add_argument("--runs", nargs="*", help="earlier run directories (bins.tsv + summary.tsv), e.g. 'out/*'")
    pn.add_argument("--out", required=True)
    pn.add_argument("--bin", type=int, default=1_000_000)
    pn.add_argument("--thin", type=int, default=1)
    pn.add_argument("--genome", default="grch38")
    pn.set_defaults(fn=cmd_panel)
    g = sub.add_parser("gc-track", help="bin GC from a reference FASTA")
    g.add_argument("--fasta", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("--bin", type=int, default=1_000_000)
    g.add_argument("--genome", default="grch38")
    g.set_defaults(fn=cmd_gc_track)
    m = sub.add_parser("mock", help="a mock trio VCF with planted events (no real data)")
    m.add_argument("--out", required=True)
    m.add_argument("--seed", type=int, default=1)
    m.add_argument("--sites-per-mb", type=int, default=60)
    m.add_argument("--no-events", action="store_true")
    m.add_argument("--xxy", action="store_true", help="the child 47,XXY (two X copies, a male)")
    m.add_argument("--prefix", default="", help="a tag before the sample names KID, DAD, MOM (several mock trios in one cohort)")
    m.add_argument("--contigs", default="", help="only these chromosomes, comma-separated (a dense small mock)")
    m.add_argument("--low-share", action="store_true", help="plant the low-share events (under the depth's threshold) instead of the default ones")
    m.set_defaults(fn=cmd_mock)
    c = sub.add_parser("cohort", help="the cohort report over many trios' runs: counts, the landscape figure, every event, the trios, the concordance, the guide")
    c.add_argument("--runs", nargs="+", required=True, help="run directories (globs)")
    c.add_argument("--out", required=True)
    c.add_argument("--events", help="events from elsewhere to match (as for run)")
    c.add_argument("--genome", default="grch38")
    c.set_defaults(fn=cmd_cohort)
    g = sub.add_parser("guide", help="the guide: how to read every figure, colour, call and column (one self-contained page with pattern cards)")
    g.add_argument("--out", required=True, help="the HTML to write")
    g.add_argument("--figures", help="keep the pattern figures (PNG, SVG, PDF, sidecars, legends) in this directory")
    g.set_defaults(fn=cmd_guide)
    rp = sub.add_parser("report", help="a run's page, sidecars, legends and guide again from its tables and figures (no VCF needed)")
    rp.add_argument("--runs", nargs="+", required=True, help="run directories (globs)")
    rp.set_defaults(fn=cmd_report)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
