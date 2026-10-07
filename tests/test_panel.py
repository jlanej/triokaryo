"""A reference panel from other genomes takes the shared depth structure out and leaves the planted events."""
import os

import numpy as np

from triokaryo.cli import main
from triokaryo.mock import write_mock
from triokaryo.panel import load_panel
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio


def test_panel_from_vcfs_and_runs_keeps_the_events(mock, tmp_path):
    # the panel: a null trio (three genomes) read from its VCF, and the same three again from a run's bins.tsv
    m0 = write_mock(str(tmp_path / "null"), seed=6, no_events=True)
    trio0 = read_trios(m0["trios"])[0]
    run_trio(m0["vcf"], trio0, str(tmp_path / "null_run"), figures=False, log=lambda s: None)
    panel = tmp_path / "panel.tsv"
    assert main(["panel", "--vcfs", m0["vcf"], "--runs", str(tmp_path / "null_run"), "--out", str(panel), "--thin", "2"]) == 0
    p = load_panel(str(panel))
    assert len(p) > 2500 and all(e[0] == 6 for e in list(p.values())[:50])
    # applied to the trio with events: every planted event still found, the X copies still read
    trio = read_trios(mock["trios"])[0]
    res = run_trio(mock["vcf"], trio, str(tmp_path / "with_panel"), gc_track=mock["gc"], panel=str(panel), figures=False, log=lambda s: None)
    found = {(e.sample, e.chrom, e.type) for e in res["events"]}
    for t in mock["truth_data"]["events"]:
        assert (t["sample"], t["chrom"], t["type"]) in found, t["label"]
    assert res["x_copies"] == {"child": 1, "father": 1, "mother": 2}
    assert res["summary"]["panel"] and res["summary"]["bins_masked"] < 100
    # bins.tsv carries the panel's columns; the summary knows a panel was used
    head = open(os.path.join(str(tmp_path / "with_panel"), "bins.tsv")).readline().split("\t")
    assert "panel_median" in head and "masked" in head


def test_shipped_panel_resolves_by_name(mock, tmp_path):
    """--panel 1kg-dragen: the twelve-genome panel shipped as package data (built from public 1000 Genomes DRAGEN calls)."""
    out = tmp_path / "shipped"
    rc = main(["run", "--vcf", mock["vcf"], "--pedigree", mock["trios"], "--child", "KID", "--panel", "1kg-dragen", "--out", str(out), "--no-figures", "--thin", "3"])
    assert rc == 0
    head = open(os.path.join(str(out), "bins.tsv")).readline().split("\t")
    rows = open(os.path.join(str(out), "bins.tsv")).read().splitlines()[1:]
    assert "panel_median" in head and any(r.split("\t")[head.index("panel_median")] != "NA" for r in rows)


def test_panel_masks_bins_where_its_genomes_disagree(mock, tmp_path):
    """A bin whose LRR spreads across the panel's genomes by more than the smallest reportable event (robust SD above 0.10) is masked;
    one within it is not; and the shipped GC track, the default at 1 Mb on GRCh38, is loaded unless 'none' is asked for."""
    from triokaryo.genome import genome
    from triokaryo.model import PANEL_MAX_RSD, load_gc_track
    from triokaryo.cli import resolve_gc_track
    G = genome()
    rows = ["chrom\tstart\tend\tn\tlrr_median\tlrr_rsd\tn_bdev\tbdev_median\tn_het\thet_rate_median"]
    loud = {("chr1", 10_000_000): 0.15, ("chr1", 11_000_000): 0.09, ("chr2", 50_000_000): PANEL_MAX_RSD + 0.001}
    for c in G.chroms:
        for st in range(0, G.length[c], 1_000_000):
            n = 4 if c == "chrY" else 6
            rows.append("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%d\t0.0700\t%d\t0.4000" % (c, st, min(st + 1_000_000, G.length[c]), n, 0.0, loud.get((c, st), 0.03), n, n))
    panel = tmp_path / "panel.tsv"
    panel.write_text("\n".join(rows) + "\n")
    trio = read_trios(mock["trios"])[0]
    res = run_trio(mock["vcf"], trio, str(tmp_path / "masked"), gc_track=resolve_gc_track(None, "grch38", 1_000_000), panel=str(panel), thin=3, figures=False, log=lambda s: None)
    b = res["bins"]
    at = lambda c, st: int(np.flatnonzero((np.asarray(b.chrom) == c) & (np.asarray(b.start) == st))[0])
    assert b.masked[at("chr1", 10_000_000)] and b.masked[at("chr2", 50_000_000)] and not b.masked[at("chr1", 11_000_000)] and not b.masked[at("chr1", 12_000_000)]
    assert res["summary"]["gc_corrected"] and np.isfinite(b.gc).sum() > 2900
    gc = load_gc_track(resolve_gc_track("hg38", "grch38", 1_000_000))
    assert len(gc) == 2970 and abs(gc[("chr1", 1_000_000)] - 0.573) < 0.01        # 3102 bins, 132 of them all N
    assert np.mean([v for (c, _), v in gc.items() if c == "chr19"]) > np.mean([v for (c, _), v in gc.items() if c == "chr4"]) + 0.08
    assert resolve_gc_track(None, "grch38", 500_000) is None and resolve_gc_track("none", "grch38", 1_000_000) is None


def test_cohort_panel_corrects_the_y_for_a_father_of_a_daughter(tmp_path):
    """Mosaic loss of Y in a father without a son needs the Y level of other males. A panel built from the cohort's runs (its
    fathers and sons, --roles father,child) holds the male Y and X levels with their mappability deficit; a 15% loss, invisible
    without a panel (no son to compare, and under the 25% floor), is then read from the panel-corrected level, and the daughter's
    X is not read as a loss despite the male-only panel."""
    from triokaryo.mock import DAD
    deficit = (0.93, 0.95)
    runs = []
    for i in range(3):                                                    # three ordinary trios: six males for the Y rows
        m = write_mock(str(tmp_path / ("t%d" % i)), seed=20 + i, no_events=True, prefix="T%d_" % i, sex_deficit=deficit)
        out = str(tmp_path / ("run%d" % i))
        run_trio(m["vcf"], read_trios(m["trios"])[0], out, figures=False, log=lambda s: None)
        runs.append(out)
    panel = tmp_path / "cohort.panel.tsv"
    assert main(["panel", "--runs"] + runs + ["--roles", "father,child", "--out", str(panel)]) == 0
    p = load_panel(str(panel))
    y = [v for (c, _), v in p.items() if c == "chrY"]
    x = [v for (c, _), v in p.items() if c == "chrX"]
    assert len(y) >= 10 and all(v[0] == 6 for v in y) and abs(np.median([v[1] for v in y]) - np.log2(0.95)) < 0.05         # one Y copy with its deficit, on the diploid scale
    assert len(x) >= 100 and abs(np.median([v[1] for v in x]) - np.log2(0.93)) < 0.05                                       # the males' X likewise
    loy = [dict(member=DAD, chrom="chrY", start=0, end=None, f=0.15, delta={"h0": -1}, label="mosaic loss of Y in the father (15% of cells)", type="loss", origin="", inherited=False)]
    m = write_mock(str(tmp_path / "daughter"), seed=30, child_sex="F", events=loy, sex_deficit=deficit)
    trio = read_trios(m["trios"])[0]
    res0 = run_trio(m["vcf"], trio, str(tmp_path / "d0"), figures=False, log=lambda s: None)
    assert not [e for e in res0["events"] if e.chrom == "chrY"] and abs(res0["summary"]["father_y_copies_raw"] - 0.81) < 0.05   # deficit and loss confounded: under the 25% floor
    res1 = run_trio(m["vcf"], trio, str(tmp_path / "d1"), panel=str(panel), figures=False, log=lambda s: None)
    ys = [e for e in res1["events"] if e.chrom == "chrY"]
    assert len(ys) == 1 and ys[0].sample == "DAD" and ys[0].type == "loss" and abs(ys[0].f - 0.15) < 0.06 and "against the panel" in ys[0].note, [(e.sample, e.type, e.f) for e in ys]
    assert res1["summary"]["y_panel"] and not [e for e in res1["events"] if e.chrom == "chrX"], [(e.sample, e.chrom, e.type, e.f) for e in res1["events"]]
    assert abs(res1["summary"]["child_x_copies_raw"] - 2) < 0.06 and abs(res1["summary"]["mother_x_copies_raw"] - 2) < 0.06
    assert res1["summary"]["father_karyotype"] == "mos 45,X[%.2f]/46,XY" % ys[0].f
