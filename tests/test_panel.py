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
