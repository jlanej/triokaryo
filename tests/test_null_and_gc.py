"""No event where none was planted; the GC correction removes the GC bias; a 47,XXY child reads two X copies."""
import os

import numpy as np

from triokaryo.mock import write_mock
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio


def test_null_trio_calls_nothing(tmp_path):
    m = write_mock(str(tmp_path / "null"), seed=3, no_events=True)
    trio = read_trios(m["trios"])[0]
    res = run_trio(m["vcf"], trio, str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    assert res["events"] == [], [(e.sample, e.chrom, e.type, e.span, e.f) for e in res["events"]]
    assert res["x_copies"] == {"child": 1, "father": 1, "mother": 2}


def test_gc_correction_tightens_the_lrr(tmp_path):
    m = write_mock(str(tmp_path / "gc"), seed=4, no_events=True)
    trio = read_trios(m["trios"])[0]
    raw = run_trio(m["vcf"], trio, str(tmp_path / "raw"), gc_track=None, figures=False, log=lambda s: None)["bins"]
    cor = run_trio(m["vcf"], trio, str(tmp_path / "cor"), gc_track=m["gc"], figures=False, log=lambda s: None)["bins"]
    for mem in range(3):                   # the correction never adds spread, and takes the GC part out (the rest is counting noise)
        a = raw.lrr[mem][raw.autosomal]
        b = cor.lrr_gc[mem][cor.autosomal]
        assert np.nanstd(b) < np.nanstd(a)
    # the raw LRR follows GC in the mock (a bias per sample); the corrected one does not
    g = cor.gc[cor.autosomal]                  # the raw run loaded no track: take the bins' GC from the corrected run (the same bins)
    ok = np.isfinite(g) & np.isfinite(raw.lrr[0][raw.autosomal])
    r_raw = np.corrcoef(g[ok], raw.lrr[0][raw.autosomal][ok])[0, 1]
    r_cor = np.corrcoef(g[ok], cor.lrr_gc[0][cor.autosomal][ok])[0, 1]
    assert abs(r_raw) > 0.5 and abs(r_cor) < 0.2


def test_xxy_child_reads_two_x_copies(tmp_path):
    m = write_mock(str(tmp_path / "xxy"), seed=5, no_events=True, xxy=True)
    trio = read_trios(m["trios"])[0]
    res = run_trio(m["vcf"], trio, str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    assert res["x_copies"]["child"] == 2 and res["summary"]["child_x_check"].startswith("X copies 2 in a reported male")
    assert not [e for e in res["events"] if e.chrom == "chrX"]
