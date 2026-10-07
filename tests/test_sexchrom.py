"""Sex chromosomes: the X and Y copy numbers and the complement against the pedigree sex; whole-X and whole-Y mosaics with their
cell fraction and, on the X, the parent of origin with a hemizygous baseline; the panel's Y convention."""
import json
import os

import numpy as np

from triokaryo.mock import SEX_EVENTS, write_mock
from triokaryo.panel import load_panel
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio
from triokaryo.sexchrom import karyotype, sex_check


def test_karyotype_and_check_strings():
    assert karyotype(1, 1) == "XY" and karyotype(2, 0) == "XX" and karyotype(2, 1) == "XXY" and karyotype(1, 0) == "X" and karyotype(1, 2) == "XYY"
    assert karyotype(float("nan"), 1) == "" and karyotype(2, float("nan")) == ""
    st = dict(x_copies=2, y_copies=1, x_copies_raw=2.0, y_copies_raw=1.0)
    assert sex_check("M", st).startswith("XXY in a reported male") and sex_check("F", dict(x_copies=1, y_copies=0)).startswith("X in a reported female (45,X)")
    assert sex_check("F", dict(x_copies=2, y_copies=0)) == "agrees" and sex_check("", st) == "" and sex_check("M", dict(x_copies=1, y_copies=float("nan"))) == "agrees (Y not in the VCF)"


def test_sex_chromosome_mosaics(tmp_path):
    """A 30% loss of Y in the father (from the father/son Y ratio: no panel), a 40% 45,X/46,XX mosaic in the mother with the lost X
    identified as the transmitted one or not, and a 40% 46,XY/47,XXY mosaic in the son with the extra X paternal, from the
    informative sites with a hemizygous baseline and from the phased track; nothing else called."""
    m = write_mock(str(tmp_path / "sex"), seed=5, events=SEX_EVENTS)
    trio = read_trios(m["trios"])[0]
    res = run_trio(m["vcf"], trio, str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    ev = {(e.sample, e.chrom, e.type): e for e in res["events"]}
    assert sorted(ev) == [("DAD", "chrY", "loss"), ("KID", "chrX", "gain"), ("MOM", "chrX", "loss")], [(e.sample, e.chrom, e.type, e.span, e.f) for e in res["events"]]
    assert all(e.span == "whole" for e in res["events"])
    kid, dad, mom = ev[("KID", "chrX", "gain")], ev[("DAD", "chrY", "loss")], ev[("MOM", "chrX", "loss")]
    assert abs(kid.f - 0.40) < 0.08 and kid.origin == "extra copy paternal" and kid.origin_phase == "extra copy paternal" and abs(kid.f_phase - 0.40) < 0.1
    assert kid.stage.startswith("meiosis I (paternal") and kid.mie_rate > 0.1 and "46,XY/47,XXY" in kid.note
    assert abs(dad.f - 0.30) < 0.08 and "father/son" in dad.note and not np.isfinite(dad.phase_shift)
    assert abs(mom.f - 0.40) < 0.08 and abs(mom.f_phase - 0.40) < 0.1 and "45,X/46,XX" in mom.note
    transmitted_lost = json.load(open(m["truth"]))["transmitted_maternal"]["chrX"] == 0
    assert mom.origin_phase == ("the lost homologue is the one passed to the child" if transmitted_lost else "the lost homologue is the one not passed to the child")
    s = res["summary"]
    assert abs(s["child_x_copies_raw"] - 1.4) < 0.08 and abs(s["mother_x_copies_raw"] - 1.6) < 0.08 and abs(s["father_y_copies_raw"] - 0.7) < 0.08
    assert abs(s["y_father_son_log2"] - np.log2(0.7)) < 0.15
    assert [s["%s_sex_karyotype" % r] for r in ("child", "father", "mother")] == ["XY", "XY", "XX"]


def test_panel_y_rows_are_read_on_the_diploid_scale(tmp_path):
    """A panel written with the Y at the male's one-copy level is lifted by one unit on loading; one already on the diploid scale is not."""
    for level, expect in ((-1.1, -0.07), (-0.05, -0.02)):                      # the row checked is the first Y row, level + 0.03
        p = tmp_path / ("panel_%s.tsv" % str(level).replace("-", "m"))
        rows = ["chrom\tstart\tend\tn\tlrr_median\tlrr_rsd\tn_bdev\tbdev_median\tn_het\thet_rate_median"]
        rows += ["chr1\t%d\t%d\t12\t0.0100\t0.0500\t12\t0.0700\t12\t0.4000" % (i * 10 ** 6, (i + 1) * 10 ** 6) for i in range(5)]
        rows += ["chrY\t%d\t%d\t4\t%.4f\t0.0500\t0\tNA\t4\t0.0000" % (i * 10 ** 6, (i + 1) * 10 ** 6, level + 0.01 * i) for i in range(3, 8)]
        p.write_text("\n".join(rows) + "\n")
        pan = load_panel(str(p))
        assert abs(pan[("chr1", 0)][1] - 0.01) < 1e-9 and abs(pan[("chrY", 3_000_000)][1] - expect) < 1e-9, (level, pan[("chrY", 3_000_000)])
