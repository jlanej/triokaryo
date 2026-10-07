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


def test_karyotype_string_grammar():
    """The ISCN-like string: modal number with the complement, constitutional terms in chromosome order, each mosaic event its own
    line with its cell fraction, the parent of origin and the meiotic stage abbreviated."""
    from triokaryo.genome import genome
    from triokaryo.segment import Event
    from triokaryo.sexchrom import karyotype_string
    G = genome()
    order = {c: i for i, c in enumerate(G.chroms)}
    L = G.length
    ev = [Event("K", "child", "chr21", 0, L["chr21"], "whole", "gain", f_lrr=1.0, origin_phase="extra copy maternal", stage="meiosis I"),
          Event("K", "child", "chr12", 0, L["chr12"], "whole", "gain", f_lrr=0.3, origin_phase="extra copy paternal", stage="mitotic, or meiosis II without a crossover"),
          Event("K", "child", "chr7", 0, L["chr7"], "whole", "LOH", f_baf=1.0, origin="maternal copy retained (paternal replaced)", note="Mendelian errors at the informative sites: a uniparental isodisomy"),
          Event("K", "child", "chr15", 0, L["chr15"], "whole", "UPD", f_baf=1.0, origin_phase="both copies maternal (heterodisomy)"),
          Event("K", "child", "chr18", 55_000_000, L["chr18"], "stretch", "loss", f_lrr=1.0, origin="paternal copy lost"),
          Event("K", "child", "chr6", 0, 59_800_000, "p", "LOH", f_baf=0.4, origin="maternal copy retained (paternal replaced)"),
          Event("K", "child", "chr1", 0, 5_000_000, "stretch", "gain", f_lrr=0.12, note="the phased bands read a share of 4% against the depth's 12%: the depth's call may be an artefact"),
          Event("K", "child", "chrX", 0, L["chrX"], "whole", "gain", f_lrr=1.0, origin_phase="extra copy maternal", stage="meiosis I")]
    s = karyotype_string(ev, "XXY", order)
    assert s == ("mos 48,XXY,dup(1)(0.0-5.0Mb)?[0.12]/48,XXY,loh(6)(p:0.0-59.8Mb)mat[0.40]/49,XXY,+12pat(MII/mit)[0.30]"
                 "/48,XXY(mat,MI),upd(7)mat(iso),upd(15)mat(hetero),del(18)(55.0-80.4Mb)pat,+21mat(MI)"), s
    assert karyotype_string([], "XX", order) == "46,XX" and karyotype_string([], "", order) == "46,?"
    assert karyotype_string([Event("M", "mother", "chrX", 0, L["chrX"], "whole", "loss", f_lrr=0.4)], "XX", order) == "mos 45,X[0.40]/46,XX"
    assert karyotype_string([Event("D", "father", "chrY", 0, L["chrY"], "whole", "loss", f_lrr=0.3)], "XY", order) == "mos 45,X[0.30]/46,XY"
    assert karyotype_string([Event("K", "child", "chrX", 0, L["chrX"], "whole", "gain", f_lrr=0.4, origin_phase="extra copy paternal", stage="meiosis I (paternal: X and Y transmitted together)")],
                            "XY", order) == "mos 47,XXY(pat,MI)[0.40]/46,XY"


def test_daughter_xxx_and_turner(tmp_path):
    """A 47,XXX from a maternal meiosis I error (both maternal homologues: heterodisomic at the centromere), a paternal 47,XXX (the
    father's single X twice: isodisomic, staged meiosis II or post-zygotic) and a 45,X with the paternal X lost, each with the
    complement, its check against the pedigree sex and the karyotype string; no other event called."""
    from triokaryo.mock import DAUGHTER_TURNER, DAUGHTER_XXX_PATERNAL
    for name, kw, expect in (("xxx_mat", dict(no_events=True, xxx=True), ("gain", "XXX", "extra copy maternal", "meiosis I", "heterodisomic", "47,XXX(mat,MI)")),
                             ("xxx_pat", dict(events=DAUGHTER_XXX_PATERNAL), ("gain", "XXX", "extra copy paternal", "meiosis II, or post-zygotic", "isodisomic", "47,XXX(pat,MII)")),
                             ("turner", dict(events=DAUGHTER_TURNER), ("loss", "X", "paternal copy lost", "", "", "45,X(pat)"))):
        m = write_mock(str(tmp_path / name), seed=7, child_sex="F", **kw)
        trio = read_trios(m["trios"])[0]
        assert trio.kid_sex == "F"
        res = run_trio(m["vcf"], trio, str(tmp_path / (name + "_out")), gc_track=m["gc"], figures=False, log=lambda s: None)
        kind, comp, origin, stage, centro, kar = expect
        (e,) = res["events"]
        s = res["summary"]
        assert e.sample == "KID" and e.chrom == "chrX" and e.type == kind and e.span == "whole" and abs(e.f - 1.0) < 0.1, (name, e.type, e.f)
        assert e.origin == origin and e.origin_phase == origin and e.stage.startswith(stage) and e.centromere == centro, (name, e.origin, e.origin_phase, e.stage, e.centromere)
        assert s["child_sex_karyotype"] == comp and s["child_sex_check"].startswith("%s in a reported female" % comp) and s["child_karyotype"] == kar, (name, s["child_sex_check"], s["child_karyotype"])
        assert s["child_y_copies"] == 0 and abs(e.f_phase - 1.0) < 0.1


def test_whole_x_event_without_pedigree_sex(tmp_path):
    """With no pedigree sex, the Y implies the sex: a 47,XXY son is still a whole-X gain against one expected copy."""
    from triokaryo.pedigree import trio_from_args
    m = write_mock(str(tmp_path / "xxy"), seed=5, no_events=True, xxy=True)
    res = run_trio(m["vcf"], trio_from_args("KID", "DAD", "MOM"), str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    (e,) = res["events"]
    assert e.chrom == "chrX" and e.type == "gain" and abs(e.f - 1.0) < 0.1 and "implied by the Y" in e.note and e.origin == "extra copy maternal" and e.stage == "meiosis I"
    assert res["summary"]["child_sex_check"] == "" and res["summary"]["child_sex_karyotype"] == "XXY" and res["summary"]["child_karyotype"] == "47,XXY(mat,MI)"


def test_x_level_corrected_within_the_trio_without_a_panel(tmp_path):
    """Without a panel, a mappability deficit on the X (7% here) would read as a mosaic X loss in every female; the median deviation
    of the members' X levels from their pedigree expectation corrects it, the Y (not corrected) still rounds to one copy, and a
    47,XXY under the same deficit is still a whole-X gain of one copy."""
    m = write_mock(str(tmp_path / "deficit"), seed=12, no_events=True, sex_deficit=(0.93, 0.90))
    trio = read_trios(m["trios"])[0]
    res = run_trio(m["vcf"], trio, str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    s = res["summary"]
    assert res["events"] == [], [(e.sample, e.chrom, e.type, e.f) for e in res["events"]]
    assert abs(s["x_offset_trio"] - np.log2(0.93)) < 0.03, s["x_offset_trio"]
    assert abs(s["mother_x_copies_raw"] - 2) < 0.06 and abs(s["child_x_copies_raw"] - 1) < 0.03 and abs(s["father_x_copies_raw"] - 1) < 0.03
    assert abs(s["father_y_copies_raw"] - 0.90) < 0.05 and s["father_y_copies"] == 1 and abs(s["y_father_son_log2"]) < 0.1
    m2 = write_mock(str(tmp_path / "xxy_deficit"), seed=12, no_events=True, xxy=True, sex_deficit=(0.93, 0.90))
    res2 = run_trio(m2["vcf"], read_trios(m2["trios"])[0], str(tmp_path / "out2"), gc_track=m2["gc"], figures=False, log=lambda s: None)
    (e,) = res2["events"]
    assert e.chrom == "chrX" and e.type == "gain" and abs(e.f - 1.0) < 0.08 and res2["summary"]["child_karyotype"] == "47,XXY(mat,MI)", (e.f, res2["summary"]["child_karyotype"])


def test_mosaic_maternal_xxy_is_staged_from_a_low_cell_fraction(tmp_path):
    """A 46,XY/47,XXY mosaic with both maternal X homologues (a maternal meiosis I error, the extra X lost in most cells) in 20% of
    cells: the whole-X gain is maternal, and the stage reads meiosis I from the auxiliary track at the father's sites, which does
    not depend on the child being called heterozygous."""
    from triokaryo.mock import KID
    ev = [dict(member=KID, chrom="chrX", start=0, end=None, f=0.2, delta={"mat_other": +1}, label="mosaic maternal XXY, 20%", type="gain", origin="extra copy maternal", inherited=False)]
    m = write_mock(str(tmp_path / "xxy20"), seed=13, events=ev)
    res = run_trio(m["vcf"], read_trios(m["trios"])[0], str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    (e,) = [e for e in res["events"] if e.chrom == "chrX"]
    assert e.type == "gain" and abs(e.f - 0.2) < 0.06 and e.origin == "extra copy maternal" and e.origin_phase == "extra copy maternal", (e.f, e.origin, e.origin_phase)
    assert e.stage == "meiosis I" and e.centromere == "heterodisomic" and e.n_crossovers == 0, (e.stage, e.centromere, e.n_crossovers, e.hetero_share)
    assert res["summary"]["child_karyotype"] == "mos 47,XXY(mat,MI)[%.2f]/46,XY" % e.f
