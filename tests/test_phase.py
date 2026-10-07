"""The phased tracks: the sign of the shift names the parent, the auxiliary tracks tell one homologue from two, the
heterodisomy is found, the edges come at site resolution, the step fit is exact, and the dense mock's low-share events
come out of the phased scan."""
import os

import numpy as np

from triokaryo.mock import LOW_SHARE_EVENTS, write_mock
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio
from triokaryo.smooth import folded_null, tv_denoise, windows_by_count


def _tv_dual(y, lam, iters=30000):
    n = len(y)
    u = np.zeros(n - 1)
    def x_of(u):                                                # x = y - D^T u, with (D x)_i = x_{i+1} - x_i
        x = y.copy()
        x[:-1] += u
        x[1:] -= u
        return x
    for _ in range(iters):
        x = x_of(u)
        u = np.clip(u + 0.25 * (x[1:] - x[:-1]), -lam, lam)    # a projected gradient step on the dual
    return x_of(u)


def test_step_fit_is_the_exact_total_variation_solution():
    rng = np.random.default_rng(0)
    for _ in range(5):
        y = rng.normal(0, 0.3, 40) + np.repeat(rng.normal(0, 1, 4), 10)
        lam = float(rng.uniform(0.1, 1.0))
        assert np.max(np.abs(tv_denoise(y, lam) - _tv_dual(y, lam))) < 1e-6
    assert np.allclose(tv_denoise(np.array([1.0, 1.0, 5.0, 1.0, 1.0]), 3.0), 1.8)      # a lone spike of 4 under 2 lambda: flattened to the mean
    assert abs(folded_null(np.array([30]))[0] - 0.0722) < 1e-3
    assert windows_by_count(np.array([1, 2, 3, 4, 5, 6, 7, 10_000_000, 10_000_001, 10_000_002]), 3) == [(0, 3), (3, 7), (7, 10)]


def _by(run, sample, chrom, etype):
    (e,) = [e for e in run["events"] if e.sample == sample and e.chrom == chrom and e.type == etype]
    return e


def test_phased_shift_names_the_parent_and_the_share(run, mock):
    truth = mock["truth_data"]
    e21 = _by(run, "KID", "chr21", "gain")                      # a maternal meiotic trisomy: the maternal fraction 2/3 at every site
    assert abs(e21.phase_shift - 1 / 6) < 0.03 and e21.origin_phase == "extra copy maternal" and abs(e21.f_phase - 1.0) < 0.15
    e12 = _by(run, "KID", "chr12", "gain")                      # a 30% paternal mitotic gain: 1/2 - 0.065
    assert abs(e12.phase_shift + 0.065) < 0.02 and e12.origin_phase == "extra copy paternal" and abs(e12.f_phase - 0.30) < 0.08
    e18 = _by(run, "KID", "chr18", "loss")                      # the paternal copy lost: the maternal fraction 1
    assert abs(e18.phase_shift - 0.5) < 0.02 and e18.origin_phase == "paternal copy lost"
    e7 = _by(run, "KID", "chr7", "LOH")                         # a maternal isodisomy: 1
    assert abs(e7.phase_shift - 0.5) < 0.02 and e7.origin_phase.startswith("maternal copy retained")
    e6 = _by(run, "KID", "chr6", "LOH")                         # a 40% copy-neutral LOH, the maternal copy retained: 1/2 + 0.2
    assert abs(e6.phase_shift - 0.2) < 0.03 and abs(e6.f_phase - 0.40) < 0.08
    d10 = _by(run, "DAD", "chr10", "gain")                      # the father's gain lies on the homologue he passed on (the mock fixes it so)
    assert abs(d10.phase_shift - 1 / 6) < 0.03 and d10.origin_phase == "the duplicated homologue is the one passed to the child"
    m8 = _by(run, "MOM", "chr8", "gain")                        # the mother's 15% gain on homologue h0: the sign says whether h0 went to the child
    expect = 1 if truth["transmitted_maternal"]["chr8"] == 0 else -1
    assert np.sign(m8.phase_shift) == expect and abs(abs(m8.phase_shift) - 0.035) < 0.015
    # the opposite-homozygote reading and the phased one agree everywhere
    for e in run["events"]:
        if e.role == "child" and e.origin and e.origin_phase and e.type != "UPD":
            assert e.origin == e.origin_phase, (e.chrom, e.origin, e.origin_phase)


def test_one_homologue_or_two(run):
    assert _by(run, "KID", "chr21", "gain").homologues.startswith("the two maternal copies are different homologues")
    assert _by(run, "KID", "chr12", "gain").homologues.startswith("the two paternal copies are one homologue")
    assert _by(run, "KID", "chr10", "gain").homologues.startswith("the two paternal copies are one homologue")
    assert _by(run, "KID", "chr15", "UPD").homologues.startswith("the two maternal copies are different homologues")


def test_heterodisomy_is_found_by_the_phased_scan_alone(run):
    e = _by(run, "KID", "chr15", "UPD")
    assert e.source == "phased" and e.span == "whole" and e.origin_phase == "both copies maternal (heterodisomy)"
    # Mendel breaks only where the parents are opposite homozygotes (a twentieth of the sites), unlike an isodisomy's hom child everywhere
    assert abs(e.phase_shift - 0.5) < 0.02 and e.het_rate_rel > 0.5 and abs(e.lrr) < 0.05 and 0.02 < e.mie_rate < 0.12
    assert "heterodisomy" in e.note


def test_edges_at_site_resolution(run):
    e18 = _by(run, "KID", "chr18", "loss")
    assert np.isfinite(e18.start_fine) and abs(e18.start_fine - 55_000_000) < 600_000 and not np.isfinite(e18.end_fine)
    e2 = _by(run, "KID", "chr2", "loss")
    assert np.isfinite(e2.start_fine) and np.isfinite(e2.end_fine) and abs(e2.start_fine - 100_000_000) < 1_200_000 and abs(e2.end_fine - 106_000_000) < 1_200_000
    e21 = _by(run, "KID", "chr21", "gain")
    assert not np.isfinite(e21.start_fine) and not np.isfinite(e21.end_fine)            # a whole chromosome has no edge to refine


def test_phased_tables_and_copies(run):
    out = run["out"]
    head = open(os.path.join(out, "phased.tsv")).readline().rstrip("\n").split("\t")
    rows = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in open(os.path.join(out, "phased.tsv")).readlines()[1:]]
    assert {r["role"] for r in rows} == {"child", "father", "mother"}
    fr = np.array([float(r["frac"]) for r in rows if r["frac"] != "NA"])
    assert fr.min() >= -0.05 and fr.max() <= 1.05
    # a normal chromosome: one maternal and one paternal copy; the trisomy's extra copy is paternal on chr12 (2 and 1)
    c3 = [r for r in rows if r["role"] == "child" and r["chrom"] == "chr3" and r["copies_tagged"] != "NA"]
    assert abs(np.median([float(r["copies_tagged"]) for r in c3]) - 1) < 0.1 and abs(np.median([float(r["copies_other"]) for r in c3]) - 1) < 0.1
    c12 = [r for r in rows if r["role"] == "child" and r["chrom"] == "chr12" and r["copies_tagged"] != "NA"]
    assert abs(np.median([float(r["copies_other"]) for r in c12]) - 1.3) < 0.12 and abs(np.median([float(r["copies_tagged"]) for r in c12]) - 1.0) < 0.1
    assert os.path.exists(os.path.join(out, "phased_rejected.tsv"))
    assert "child_phased_sites" in run["summary"] and run["summary"]["child_phased_sites"] > 1000


def test_dense_mock_low_share_events_come_from_the_phased_scan(tmp_path):
    """At a real genome's site density, an 8% gain of a whole chromosome and an 8% loss of an arm - both under the depth's
    threshold - come out of the phased scan, typed by the depth's lean, with the parent of origin; nothing else is called."""
    paths = write_mock(str(tmp_path / "dense"), seed=3, sites_per_mb=1000, contigs=["chr15", "chr16", "chr17"], events=LOW_SHARE_EVENTS)
    trio = read_trios(paths["trios"])[0]
    res = run_trio(paths["vcf"], trio, str(tmp_path / "out"), gc_track=paths["gc"], figures=False, log=lambda s: None)
    ev = res["events"]
    assert sorted((e.sample, e.chrom, e.type, e.span) for e in ev) == [("KID", "chr16", "gain", "whole"), ("MOM", "chr16", "loss", "p")], [(e.sample, e.chrom, e.type, e.span, e.f) for e in ev]
    g = [e for e in ev if e.sample == "KID"][0]
    assert g.source == "phased" and abs(g.f - 0.08) < 0.04 and g.origin_phase == "extra copy paternal" and "under its own threshold" in g.note
    l = [e for e in ev if e.sample == "MOM"][0]
    assert l.source == "phased" and abs(l.f - 0.08) < 0.04


def test_meiotic_stage_from_the_centromere_and_the_crossovers(tmp_path):
    """Three whole-chromosome trisomies at a real genome's site density: a meiosis II error (isodisomic at the centromere, one
    crossover), a meiosis I error (heterodisomic at the centromere, one crossover) and a mitotic duplication (one homologue
    throughout), each classified from the auxiliary track's state along the chromosome, with the crossover placed."""
    from triokaryo.mock import MEIOSIS_EVENTS
    # three normal chromosomes keep the autosomal median diploid; no GC bias, since three trisomies covering 40% of the bins would
    # confound a GC running median (no real child carries that)
    paths = write_mock(str(tmp_path / "meiosis"), seed=11, sites_per_mb=300, contigs=["chr10", "chr11", "chr12", "chr13", "chr16", "chr17"], events=MEIOSIS_EVENTS,
                       gc_beta=(0.0, 0.0, 0.0))
    trio = read_trios(paths["trios"])[0]
    res = run_trio(paths["vcf"], trio, str(tmp_path / "out"), figures=False, log=lambda s: None)
    ev = {(e.sample, e.chrom): e for e in res["events"]}
    assert sorted(ev) == [("KID", "chr13"), ("KID", "chr16"), ("KID", "chr17")], [(e.sample, e.chrom, e.type, e.span) for e in res["events"]]
    e13, e16, e17 = ev[("KID", "chr13")], ev[("KID", "chr16")], ev[("KID", "chr17")]
    assert all(e.type == "gain" and e.span == "whole" for e in (e13, e16, e17))
    assert e13.origin_phase == "extra copy maternal" and e13.centromere == "isodisomic" and e13.stage == "meiosis II"
    assert e13.n_crossovers == 1 and abs(float(e13.crossovers) - 60) < 3, (e13.n_crossovers, e13.crossovers)
    assert e16.origin_phase == "extra copy paternal" and e16.centromere == "heterodisomic" and e16.stage == "meiosis I"
    assert e16.n_crossovers == 1 and abs(float(e16.crossovers) - 70) < 3, (e16.n_crossovers, e16.crossovers)
    assert e17.centromere == "isodisomic" and e17.stage.startswith("mitotic") and e17.n_crossovers == 0 and e17.hetero_share < 0.05
    assert e13.homologues.startswith("the two maternal copies differ over") and e17.homologues.startswith("the two maternal copies are one homologue")
    assert e13.crossover_states == "iso>hetero" and e16.crossover_states == "hetero>iso"
    xo = [l.split("\t") for l in open(os.path.join(str(tmp_path / "out"), "crossovers.tsv")).read().splitlines()]
    rows = [dict(zip(xo[0], r)) for r in xo[1:]]
    assert [(r["chrom"], r["from_state"], r["to_state"], r["parent"]) for r in rows] == [("chr13", "isodisomic", "heterodisomic", "maternal"), ("chr16", "heterodisomic", "isodisomic", "paternal")]


def test_parent_event_is_read_folded_across_the_childs_crossover(tmp_path):
    """The transmitted-allele sign flips at the child's crossover in that parent's meiosis. A mother's mosaic gain of the whole of
    chromosome 8 spanning a crossover at 70 Mb is read folded onto the sign runs: its phased cell fraction matches the planted one
    instead of cancelling, the crossover is placed, and the depth call is not doubted."""
    from triokaryo.mock import DEFAULT_EVENTS, MOM, write_mock
    from triokaryo.pedigree import read_trios
    from triokaryo.phase import parent_origin_text, transmitted_sign_runs
    from triokaryo.pipeline import run_trio
    rng = np.random.default_rng(5)
    y = np.concatenate([0.05 + rng.normal(0, 0.01, 30), -0.05 + rng.normal(0, 0.01, 30)])
    y[10] = np.nan
    sign, flips = transmitted_sign_runs(y)
    assert sign[10] == 0 and (sign[:10] == 1).all() and (sign[11:30] == 1).all() and (sign[30:] == -1).all() and flips == [30]
    assert parent_origin_text("gain", 1, [70.2]) == "the duplicated homologue is, up to 70.2 Mb, the one passed to the child; the transmitted homologue switches at 70.2 Mb (the child's crossover)"
    assert parent_origin_text("loss", -1, []) == "the lost homologue is the one passed to the child"
    ev = [dict(e, f=0.30) for e in DEFAULT_EVENTS if e["member"] == MOM and e["chrom"] == "chr8"]
    m = write_mock(str(tmp_path / "sw"), seed=11, sites_per_mb=300, events=ev, contigs=("chr8", "chr9", "chr11", "chr12"), switch_maternal={"chr8": 70_000_000})
    trio = read_trios(m["trios"])[0]
    res = run_trio(m["vcf"], trio, str(tmp_path / "out"), gc_track=m["gc"], figures=False, log=lambda s: None)
    (e,) = [e for e in res["events"] if e.sample == "MOM" and e.chrom == "chr8" and e.type == "gain"]
    assert abs(e.f_phase - 0.30) < 0.08 and "artefact" not in e.note, (e.f_phase, e.note)
    assert e.n_crossovers == 1 and abs(float(e.crossovers) - 70.0) < 3 and e.crossover_states in ("transmitted>untransmitted", "untransmitted>transmitted"), (e.crossovers, e.crossover_states)
    assert "switches at" in e.origin_phase and e.origin_phase.startswith("the duplicated homologue is, up to"), e.origin_phase
    assert np.isfinite(e.start_fine) or e.start == 0
