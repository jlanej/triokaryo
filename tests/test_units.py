import numpy as np

from triokaryo.genome import genome
from triokaryo.model import d_hat, f_from_d, f_from_lrr
from triokaryo.pedigree import read_trios, sex_of
from triokaryo.segment import binary_segmentation, merge_similar, refine_boundaries


def test_genome_names_and_arms():
    G = genome()
    assert G.normalize("1") == "chr1" and G.normalize("x") == "chrX" and G.normalize("chrUn_x") is None
    assert G.arm("chr1", 1000) == "p" and G.arm("chr1", 200_000_000) == "q"
    assert G.is_par("chrX", 100_000) and not G.is_par("chrX", 50_000_000)


def test_band_estimator_recovers_d():
    rng = np.random.default_rng(1)
    n = 400
    dp = rng.poisson(30, n)
    for d_true in (0.0, 0.1, 1 / 6):
        side = rng.random(n) < 0.5
        alt = rng.binomial(dp, np.where(side, 0.5 + d_true, 0.5 - d_true))
        d, llr = d_hat(alt, dp)
        assert abs(d - d_true) < 0.02
        if d_true:
            assert llr > 20
    assert abs(f_from_d(1 / 6, "gain") - 1.0) < 0.01 and abs(f_from_lrr(np.log2(1.5), "gain") - 1.0) < 1e-9
    assert abs(f_from_d(0.2, "LOH") - 0.4) < 1e-9 and abs(f_from_lrr(np.log2(0.5), "loss") - 1.0) < 1e-9


def test_segmentation_finds_one_step():
    rng = np.random.default_rng(2)
    y = rng.normal(0, 0.03, 120)
    y[40:80] += 0.2
    segs = merge_similar(y, binary_segmentation(y, 5, 5.0), 0.03)
    assert [s for s in segs if 38 <= s[0] <= 42 and 78 <= s[1] <= 82]
    # a six-bin step: the boundaries land on the step's edges exactly once refined
    y2 = rng.normal(0, 0.03, 60)
    y2[20:26] -= 1.0
    segs2 = merge_similar(y2, refine_boundaries(y2, binary_segmentation(y2, 5, 5.0, 0.03)), 0.03)
    assert (20, 26) in segs2, segs2


def test_trios_file_and_ped(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("#kid\tdad\tmom\tkid_sex\tdad_sex\tmom_sex\nA\tB\tC\t1\t1\t2\n")
    (t,) = read_trios(str(p))
    assert t.members == ("A", "B", "C") and t.sexes == ("M", "M", "F")
    q = tmp_path / "t.ped"
    q.write_text("F1\tA\tB\tC\t1\t2\nF1\tB\t0\t0\t1\t1\nF1\tC\t0\t0\t2\t1\n")
    (t2,) = read_trios(str(q))
    assert t2.members == ("A", "B", "C") and t2.sexes == ("M", "M", "F") and t2.family == "F1"
    assert sex_of("female") == "F" and sex_of("0") == ""


def test_merged_absent_homref_counts_as_confident_parent():
    import numpy as np
    from triokaryo.trio import informative, mie_rate
    from triokaryo.vcfscan import Sites
    n = 40
    pos = np.arange(1, n + 1) * 1000
    dp = np.array([np.full(n, 30), np.zeros(n, int), np.full(n, 30)])          # the father absent from his VCF: hom-ref by the merge
    alt = np.array([np.full(n, 10), np.zeros(n, int), np.full(n, 30)])
    gt = np.array([np.ones(n, int), np.zeros(n, int), np.full(n, 2)])           # child het, father 0/0 (merge), mother 1/1
    gq = np.array([np.full(n, 99), np.full(n, -2), np.full(n, 99)])
    sites = Sites("chr1", pos, dp, alt, gt, gq, np.zeros(n, bool))
    k, d, mat = informative(sites, 1, 10 ** 9, 8, 20)
    assert len(k) == n and mat.all()
    gq2 = gq.copy(); gq2[1] = -1                                                # absent with no convention (a joint VCF's missing call): not confident
    k2, _, _ = informative(Sites("chr1", pos, dp, alt, gt, gq2, np.zeros(n, bool)), 1, 10 ** 9, 8, 20)
    assert len(k2) == 0
    assert mie_rate(sites, 1, 10 ** 9, 8, 20) == 0.0
