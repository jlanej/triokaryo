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


def test_single_bin_spike_is_smoothed_and_a_plateau_kept():
    """A bin far above the noise (a germline CNV, collapsed repeats) is shrunk to its neighbours' level before segmentation; a plateau
    of two bins or more keeps its level, so no event of the minimum length is altered; the spike alone makes no segment."""
    from triokaryo.segment import smooth_outliers
    rng = np.random.default_rng(3)
    sd = 0.02
    y = rng.normal(0, sd, 60)
    y[20] = 0.6                                                    # one bin at a 1.5-fold depth
    y[40:45] += 0.25                                               # a five-bin gain
    s = smooth_outliers(y, sd)
    assert abs(s[20]) < 3 * sd, s[20]
    assert np.allclose(s[40:45], y[40:45]) and np.allclose(np.delete(s, 20), np.delete(y, 20))
    segs = merge_similar(s, refine_boundaries(s, binary_segmentation(s, 5, 5.0, sd)), sd)
    of = lambda i: next((a, b) for a, b in segs if a <= i < b)
    assert abs(np.mean(s[slice(*of(20))])) < 0.05 and np.mean(s[slice(*of(42))]) > 0.2, segs
    y2 = rng.normal(0, sd, 40)
    y2[30:32] += 0.5                                               # a two-bin plateau: each bin has a neighbour at its level
    assert np.allclose(smooth_outliers(y2, sd), y2)


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


def test_inheritance_needs_reciprocal_overlap_but_external_matching_does_not():
    """A small parental event inside a large event of the child is not 'inherited'; a caller's fragment inside an event still matches."""
    import types
    from triokaryo.segment import Event
    from triokaryo.trio import read_trio
    big = Event("KID", "child", "chr1", 0, 100_000_000, "stretch", "gain", f_lrr=1.0)
    small = Event("DAD", "father", "chr1", 10_000_000, 20_000_000, "stretch", "gain", f_lrr=1.0)
    assert big.overlap(small) == 1.0 and abs(big.reciprocal_overlap(small) - 0.1) < 1e-9
    scan = types.SimpleNamespace(chroms={}, sites=lambda c: None)
    read_trio([big, small], scan, genome(), 8, 20)
    assert big.inheritance.startswith("new") and small.inheritance == "not passed to the child"
    same = Event("DAD", "father", "chr1", 0, 90_000_000, "stretch", "gain", f_lrr=1.0)
    read_trio([big, same], scan, genome(), 8, 20)
    assert big.inheritance == "inherited from the father" and same.inheritance == "passed to the child"


def test_run_of_homozygosity_gets_no_parent_of_origin_from_the_phased_track():
    """An event annotated as a run of homozygosity keeps its phased statistics but receives no origin_phase or homologue count."""
    import types
    from triokaryo.phase import PhasedTrack, annotate_events
    from triokaryo.segment import Event
    from triokaryo.vcfscan import Sites
    n = 200
    pos = np.arange(1, n + 1) * 10_000
    dp = np.full((3, n), 30)
    alt = np.array([np.full(n, 24), np.zeros(n, int), np.full(n, 30)])             # the child's maternal (alt) allele at 0.8
    gt = np.array([np.ones(n, int), np.zeros(n, int), np.full(n, 2)])
    gq = np.full((3, n), 99)
    sites = Sites("chr1", pos, dp, alt, gt, gq, np.zeros(n, bool))
    scan = types.SimpleNamespace(chroms={"chr1": sites}, sites=lambda c: sites)
    idx = np.arange(n)
    frac = alt[0] / 30.0
    w = np.arange(0, n, 20)
    track = PhasedTrack("chr1", idx, frac, np.ones(n, bool), pos[w], pos[np.minimum(w + 19, n - 1)], pos[w + 10], np.full(len(w), 20),
                        np.full(len(w), 600), np.full(len(w), 0.8), np.full(len(w), 0.02), np.zeros(len(w), int), np.full(len(w), 0.8),
                        np.zeros(len(w)), np.ones(len(w)), np.ones(len(w)))
    bins = types.SimpleNamespace(bin_size=1_000_000)
    for note, expect in (("", True), ("a run of homozygosity (both copies identical by descent), not a uniparental disomy: no Mendelian errors", False)):
        ev = Event("KID", "child", "chr1", 0, 2_000_000, "stretch", "LOH", f_baf=1.0, source="bands", note=note)
        annotate_events([ev], [{"chr1": track}, {}, {}], scan, bins, genome())
        assert ev.n_phased == n and abs(ev.phase_shift - 0.3) < 1e-9
        assert bool(ev.origin_phase) is expect and bool(ev.homologues) is False


def test_trimmed_mean_and_noise_floor_keep_a_quantised_depth_track_segmentable():
    """A bin depth that is the median of hundreds of integer depths is quantised to one read, so a track of such bins can have
    first differences that are all zero and a robust noise scale of zero, which would disable the segmentation; the trimmed mean
    is continuous and the counting-noise floor bounds the scale from below."""
    from triokaryo.model import lrr_noise_floor, trimmed_mean
    from triokaryo.segment import binary_segmentation, robust_sd
    rng = np.random.default_rng(5)
    med = [float(np.median(rng.poisson(30, 600))) for _ in range(50)]
    tm = [trimmed_mean(rng.poisson(30, 600)) for _ in range(50)]
    assert len(set(med)) <= 3 and len(set(tm)) >= 40                      # the median takes a few integer values; the trimmed mean does not
    assert trimmed_mean([1, 2, 3, 4, 100]) == 3.0 and trimmed_mean([1, 2, 100]) == 2.0 and trimmed_mean([]) != trimmed_mean([])   # NaN for no values
    y = np.zeros(100)
    y[40:60] = np.log2(1.5)
    assert robust_sd(y) == 0.0                                            # the degenerate scale
    floor = lrr_noise_floor(30, 600)
    assert 0.01 < floor < 0.015
    segs = binary_segmentation(y, 5, 5.0, sd=max(robust_sd(y), floor))
    assert (40, 60) in segs, segs
    assert binary_segmentation(y, 5, 5.0, sd=0.0) == [(0, 100)]            # without the floor nothing is segmented


def test_vcf_without_gq_is_analysed_on_depth_alone(tmp_path):
    """A caller that writes no GQ: the scan says so, every called genotype with reads counts as confident, and the planted events
    are still found; a VCF without AD is refused with a clear message."""
    import pytest
    from triokaryo.mock import write_mock
    from triokaryo.pipeline import run_trio
    from triokaryo.vcfscan import scan_vcf
    m = write_mock(str(tmp_path / "nogq"), seed=1, with_gq=False)
    trio = read_trios(m["trios"])[0]
    logs = []
    scan = scan_vcf(m["vcf"], trio.members, genome(), thin=4, log=logs.append)
    assert not scan.gq_in_header and any("no GQ" in l for l in logs)
    s = scan.sites("chr21")
    assert (s.gq[0][s.dp[0] > 0] == 99).all()
    res = run_trio(m["vcf"], trio, str(tmp_path / "out"), gc_track=m["gc"], figures=False, thin=2, log=logs.append)
    found = {(e.sample, e.chrom, e.type) for e in res["events"]}
    assert ("KID", "chr21", "gain") in found and ("KID", "chr7", "LOH") in found and ("KID", "chr18", "loss") in found
    # no AD at all: refused
    import gzip
    p = tmp_path / "noad.vcf"
    with gzip.open(m["vcf"], "rt") as fh, open(p, "w") as out:
        for line in fh:
            if line.startswith("##FORMAT=<ID=AD"):
                continue
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            f[8] = "GT:DP"
            f[9:] = [":".join(x.split(":")[i] for i in (0, 2)) for x in f[9:]]
            out.write("\t".join(f) + "\n")
            break                                                     # one record is enough: the header decides
    with pytest.raises(SystemExit, match="no AD"):
        scan_vcf(str(p), trio.members, genome())


def test_reference_bias_is_the_pooled_fraction():
    """A 2% deficit of alt reads at heterozygous sites is measured (the median of k/n would read 0 at this depth)."""
    import types
    from triokaryo.phase import ref_bias
    from triokaryo.vcfscan import Sites
    rng = np.random.default_rng(3)
    n = 7 * 400
    dp = np.full((3, n), 30)
    alt = np.array([rng.binomial(30, 0.48, n), np.zeros(n, int), np.zeros(n, int)])
    gt = np.array([np.ones(n, int), np.zeros(n, int), np.zeros(n, int)])
    gq = np.full((3, n), 99)
    sites = Sites("chr1", np.arange(1, n + 1) * 1000, dp, alt, gt, gq, np.zeros(n, bool))
    scan = types.SimpleNamespace(chroms={"chr1": sites})
    b = ref_bias(scan, genome(), 0, 8, 20)
    assert -0.03 < b < -0.01, b


def test_short_state_runs_are_not_crossovers():
    """Along a whole-chromosome event, a state run (two homologues or one) shorter than the minimum run is noise near the floor of
    the homologue test, not a crossover: it takes its neighbour's state."""
    from triokaryo.phase import flatten_short_runs
    s = np.array([0] * 20 + [1] * 3 + [0] * 20, bool)
    assert not flatten_short_runs(s, 5).any()
    s = np.array([0] * 20 + [1] * 8 + [0] * 20, bool)
    assert flatten_short_runs(s, 5).sum() == 8                                # a real segment keeps both of its crossovers
    s = np.array([1] * 2 + [0] * 30, bool)
    assert not flatten_short_runs(s, 5).any()
    s = np.array([0] * 3 + [1] * 3, bool)
    assert len(set(flatten_short_runs(s, 5).tolist())) == 1                   # two short runs: one state is left
    assert flatten_short_runs(np.zeros(0, bool), 5).size == 0
