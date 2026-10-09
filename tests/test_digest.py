"""The digest: the cohort's largest events on one self-contained page, ranked by impact, the figures embedded, the long tail
counted and left in the tables."""
import os

from triokaryo.cli import main
from triokaryo.mock import write_mock
from triokaryo.pedigree import read_trios
from triokaryo.pipeline import run_trio


def _run(tmp_path, name, seed, figures, **kw):
    """A mock trio run without a panel (the mock has no shared depth structure to remove; the shipped panel, from real genomes,
    would add it): the digest is told so (--allow-first-pass)."""
    m = write_mock(str(tmp_path / name), seed=seed, **kw)
    out = str(tmp_path / (name + "_run"))
    run_trio(m["vcf"], read_trios(m["trios"])[0], out, gc_track=m["gc"], figures=figures, log=lambda s: None, events_path=m.get("events"))
    return out


def test_digest_ranks_the_largest_events_and_embeds_the_figures(tmp_path):
    run = _run(tmp_path, "m", 2, True)                         # the planted trio, with its figures
    null = _run(tmp_path, "n", 9, False, no_events=True, prefix="N")
    out = tmp_path / "digest"
    assert main(["digest", "--runs", run, null, "--out", str(out), "--top", "6", "--allow-first-pass"]) == 0
    page = open(out / "digest.html").read()
    assert "<img src=\"data:image/png;base64," in page and "src=\"http" not in page and "<script" not in page
    rows = [l.rstrip("\n").split("\t") for l in open(out / "digest.tsv")]
    hdr, rows = rows[0], rows[1:]
    col = lambda r, c: r[hdr.index(c)]
    assert len(rows) == 6 and [col(r, "rank") for r in rows] == [str(i + 1) for i in range(6)]
    # the whole chromosomes first: the trisomy 21 (a copy-number change in every cell) before the uniparental disomies, those before
    # the mosaic +12 and the mother's mosaic +8
    first = [(col(r, "sample"), col(r, "chrom"), col(r, "type"), col(r, "span")) for r in rows[:5]]
    assert first[0] == ("KID", "chr21", "gain", "whole"), first
    assert all(col(r, "span") == "whole" for r in rows[:5]) and {f[1] for f in first[1:3]} == {"chr7", "chr15"}, first
    assert {f[:2] for f in first[3:5]} == {("KID", "chr12"), ("MOM", "chr8")}, first                 # the mosaic whole chromosomes after
    assert "an extra copy of chromosome 21" in page and "New in the child" in page
    left = [l.rstrip("\n").split("\t") for l in open(out / "digest_left_out.tsv")]
    lh, left = left[0], left[1:]
    why = {(r[lh.index("sample")], r[lh.index("chrom")], r[lh.index("type")]): r[lh.index("left_out")] for r in left}
    assert any(k[1] == "chr2" and "under 10 Mb" in v for k, v in why.items()), why      # the 6-Mb deletion (father and child): below the floor
    assert not any(r[lh.index("sample")].startswith("N") and "beyond" in r[lh.index("left_out")] for r in left)
    assert "left out by the floors" in page or "beyond the first" in page


def test_digest_without_figures_draws_from_the_bins_and_refuses_the_first_pass(tmp_path):
    run = _run(tmp_path, "m", 2, False)
    out = tmp_path / "d2"
    assert main(["digest", "--runs", run, "--out", str(out), "--top", "3", "--allow-first-pass"]) == 0
    page = open(out / "digest.html").read()
    assert page.count("drawn from the bins") >= 1 and page.count("data:image/png;base64,") >= 2
    m = write_mock(str(tmp_path / "fp"), seed=3)
    fp = str(tmp_path / "fp_run")
    run_trio(m["vcf"], read_trios(m["trios"])[0], fp, gc_track=m["gc"], figures=False, log=lambda s: None)   # no panel: the first pass
    try:
        main(["digest", "--runs", fp, "--out", str(tmp_path / "d3")])
        assert False, "the first pass must be refused"
    except SystemExit as e:
        assert "first-pass" in str(e)
    assert main(["digest", "--runs", fp, "--out", str(tmp_path / "d4"), "--allow-first-pass"]) == 0
