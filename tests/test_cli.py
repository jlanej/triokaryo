import os
import subprocess
import sys

from triokaryo.cli import main


def test_cli_run_and_cohort(mock, tmp_path):
    out = tmp_path / "KID"
    rc = main(["run", "--vcf", mock["vcf"], "--pedigree", mock["trios"], "--child", "KID", "--gc-track", mock["gc"], "--events", mock["events"],
               "--out", str(out), "--no-figures", "--thin", "2"])
    assert rc == 0 and (out / "events.tsv").exists()
    text = (out / "events.tsv").read_text()
    assert "chr21" in text and "extra copy maternal" in text
    coh = tmp_path / "cohort"
    rc = main(["cohort", "--runs", str(out), "--events", mock["events"], "--out", str(coh)])
    assert rc == 0
    for f in ("events.all.tsv", "summary.all.tsv", "concordance.tsv", "index.html"):
        assert (coh / f).exists()
    conc = (coh / "concordance.tsv").read_text().splitlines()
    assert len(conc) == 9 and all("matched" in l for l in conc[1:])


def test_cli_members_given_outright(mock, tmp_path):
    out = tmp_path / "direct"
    rc = main(["run", "--vcf", mock["vcf"], "--child", "KID", "--father", "DAD", "--mother", "MOM", "--sex", "M,M,F", "--out", str(out), "--no-figures", "--thin", "3"])
    assert rc == 0 and "chr7" in (out / "events.tsv").read_text()


def test_cli_entry_point_help():
    r = subprocess.run([sys.executable, "-m", "triokaryo.cli", "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "gc-track" in r.stdout


def test_cohort_report_guide_and_rebuild(mock, run, tmp_path):
    """The cohort report (counts, the landscape figure, every event, the trios, the concordance, the segments set aside) and
    the guide beside it; a run's page rebuilt from its tables alone."""
    import os
    from triokaryo.cli import main
    out = tmp_path / "cohort"
    assert main(["cohort", "--runs", run["out"], "--events", mock["events"], "--out", str(out)]) == 0
    for f in ("index.html", "guide.html", "events.all.tsv", "summary.all.tsv", "concordance.tsv", "flags.tsv", "rejected.all.tsv", "figures/landscape.png",
              "figures/landscape.txt", "figures/patterns_copy.png", "figures/patterns_disomy.png", "figures/legends/landscape_legend.png"):
        assert os.path.exists(os.path.join(str(out), f)), f
    page = open(os.path.join(str(out), "index.html")).read()
    assert "sortTable" in page and 'id="filter"' in page and "guide.html" in page and "heterodisomies" in page and "chrom_chr21.png" in page
    guide = open(os.path.join(str(out), "guide.html")).read()
    assert "below 1/2 means the paternal homologue is in excess" in guide and "Pattern cards 2" in guide and "<code>origin_phase</code>" in guide
    head = open(os.path.join(str(out), "events.all.tsv")).readline().split("\t")
    assert "homologues" in head and "source" in head and "start_fine" in head
    # the run's page again, from its tables
    before = open(os.path.join(run["out"], "index.html")).read()
    assert main(["report", "--runs", run["out"]]) == 0
    after = open(os.path.join(run["out"], "index.html")).read()
    assert "guide.html" in after and "below 1/2 means the paternal homologue is in excess" in after and after.count("<figure>") == before.count("<figure>")
    assert os.path.exists(os.path.join(run["out"], "guide.html"))
    sc = open(os.path.join(run["out"], "figures", "genome.txt")).read()
    assert "keys: depth,step,baf,phased," in sc
    g = tmp_path / "guide.html"
    assert main(["guide", "--out", str(g)]) == 0 and g.stat().st_size > 100_000
