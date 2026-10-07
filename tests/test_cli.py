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
