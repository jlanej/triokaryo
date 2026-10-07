"""The calibration grid: planted gains, losses and CN-LOH matched to the child's calls, with the table, summary and figure."""
import os

from triokaryo.calibrate import run_calibration


def test_calibration_grid_runs_and_reports(tmp_path):
    rows = run_calibration(str(tmp_path), fractions=(1.0,), depths=(30,), sizes=(30,), replicates=1, seed=3, log=lambda s: None)
    assert len(rows) == 3 and {r["type"] for r in rows} == {"gain", "loss", "LOH"}
    assert all(r["detected"] and r["origin_ok"] and abs(r["f_called"] - 1.0) < 0.1 for r in rows), [(r["type"], r["detected"], r["origin"], r["f_called"]) for r in rows]
    for f in ("calibration.tsv", "calibration.md", "calibration.png", "calibration.svg", "calibration.pdf"):
        assert os.path.getsize(os.path.join(str(tmp_path), f)) > 0, f
    md = open(os.path.join(str(tmp_path), "calibration.md")).read()
    assert "## gain, 30×" in md and "| 1.00 | 100% |" in md and "parent of origin correct" in md
    head = open(os.path.join(str(tmp_path), "calibration.tsv")).readline().split("\t")
    assert head[:8] == ["replicate", "depth", "f", "type", "size_mb", "chrom", "detected", "source"]
