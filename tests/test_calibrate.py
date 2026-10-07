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


def test_calibration_whole_chromosome_events_and_false_positive_table(tmp_path):
    """--whole adds a maternal meiosis I trisomy and a heterodisomy: found, staged and attributed at f = 1; the false-positive table and
    the replicate interval are written."""
    rows = run_calibration(str(tmp_path), fractions=(1.0,), depths=(30,), sizes=(10,), replicates=2, seed=4, log=lambda s: None, figures=False, whole=True,
                           contigs=["chr1", "chr2", "chr3", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18"])
    whole = [r for r in rows if r["size_mb"] == "whole"]
    assert len(whole) == 4 and all(r["detected"] for r in whole), [(r["type"], r["detected"], r["type_called"]) for r in whole]
    tri = [r for r in whole if r["type"] == "gain"]
    assert all(r["stage_ok"] is True and r["origin_ok"] and r["origin"] == "extra copy maternal" for r in tri), [(r["stage"], r["origin"]) for r in tri]
    upd = [r for r in whole if r["type"] == "UPD"]
    assert all(r["type_called"] == "UPD" and r["origin_ok"] for r in upd), [(r["type_called"], r["origin"]) for r in upd]
    md = open(os.path.join(str(tmp_path), "calibration.md")).read()
    assert "## Whole-chromosome events" in md and "Wilson 95% interval" in md and "## Calls on unaffected chromosomes" in md
    assert os.path.exists(os.path.join(str(tmp_path), "false_positives.tsv"))
