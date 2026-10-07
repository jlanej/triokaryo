"""Calibration on the simulator: detection and cell-fraction accuracy as a function of cell fraction, depth and event size.

For each cell fraction and depth, one simulated trio carries a gain, a loss and a copy-neutral LOH of each size on
separate autosomes (the duplicated, lost or replaced copy paternal). Each planted event is matched to the child's calls
(same type, intersection at least half the planted length) and the run records whether it was found, by which source,
the three cell-fraction estimates and whether the parent of origin was named correctly. The table and the figure show
where the method's floor lies for a given design: detection per type, size and cell fraction at each depth, and the
error of each estimate. No method is changed here; the grid is a description of performance."""
import os
import time

import numpy as np

from .genome import genome as load_genome
from .mock import KID, write_mock
from .pedigree import read_trios
from .pipeline import run_trio

DEFAULT_FRACTIONS = (0.05, 0.1, 0.2, 0.3, 0.5, 1.0)
DEFAULT_DEPTHS = (30,)
DEFAULT_SIZES = (5, 10, 20, 50)
TYPES = (("gain", {"pat": +1}, "extra copy paternal"), ("loss", {"pat": -1}, "paternal copy lost"), ("LOH", {"mat": -1, "pat": +1}, "paternal copy retained (maternal replaced)"))
# host chromosomes for the planted events, one per (type, size): the event starts 20 Mb into the q arm, clear of the centromere
HOSTS = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16")


def planted_events(genome, fraction, sizes):
    """The events of one calibration trio: a gain, a loss and a CN-LOH of each size, each on its own chromosome."""
    out = []
    hosts = iter(HOSTS)
    for size in sizes:
        for kind, delta, origin in TYPES:
            c = next(hosts)
            start = genome.p_end[c] + 20_000_000
            end = min(start + int(size * 1e6), genome.length[c])
            out.append(dict(member=KID, chrom=c, start=start, end=end, f=float(fraction), delta=delta, type=kind, origin=origin, inherited=False,
                            label="%s of %g Mb at f = %.2f" % (kind, size, fraction), size=size))
    return out


def run_calibration(out, fractions=DEFAULT_FRACTIONS, depths=DEFAULT_DEPTHS, sizes=DEFAULT_SIZES, replicates=1, seed=1, sites_per_mb=60, log=None, figures=True,
                    contigs=None):
    """Runs the grid, writes calibration.tsv (one row per planted event), calibration.md and calibration.{png,svg,pdf}; returns the
    rows. sites_per_mb: the simulated site density (real WGS has about 1,000 PASS SNVs per Mb; 60 is quick but understates the phased
    scan); contigs: restrict the simulation to these chromosomes, which must include the hosts (chr1 onwards, one per type and
    size) and some unaffected ones for the autosomal median."""
    log = log or (lambda s: None)
    os.makedirs(out, exist_ok=True)
    G = load_genome("grch38")
    rows = []
    t0 = time.time()
    k = 0
    for rep in range(replicates):
        for depth in depths:
            for f in fractions:
                k += 1
                tag = "d%g_f%g_r%d" % (depth, f, rep)
                events = planted_events(G, f, sizes)
                d = os.path.join(out, "runs", tag)
                m = write_mock(os.path.join(d, "mock"), seed=seed + 1000 * rep + k, sites_per_mb=sites_per_mb, events=events, depth=(float(depth),) * 3,
                               contigs=contigs)
                trio = read_trios(m["trios"])[0]
                res = run_trio(m["vcf"], trio, os.path.join(d, "run"), gc_track=m["gc"], figures=False, log=lambda s: None)
                calls = [e for e in res["events"] if e.role == "child"]
                for ev in events:
                    hits = [e for e in calls if e.chrom == ev["chrom"] and e.type == ev["type"]
                            and min(e.end, ev["end"]) - max(e.start, ev["start"]) >= 0.5 * (ev["end"] - ev["start"])]
                    others = [e for e in calls if e.chrom == ev["chrom"] and e not in hits]
                    hit = max(hits, key=lambda e: e.end - e.start) if hits else None
                    rows.append(dict(replicate=rep, depth=depth, f=f, type=ev["type"], size_mb=ev["size"], chrom=ev["chrom"], detected=bool(hit),
                                     source=hit.source if hit else "", f_lrr=hit.f_lrr if hit else np.nan, f_baf=hit.f_baf if hit else np.nan,
                                     f_phase=hit.f_phase if hit else np.nan, f_called=hit.f if hit else np.nan,
                                     origin_ok=bool(hit and (hit.origin_phase or hit.origin) == ev["origin"]), origin=(hit.origin_phase or hit.origin) if hit else "",
                                     span_called=(hit.end - hit.start) / 1e6 if hit else np.nan, other_calls_on_chrom=len(others)))
                extra = [e for e in calls if not any(e.chrom == ev["chrom"] for ev in events)]
                log("%s: %d of %d planted events found, %d other child call(s); %.0f s" % (tag, sum(1 for r in rows[-len(events):] if r["detected"]), len(events), len(extra), time.time() - t0))
    from .report import write_tsv
    cols = ["replicate", "depth", "f", "type", "size_mb", "chrom", "detected", "source", "f_called", "f_lrr", "f_baf", "f_phase", "origin_ok", "origin", "span_called", "other_calls_on_chrom"]
    write_tsv(os.path.join(out, "calibration.tsv"), cols, rows)
    write_summary(out, rows, fractions, depths, sizes)
    if figures:
        fig_calibration(out, rows, fractions, depths, sizes)
    return rows


def write_summary(out, rows, fractions, depths, sizes):
    """calibration.md: detection rate per type, depth, size and cell fraction, and the median absolute error of each estimate."""
    lines = ["# Calibration", "", "Detection rate of planted events in the child (rows: cell fraction; columns: event size in Mb), per type and depth.", ""]
    for depth in depths:
        for kind, _, _ in TYPES:
            lines.append("## %s, %g×" % (kind, depth))
            lines.append("")
            lines.append("| f | " + " | ".join("%g Mb" % s for s in sizes) + " |")
            lines.append("| --- | " + " | ".join("---" for _ in sizes) + " |")
            for f in fractions:
                cells = []
                for s in sizes:
                    sel = [r for r in rows if r["depth"] == depth and r["type"] == kind and r["f"] == f and r["size_mb"] == s]
                    cells.append("%.0f%%" % (100 * np.mean([r["detected"] for r in sel])) if sel else "")
                lines.append("| %.2f | " % f + " | ".join(cells) + " |")
            lines.append("")
    lines.append("## Cell-fraction error of detected events (median |estimate − f|)")
    lines.append("")
    lines.append("| type | depth | f_called | f_lrr | f_baf | f_phase | parent of origin correct |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for depth in depths:
        for kind, _, _ in TYPES:
            sel = [r for r in rows if r["depth"] == depth and r["type"] == kind and r["detected"]]
            if not sel:
                continue
            err = lambda key: np.nanmedian([abs(r[key] - r["f"]) for r in sel if np.isfinite(r[key])]) if any(np.isfinite(r[key]) for r in sel) else np.nan  # noqa: E731
            lines.append("| %s | %g× | %.3f | %.3f | %.3f | %.3f | %.0f%% |" % (kind, depth, err("f_called"), err("f_lrr"), err("f_baf"), err("f_phase"),
                                                                              100 * np.mean([r["origin_ok"] for r in sel])))
    with open(os.path.join(out, "calibration.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def fig_calibration(out, rows, fractions, depths, sizes):
    """calibration.{png,svg,pdf}: per depth (columns) and type (rows), the detection rate over size and cell fraction, with
    the source of the call marked (D depth, B bands, P phased) where every replicate agrees."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .plots import MM
    nd, nt = len(depths), len(TYPES)
    fig, axes = plt.subplots(nt, nd, figsize=(max(60 * nd, 80) * MM, (30 * nt + 15) * MM), squeeze=False)
    for j, depth in enumerate(depths):
        for i, (kind, _, _) in enumerate(TYPES):
            ax = axes[i][j]
            grid = np.full((len(fractions), len(sizes)), np.nan)
            src = [["" for _ in sizes] for _ in fractions]
            for a, f in enumerate(fractions):
                for b, s in enumerate(sizes):
                    sel = [r for r in rows if r["depth"] == depth and r["type"] == kind and r["f"] == f and r["size_mb"] == s]
                    if sel:
                        grid[a, b] = np.mean([r["detected"] for r in sel])
                        srcs = {r["source"] for r in sel if r["detected"]}
                        src[a][b] = {"depth": "D", "bands": "B", "phased": "P"}.get(srcs.pop(), "") if len(srcs) == 1 else ("mixed" if srcs else "")
            ax.imshow(grid, vmin=0, vmax=1, cmap="Greys", origin="lower", aspect="auto")
            for a in range(len(fractions)):
                for b in range(len(sizes)):
                    if np.isfinite(grid[a, b]):
                        ax.text(b, a, "%.0f%%\n%s" % (100 * grid[a, b], src[a][b]), ha="center", va="center", fontsize=5, color="white" if grid[a, b] > 0.5 else "black")
            ax.set_xticks(range(len(sizes)))
            ax.set_xticklabels(["%g" % s for s in sizes], fontsize=6)
            ax.set_yticks(range(len(fractions)))
            ax.set_yticklabels(["%.2f" % f for f in fractions], fontsize=6)
            ax.set_title("%s, %g×" % (kind, depth), fontsize=7, loc="left")
            if i == nt - 1:
                ax.set_xlabel("event size (Mb)", fontsize=6)
            if j == 0:
                ax.set_ylabel("cell fraction f", fontsize=6)
    fig.tight_layout()
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(out, "calibration.%s" % ext), bbox_inches="tight")
    plt.close(fig)
