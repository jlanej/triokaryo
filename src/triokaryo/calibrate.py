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
# whole-chromosome events (--whole): a maternal meiosis I trisomy (both maternal homologues) and a maternal heterodisomy, with their expected stage
WHOLE = (("chr17", "gain", {"mat_other": +1}, "extra copy maternal", "meiosis I"), ("chr18", "UPD", {"pat": -1, "mat_other": +1}, "both copies maternal (heterodisomy)", "meiosis I"))


def planted_events(genome, fraction, sizes, whole=False):
    """The events of one calibration trio: a gain, a loss and a CN-LOH of each size, each on its own chromosome; with `whole`, a
    whole-chromosome maternal meiosis I trisomy and a maternal heterodisomy as well (mosaic below f = 1)."""
    out = []
    hosts = iter(HOSTS)
    for size in sizes:
        for kind, delta, origin in TYPES:
            c = next(hosts)
            start = genome.p_end[c] + 20_000_000
            end = min(start + int(size * 1e6), genome.length[c])
            out.append(dict(member=KID, chrom=c, start=start, end=end, f=float(fraction), delta=delta, type=kind, origin=origin, inherited=False,
                            label="%s of %g Mb at f = %.2f" % (kind, size, fraction), size=size, stage=""))
    if whole:
        for c, kind, delta, origin, stage in WHOLE:
            out.append(dict(member=KID, chrom=c, start=0, end=None, f=float(fraction), delta=delta, type=kind, origin=origin, inherited=False,
                            label="whole-chromosome %s at f = %.2f" % (kind, fraction), size="whole", stage=stage))
    return out


def _matches(call, ev, chrom_length):
    """A call matches a planted event when it has the type (a heterodisomy planted at f < 1 is read as copy-neutral LOH, which the
    method cannot tell apart from a mosaic heterodisomy) and covers at least half of it."""
    same_type = call.type == ev["type"] or (ev["type"] == "UPD" and call.type == "LOH")
    end = ev["end"] if ev["end"] is not None else chrom_length
    return call.chrom == ev["chrom"] and same_type and min(call.end, end) - max(call.start, ev["start"]) >= 0.5 * (end - ev["start"])


def _parent(label):
    return "maternal" if "maternal" in label else "paternal" if "paternal" in label else ""


def wilson(k, n, z=1.96):
    """The Wilson 95% interval of a rate k/n, as (low, high)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def run_calibration(out, fractions=DEFAULT_FRACTIONS, depths=DEFAULT_DEPTHS, sizes=DEFAULT_SIZES, replicates=1, seed=1, sites_per_mb=60, log=None, figures=True,
                    contigs=None, whole=False):
    """Runs the grid, writes calibration.tsv (one row per planted event), false_positives.tsv (the child's calls on unaffected
    chromosomes), calibration.md and calibration.{png,svg,pdf}; returns the rows. sites_per_mb: the simulated site density (real WGS
    has about 1,000 PASS SNVs per Mb; 60 is quick but understates the phased scan); contigs: restrict the simulation to these
    chromosomes, which must include the hosts (chr1 onwards, one per type and size; chr17 and chr18 with `whole`) and some
    unaffected ones for the autosomal median; whole: add a whole-chromosome maternal meiosis I trisomy and a heterodisomy."""
    log = log or (lambda s: None)
    os.makedirs(out, exist_ok=True)
    G = load_genome("grch38")
    rows, fps = [], []
    t0 = time.time()
    k = 0
    for rep in range(replicates):
        for depth in depths:
            for f in fractions:
                k += 1
                tag = "d%g_f%g_r%d" % (depth, f, rep)
                events = planted_events(G, f, sizes, whole=whole)
                d = os.path.join(out, "runs", tag)
                m = write_mock(os.path.join(d, "mock"), seed=seed + 1000 * rep + k, sites_per_mb=sites_per_mb, events=events, depth=(float(depth),) * 3,
                               contigs=contigs)
                trio = read_trios(m["trios"])[0]
                res = run_trio(m["vcf"], trio, os.path.join(d, "run"), gc_track=m["gc"], figures=False, log=lambda s: None)
                calls = [e for e in res["events"] if e.role == "child"]
                for ev in events:
                    hits = [e for e in calls if _matches(e, ev, G.length[ev["chrom"]])]
                    others = [e for e in calls if e.chrom == ev["chrom"] and e not in hits]
                    hit = max(hits, key=lambda e: e.end - e.start) if hits else None
                    origin = (hit.origin_phase or hit.origin) if hit else ""
                    rows.append(dict(replicate=rep, depth=depth, f=f, type=ev["type"], size_mb=ev["size"], chrom=ev["chrom"], detected=bool(hit),
                                     source=hit.source if hit else "", type_called=hit.type if hit else "", f_lrr=hit.f_lrr if hit else np.nan, f_baf=hit.f_baf if hit else np.nan,
                                     f_phase=hit.f_phase if hit else np.nan, f_called=hit.f if hit else np.nan,
                                     origin_ok=bool(hit and _parent(origin) == _parent(ev["origin"]) and _parent(origin)), origin=origin,
                                     stage=hit.stage if hit else "", stage_ok=bool(hit and ev["stage"] and hit.stage.startswith(ev["stage"])) if ev["stage"] else "",
                                     span_called=(hit.end - hit.start) / 1e6 if hit else np.nan, other_calls_on_chrom=len(others)))
                planted_chroms = {ev["chrom"] for ev in events}
                for e in calls:
                    if e.chrom not in planted_chroms:
                        fps.append(dict(replicate=rep, depth=depth, f=f, chrom=e.chrom, type=e.type, source=e.source, start=e.start, end=e.end, f_called=e.f, note=e.note))
                n_fp = sum(1 for r in fps if r["replicate"] == rep and r["depth"] == depth and r["f"] == f)
                log("%s: %d of %d planted events found, %d call(s) on unaffected chromosomes; %.0f s" % (tag, sum(1 for r in rows[-len(events):] if r["detected"]), len(events), n_fp, time.time() - t0))
    from .report import write_tsv
    cols = ["replicate", "depth", "f", "type", "size_mb", "chrom", "detected", "source", "type_called", "f_called", "f_lrr", "f_baf", "f_phase", "origin_ok", "origin", "stage",
            "stage_ok", "span_called", "other_calls_on_chrom"]
    write_tsv(os.path.join(out, "calibration.tsv"), cols, rows)
    write_tsv(os.path.join(out, "false_positives.tsv"), ["replicate", "depth", "f", "chrom", "type", "source", "start", "end", "f_called", "note"], fps)
    write_summary(out, rows, fractions, depths, sizes, fps=fps, replicates=replicates, whole=whole)
    if figures:
        fig_calibration(out, rows, fractions, depths, sizes)
    return rows


def _rate(sel, ci):
    k, n = sum(1 for r in sel if r["detected"]), len(sel)
    if not n:
        return ""
    if ci and n > 1:
        lo, hi = wilson(k, n)
        return "%.0f%% (%.0f–%.0f)" % (100 * k / n, 100 * lo, 100 * hi)
    return "%.0f%%" % (100 * k / n)


def write_summary(out, rows, fractions, depths, sizes, fps=(), replicates=1, whole=False):
    """calibration.md: detection rate per type, depth, size and cell fraction (with a Wilson 95% interval when there are
    replicates), the whole-chromosome events, the calls on unaffected chromosomes, and the median absolute error of each estimate."""
    ci = replicates > 1
    lines = ["# Calibration", "", "Detection rate of planted events in the child (rows: cell fraction; columns: event size in Mb), per type and depth%s." % (
        "; %d replicates, Wilson 95%% interval in parentheses" % replicates if ci else "; one replicate"), ""]
    for depth in depths:
        for kind, _, _ in TYPES:
            lines.append("## %s, %g×" % (kind, depth))
            lines.append("")
            lines.append("| f | " + " | ".join("%g Mb" % s for s in sizes) + " |")
            lines.append("| --- | " + " | ".join("---" for _ in sizes) + " |")
            for f in fractions:
                cells = [_rate([r for r in rows if r["depth"] == depth and r["type"] == kind and r["f"] == f and r["size_mb"] == s], ci) for s in sizes]
                lines.append("| %.2f | " % f + " | ".join(cells) + " |")
            lines.append("")
    if whole:
        lines.append("## Whole-chromosome events")
        lines.append("")
        lines.append("A maternal meiosis I trisomy (both maternal homologues) and a maternal heterodisomy; a heterodisomy below f = 1 is read as copy-neutral LOH, "
                     "which counts as found. Stage: the share of found trisomies staged meiosis I.")
        lines.append("")
        lines.append("| depth | f | trisomy found | staged meiosis I | origin correct | heterodisomy found | read as |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for depth in depths:
            for f in fractions:
                tri = [r for r in rows if r["depth"] == depth and r["f"] == f and r["size_mb"] == "whole" and r["type"] == "gain"]
                upd = [r for r in rows if r["depth"] == depth and r["f"] == f and r["size_mb"] == "whole" and r["type"] == "UPD"]
                if not tri and not upd:
                    continue
                found_tri = [r for r in tri if r["detected"]]
                read_as = sorted({r["type_called"] for r in upd if r["detected"]})
                lines.append("| %g× | %.2f | %s | %s | %s | %s | %s |" % (
                    depth, f, _rate(tri, ci), ("%.0f%%" % (100 * np.mean([r["stage_ok"] is True for r in found_tri]))) if found_tri else "",
                    ("%.0f%%" % (100 * np.mean([r["origin_ok"] for r in found_tri]))) if found_tri else "", _rate(upd, ci), ", ".join(read_as)))
        lines.append("")
    lines.append("## Calls on unaffected chromosomes")
    lines.append("")
    n_runs = len({(r["replicate"], r["depth"], r["f"]) for r in rows})
    if fps:
        lines.append("%d call(s) in %d of %d runs (false_positives.tsv):" % (len(fps), len({(r["replicate"], r["depth"], r["f"]) for r in fps}), n_runs))
        lines.append("")
        for r in fps:
            lines.append("- %g×, f = %.2f: %s %s %.1f–%.1f Mb, f = %.2f (%s)" % (r["depth"], r["f"], r["chrom"], r["type"], r["start"] / 1e6, r["end"] / 1e6, r["f_called"], r["source"]))
    else:
        lines.append("None in %d runs." % n_runs)
    lines.append("")
    lines.append("## Cell-fraction error of detected events (median |estimate − f|)")
    lines.append("")
    lines.append("| type | depth | f_called | f_lrr | f_baf | f_phase | parent of origin correct |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for depth in depths:
        for kind, _, _ in TYPES:
            sel = [r for r in rows if r["depth"] == depth and r["type"] == kind and r["detected"] and r["size_mb"] != "whole"]
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
            rows = [r for r in rows if r["size_mb"] != "whole"]
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
