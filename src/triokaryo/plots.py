"""Figures for review. No legend and no title in the image: each figure writes a sidecar (title, caption, key) and a
legend image of its own (PNG, SVG, PDF), so the figure can go into a manuscript as it is. Okabe-Ito colours."""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PAL = dict(child="#D55E00", father="#E69F00", mother="#CC79A7", gain="#D55E00", loss="#0072B2", loh="#CC79A7", depth="#444444",
           baf="#9C9C9C", ref="#7F7F7F", ext="#000000", mat="#CC79A7", pat="#E69F00", band="#F2F2F2")
KEY = {
    "depth": (PAL["depth"], "point", "LRR: log2 of the bin's median depth over the member's autosomal median (GC-corrected where a GC track was given)"),
    "baf": (PAL["baf"], "point", "BAF: the alt-allele fraction at the member's heterozygous sites (a sample of them)"),
    "gain": (PAL["gain"], "line", "a called gain (line at the segment's mean LRR)"),
    "loss": (PAL["loss"], "line", "a called loss"),
    "loh": (PAL["loh"], "line", "a called copy-neutral loss of heterozygosity (drawn at zero)"),
    "mat": (PAL["mat"], "point", "an informative site whose alt allele came from the mother (father 0/0, mother 1/1)"),
    "pat": (PAL["pat"], "point", "an informative site whose alt allele came from the father (father 1/1, mother 0/0)"),
    "ext": (PAL["ext"], "bracket", "an event given from elsewhere (--events), drawn above the track"),
    "ref": (PAL["ref"], "line", "reference lines: zero LRR, BAF at 1/2, 1/3 and 2/3; the centromere"),
    "trio": (PAL["child"], "point", "the child over the parents' mean depth, site by site, per bin (the within-family difference)"),
}
MARK = {"point": "filled circle", "line": "line", "bracket": "bracket"}
MM = 1 / 25.4
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                     "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 300})


def _save(fig, out_dir, name, title, caption, keys):
    os.makedirs(os.path.join(out_dir, "legends"), exist_ok=True)
    paths = {}
    for ext in ("png", "svg", "pdf"):
        p = os.path.join(out_dir, "%s.%s" % (name, ext))
        fig.savefig(p, bbox_inches="tight")
        paths[ext] = p
    plt.close(fig)
    with open(os.path.join(out_dir, name + ".txt"), "w") as fh:
        fh.write("title: %s\ncaption: %s\nkey:\n" % (title, caption))
        for k in keys:
            col, mk, text = KEY[k]
            fh.write("  %s  %s  %s\n" % (col, MARK[mk], text))
    # the legend as its own image
    fig, ax = plt.subplots(figsize=(120 * MM, (6 + 5 * len(keys)) * MM))
    ax.axis("off")
    for i, k in enumerate(keys):
        col, mk, text = KEY[k]
        y = 1 - (i + 0.5) / len(keys)
        if mk == "point":
            ax.plot([0.03], [y], "o", color=col, ms=5, transform=ax.transAxes)
        elif mk == "line":
            ax.plot([0.01, 0.05], [y, y], color=col, lw=2, transform=ax.transAxes)
        else:
            ax.plot([0.01, 0.05], [y, y], color=col, lw=1.2, transform=ax.transAxes)
            ax.plot([0.01, 0.01], [y - 0.03, y], color=col, lw=1.2, transform=ax.transAxes)
            ax.plot([0.05, 0.05], [y - 0.03, y], color=col, lw=1.2, transform=ax.transAxes)
        ax.text(0.08, y, text, va="center", fontsize=6, transform=ax.transAxes)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(out_dir, "legends", "%s_legend.%s" % (name, ext)), bbox_inches="tight", transparent=True)
    plt.close(fig)
    return dict(name=name, title=title, caption=caption, keys=keys, **paths)


def _genome_axis(bins, genome):
    """Cumulative coordinates: offset per chromosome, in the order of the bins."""
    off, x = {}, 0
    for c in bins.index:
        off[c] = x
        x += genome.length[c]
    return off, x


def _draw_events(ax, evs, off, y_of, lw=2.0):
    for e in evs:
        col = PAL["gain" if e.type == "gain" else "loss" if e.type == "loss" else "loh"]
        x0 = off.get(e.chrom, 0) + e.start
        x1 = off.get(e.chrom, 0) + e.end
        ax.plot([x0, x1], [y_of(e), y_of(e)], color=col, lw=lw, solid_capstyle="butt")


def _draw_external(ax, ext, off, y, sample):
    for x in ext:
        if x.sample != sample or x.chrom not in off:
            continue
        x0, x1 = off[x.chrom] + x.start, off[x.chrom] + x.end
        ax.plot([x0, x1], [y, y], color=PAL["ext"], lw=0.9)
        ax.plot([x0, x0], [y - 0.08, y], color=PAL["ext"], lw=0.9)
        ax.plot([x1, x1], [y - 0.08, y], color=PAL["ext"], lw=0.9)


def fig_genome(trio, bins, scan, events, external, genome, out_dir, name="genome", max_points=30000, min_dp=8, min_gq=20):
    """Every chromosome: LRR and BAF per member, the child over the parents' mean, the calls, the events given from elsewhere."""
    off, total = _genome_axis(bins, genome)
    roles = ("child", "father", "mother")
    fig, axes = plt.subplots(7, 1, figsize=(180 * MM, 150 * MM), sharex=True, gridspec_kw=dict(height_ratios=[1, 0.8, 1, 0.8, 1, 0.8, 1], hspace=0.12))
    xs = np.array([off[c] for c in bins.chrom]) + (bins.start + bins.end) / 2.0
    rng = np.random.default_rng(1)
    for m, role in enumerate(roles):
        ax = axes[2 * m]
        for k, c in enumerate(bins.index):
            if k % 2:
                ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
        ax.axhline(0, color=PAL["ref"], lw=0.5, ls=":")
        ax.scatter(xs, bins.lrr_gc[m], s=1.2, c=PAL["depth"], lw=0, rasterized=True)
        _draw_events(ax, [e for e in events if e.role == role and e.type in ("gain", "loss")], off, lambda e: e.lrr)
        _draw_events(ax, [e for e in events if e.role == role and e.type == "LOH"], off, lambda e: 0.0)
        _draw_external(ax, external, off, 1.05, trio.members[m])
        ax.set_ylim(-1.5, 1.2)
        ax.set_ylabel("%s\nLRR" % role, fontsize=6)
        ax.text(-0.06, 1.0, "abcdefg"[2 * m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        ax = axes[2 * m + 1]
        for k, c in enumerate(bins.index):
            if k % 2:
                ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
        # the heterozygous sites, a sample
        px, py = [], []
        for c, sites in scan.chroms.items():
            if c not in off:
                continue
            het = (sites.gt[m] == 1) & (sites.dp[m] >= min_dp) & (sites.gq[m] >= min_gq)
            idx = np.flatnonzero(het)
            if len(idx) > max_points // max(1, len(scan.chroms)):
                idx = rng.choice(idx, max_points // max(1, len(scan.chroms)), replace=False)
            px.append(off[c] + sites.pos[idx])
            py.append(sites.alt[m][idx] / np.maximum(sites.dp[m][idx], 1))
        if px:
            ax.scatter(np.concatenate(px), np.concatenate(py), s=0.5, c=PAL["baf"], lw=0, rasterized=True)
        for yv in (1 / 3, 0.5, 2 / 3):
            ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("BAF", fontsize=6)
        ax.text(-0.06, 1.0, "abcdefg"[2 * m + 1], transform=ax.transAxes, fontweight="bold", fontsize=8)
    ax = axes[6]
    for k, c in enumerate(bins.index):
        if k % 2:
            ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
    ax.axhline(0, color=PAL["ref"], lw=0.5, ls=":")
    ax.scatter(xs, bins.child_vs_mid, s=1.2, c=PAL["child"], lw=0, rasterized=True)
    ax.set_ylim(-1.5, 1.2)
    ax.set_ylabel("child over\nparents' mean", fontsize=6)
    ax.text(-0.06, 1.0, "g", transform=ax.transAxes, fontweight="bold", fontsize=8)
    ax.set_xlim(0, total)
    ax.set_xticks([off[c] + genome.length[c] / 2 for c in bins.index])
    ax.set_xticklabels([c[3:] for c in bins.index], fontsize=5)
    ax.set_xlabel("chromosome")
    caption = ("For the child (a, b), the father (c, d) and the mother (e, f): the LRR of every %d-kb bin (log2 of its median depth over the member's autosomal "
               "median%s) with the called gains, losses and copy-neutral losses of heterozygosity as lines, and the B-allele fraction at the member's "
               "heterozygous sites (a sample; dotted lines at 1/2, 1/3 and 2/3). (g) The child's depth over the parents' mean, site by site, per bin: "
               "the within-family difference, zero where the child inherited what the parents carry. Events given from elsewhere are brackets above "
               "the LRR." % (bins.bin_size // 1000, ", GC-corrected" if np.isfinite(bins.gc).any() else ""))
    return _save(fig, out_dir, name, "Large chromosomal events in trio %s, genome-wide" % trio.name, caption,
                 ["depth", "baf", "gain", "loss", "loh", "trio", "ext", "ref"])


def fig_chrom(trio, chrom, bins, scan, events, external, genome, out_dir, min_dp=8, min_gq=20, max_points=6000):
    """One chromosome: for each member the LRR with its calls, the BAF with the informative sites coloured by the parent of the
    alt allele (the child), and the heterozygosity rate."""
    sl = bins.of(chrom)
    sites = scan.sites(chrom)
    L = genome.length[chrom]
    pe = genome.p_end.get(chrom, 0)
    roles = ("child", "father", "mother")
    fig, axes = plt.subplots(3, 3, figsize=(180 * MM, 120 * MM), sharex=True, gridspec_kw=dict(height_ratios=[1, 1, 0.6], hspace=0.15, wspace=0.25))
    x = (bins.start[sl] + bins.end[sl]) / 2e6
    rng = np.random.default_rng(2)
    # the informative sites of the child: the parent of the alt allele
    inf_mat = inf_pat = None
    if sites is not None:
        from .trio import _confident
        ok = (sites.dp[0] >= min_dp) & (sites.gq[0] >= min_gq)
        for mm in (1, 2):
            ok &= _confident(sites.dp[mm], sites.gq[mm], sites.gt[mm], min_dp, min_gq)
        inf_mat = ok & (sites.gt[1] == 0) & (sites.gt[2] == 2)
        inf_pat = ok & (sites.gt[1] == 2) & (sites.gt[2] == 0)
    for m, role in enumerate(roles):
        ax = axes[0, m]
        ax.axhline(0, color=PAL["ref"], lw=0.5, ls=":")
        ax.axvline(pe / 1e6, color=PAL["ref"], lw=0.5, ls="--")
        ax.scatter(x, bins.lrr_gc[m][sl], s=4, c=PAL["depth"], lw=0)
        evs = [e for e in events if e.role == role and e.chrom == chrom]
        for e in evs:
            col = PAL["gain" if e.type == "gain" else "loss" if e.type == "loss" else "loh"]
            yv = e.lrr if e.type in ("gain", "loss") else 0.0
            ax.plot([e.start / 1e6, e.end / 1e6], [yv, yv], color=col, lw=2.2, solid_capstyle="butt")
        for xe in [v for v in external if v.sample == trio.members[m] and v.chrom == chrom]:
            ax.plot([xe.start / 1e6, xe.end / 1e6], [1.05, 1.05], color=PAL["ext"], lw=0.9)
            ax.plot([xe.start / 1e6] * 2, [0.97, 1.05], color=PAL["ext"], lw=0.9)
            ax.plot([xe.end / 1e6] * 2, [0.97, 1.05], color=PAL["ext"], lw=0.9)
        ax.set_ylim(-1.5, 1.2)
        ax.set_ylabel("LRR" if m == 0 else "")
        ax.set_title(role, fontsize=7, loc="left")
        ax.text(-0.18, 1.02, "abc"[m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        ax = axes[1, m]
        if sites is not None:
            het = (sites.gt[m] == 1) & (sites.dp[m] >= min_dp) & (sites.gq[m] >= min_gq)
            idx = np.flatnonzero(het)
            if len(idx) > max_points:
                idx = rng.choice(idx, max_points, replace=False)
            ax.scatter(sites.pos[idx] / 1e6, sites.alt[m][idx] / np.maximum(sites.dp[m][idx], 1), s=1.5, c=PAL["baf"], lw=0, rasterized=True)
            if m == 0 and inf_mat is not None:
                for mask, key in ((inf_mat, "mat"), (inf_pat, "pat")):
                    j = np.flatnonzero(mask)
                    if len(j) > max_points // 2:
                        j = rng.choice(j, max_points // 2, replace=False)
                    ax.scatter(sites.pos[j] / 1e6, sites.alt[0][j] / np.maximum(sites.dp[0][j], 1), s=2.5, c=PAL[key], lw=0, rasterized=True)
        for yv in (1 / 3, 0.5, 2 / 3):
            ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
        ax.axvline(pe / 1e6, color=PAL["ref"], lw=0.5, ls="--")
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("BAF" if m == 0 else "")
        ax.text(-0.18, 1.02, "def"[m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        ax = axes[2, m]
        hr = bins.het_rate[m][sl]
        ax.scatter(x, hr, s=4, c=PAL["depth"], lw=0)
        ax.axvline(pe / 1e6, color=PAL["ref"], lw=0.5, ls="--")
        ax.set_ylim(0, max(0.05, float(np.nanmax(hr)) * 1.2 if np.isfinite(hr).any() else 0.5))
        ax.set_ylabel("het rate" if m == 0 else "")
        ax.set_xlabel("%s (Mb)" % chrom)
        ax.text(-0.18, 1.02, "ghi"[m], transform=ax.transAxes, fontweight="bold", fontsize=8)
    axes[0, 0].set_xlim(0, L / 1e6)
    evs = [e for e in events if e.chrom == chrom]
    what = "; ".join("%s: %s %s f %.2f%s" % (e.role, e.type, e.span, e.f if np.isfinite(e.f) else float("nan"), (", " + e.origin) if e.origin else "") for e in evs) or "no event called"
    caption = ("%s in trio %s. Top: the LRR of each %d-kb bin per member with the calls as lines (gain, loss, or copy-neutral LOH at zero); events given "
               "from elsewhere as brackets. Middle: the B-allele fraction at each member's heterozygous sites; in the child, the sites where the parents "
               "are opposite homozygotes are coloured by the parent of the alt allele - under a gain the duplicated parent's allele sits at 2/3, under a "
               "loss the retained parent's at 1, under a copy-neutral LOH the retained parent's at (1 + f)/2. Bottom: the heterozygosity rate per bin, "
               "which falls to zero under a loss of heterozygosity in every cell. Dashed: the centromere. Calls: %s." % (chrom, trio.name, bins.bin_size // 1000, what))
    return _save(fig, out_dir, "chrom_%s" % chrom, "Trio %s, %s" % (trio.name, chrom), caption, ["depth", "baf", "mat", "pat", "gain", "loss", "loh", "ext", "ref"])
