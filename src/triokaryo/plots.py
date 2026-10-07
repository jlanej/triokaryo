"""Figures. The image carries no legend or title; each figure writes a sidecar (title, caption, key) and a separate
legend image (PNG, SVG, PDF) so that it can be placed in a manuscript as it is. Okabe-Ito colours."""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PAL = dict(child="#D55E00", father="#E69F00", mother="#CC79A7", gain="#D55E00", loss="#0072B2", loh="#CC79A7", depth="#444444",
           baf="#9C9C9C", ref="#7F7F7F", ext="#000000", mat="#CC79A7", pat="#E69F00", band="#F2F2F2", phased="#009E73", psite="#C8C8C8", step="#000000",
           trans="#D55E00", untrans="#56B4E9")
KEY = {
    "depth": (PAL["depth"], "point", "LRR per bin: log2 of the bin's trimmed-mean depth over the member's autosomal median (GC- and panel-corrected where applicable)"),
    "baf": (PAL["baf"], "point", "BAF: alt-allele read fraction at the member's heterozygous sites (subsampled)"),
    "gain": (PAL["gain"], "line", "called gain (line at the segment's mean LRR)"),
    "loss": (PAL["loss"], "line", "called loss"),
    "loh": (PAL["loh"], "line", "called copy-neutral event: loss of heterozygosity or uniparental heterodisomy (drawn at LRR 0)"),
    "mat": (PAL["mat"], "point", "informative site whose alt allele is maternal (father 0/0, mother 1/1)"),
    "pat": (PAL["pat"], "point", "informative site whose alt allele is paternal (father 1/1, mother 0/0)"),
    "ext": (PAL["ext"], "bracket", "event supplied with --events (bracket above the track)"),
    "ref": (PAL["ref"], "line", "reference lines: LRR 0; BAF 1/2, 1/3 and 2/3; the centromere"),
    "trio": (PAL["child"], "point", "child's depth over the parents' mean depth per site, median per bin (within-family difference)"),
    "psite": (PAL["psite"], "point", "phased fraction at one site: maternal allele (child) or transmitted allele (parent), reference bias removed"),
    "phased": (PAL["phased"], "point", "phased fraction pooled over a window of sites (depth-weighted): 1/2 with equal homologue copy, 1/2 +/- d along an event with the sign "
                                       "giving the parent (child: below 1/2 the paternal homologue is in excess, above 1/2 the maternal; parent: below 1/2 the event lies on the "
                                       "untransmitted homologue, above 1/2 on the transmitted one); hollow where the window deviates in two or more members (paralogy; excluded from "
                                       "fits and scan)"),
    "step": (PAL["step"], "line", "step fit of the track (exact total-variation denoising)"),
    "aux_mat": (PAL["mat"], "line", "child's maternal fraction at sites where the father is homozygous and the mother heterozygous; departs from the main track where the child "
                                    "carries two different maternal homologues (maternal meiotic trisomy or heterodisomy) and rejoins it after a crossover"),
    "aux_pat": (PAL["pat"], "line", "the same at sites where the mother is homozygous and the father heterozygous; departs where the child carries two different paternal homologues"),
    "aux_par": (PAL["untrans"], "line", "parent's transmitted-allele fraction at sites where the child is heterozygous and the other parent homozygous; departs from the main "
                                        "track where the child carries two different homologues of this parent"),
    "cmat": (PAL["mat"], "line", "child's maternal copy number: the LRR step fit's copy number (2 x 2^LRR) times the maternal fraction's step fit"),
    "cpat": (PAL["pat"], "line", "child's paternal copy number, likewise"),
    "ctrans": (PAL["trans"], "line", "parent's copy number of the transmitted homologue, likewise"),
    "cuntrans": (PAL["untrans"], "line", "parent's copy number of the untransmitted homologue"),
    "xover": (PAL["step"], "line", "a crossover: where the auxiliary track's state (two different homologues, or one) changes along a whole-chromosome event (dashed vertical line)"),
}
MARK = {"point": "filled circle", "line": "line", "bracket": "bracket"}
GENOME_KEYS = ["depth", "step", "baf", "phased", "gain", "loss", "loh", "trio", "cmat", "cpat", "ext", "ref"]
CHROM_KEYS = ["depth", "step", "baf", "mat", "pat", "psite", "phased", "aux_mat", "aux_pat", "aux_par", "xover", "cmat", "cpat", "ctrans", "cuntrans", "gain", "loss", "loh", "ext", "ref"]
LANDSCAPE_KEYS = ["gain", "loss", "loh", "ext"]
DIRECTION = ("Sign convention: in the child, a phased fraction below 1/2 means the paternal homologue is in excess (paternal gain or maternal loss) and above 1/2 the "
             "maternal homologue; in a parent, below 1/2 means the event lies on the untransmitted homologue and above 1/2 on the transmitted one.")
MM = 1 / 25.4
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
                     "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 300})


def write_sidecar(out_dir, name, title, caption, keys):
    """The figure's title, caption and key (colour, mark, text), and the key names, beside the image."""
    with open(os.path.join(out_dir, name + ".txt"), "w") as fh:
        fh.write("title: %s\ncaption: %s\nkeys: %s\nkey:\n" % (title, caption, ",".join(keys)))
        for k in keys:
            col, mk, text = KEY[k]
            fh.write("  %s  %s  %s\n" % (col, MARK[mk], text))


def read_sidecar(path):
    """(title, caption, key names) back from a sidecar; the names from its 'keys:' line, else by the figure's kind."""
    title = caption = ""
    keys = None
    for line in open(path):
        if line.startswith("title: "):
            title = line[7:].rstrip("\n")
        elif line.startswith("caption: "):
            caption = line[9:].rstrip("\n")
        elif line.startswith("keys: "):
            keys = [k for k in line[6:].strip().split(",") if k]
    if keys is None:
        base = os.path.basename(path)[:-4]
        keys = GENOME_KEYS if base == "genome" else LANDSCAPE_KEYS if base == "landscape" else CHROM_KEYS
    return title, caption, keys


def _save(fig, out_dir, name, title, caption, keys):
    os.makedirs(os.path.join(out_dir, "legends"), exist_ok=True)
    paths = {}
    for ext in ("png", "svg", "pdf"):
        p = os.path.join(out_dir, "%s.%s" % (name, ext))
        fig.savefig(p, bbox_inches="tight")
        paths[ext] = p
    plt.close(fig)
    write_sidecar(out_dir, name, title, caption, keys)
    write_legend(out_dir, name, keys)
    return dict(name=name, title=title, caption=caption, keys=keys, **paths)


def write_legend(out_dir, name, keys):
    """The legend as its own image (PNG, SVG, PDF) under legends/."""
    os.makedirs(os.path.join(out_dir, "legends"), exist_ok=True)
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
        ax.text(0.08, y, text, va="center", fontsize=6, transform=ax.transAxes, wrap=True)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(out_dir, "legends", "%s_legend.%s" % (name, ext)), bbox_inches="tight", transparent=True)
    plt.close(fig)


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


def _lrr_step(bins, m):
    from .phase import lrr_step
    return lrr_step(bins, m)


def _phased_points(tracks, m, chrom=None, off=None):
    """(x, frac) of the member's pooled windows, genome-wide (off given) or for one chromosome, and (x, step)."""
    xs, ys, st = [], [], []
    if not tracks:
        return np.array([]), np.array([]), np.array([])
    for c, t in tracks[m].items():
        if chrom and c != chrom:
            continue
        if off is not None and c not in off:
            continue
        x = (off[c] + t.w_mid) if off is not None else t.w_mid / 1e6
        xs.append(x)
        ys.append(t.w_frac)
        st.append(t.w_step)
    if not xs:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(xs), np.concatenate(ys), np.concatenate(st)


def _copies_lines(ax, tracks, m, chrom=None, off=None, lw=0.9):
    """Per-homologue copy number along the member: maternal and paternal (child), transmitted and untransmitted (parent)."""
    if not tracks or m >= len(tracks):
        return
    c1, c2 = (PAL["mat"], PAL["pat"]) if m == 0 else (PAL["trans"], PAL["untrans"])
    for c, t in tracks[m].items():
        if chrom and c != chrom:
            continue
        if off is not None and c not in off:
            continue
        x = (off[c] + t.w_mid) if off is not None else t.w_mid / 1e6
        ax.plot(x, t.w_copies_tag, color=c1, lw=lw, rasterized=True)
        ax.plot(x, t.w_copies_other, color=c2, lw=lw, rasterized=True)


def genome_caption(bin_size, gc_corrected):
    """The genome figure's caption; also used by `triokaryo report` to rebuild the sidecar from the tables."""
    return ("Child (a, b), father (c, d) and mother (e, f): LRR per %d-kb bin (log2 of the bin's trimmed-mean depth over the member's autosomal median%s) with its "
            "step fit and the called gains, losses and copy-neutral events as lines; and BAF at the member's heterozygous sites (subsampled; dotted lines at "
            "1/2, 1/3 and 2/3) overlaid with the phased fraction pooled by window (maternal allele along the child, transmitted allele along a parent), which "
            "is 1/2 with equal homologue copy and deviates along an event with the sign giving the parent of origin (%s). (g) The child's depth over the "
            "parents' mean per site, median per bin: zero where the child's copy number equals the parental mean. (h) The child's maternal and paternal copy "
            "number: the LRR-derived copy number split by the phased fraction. Brackets above the LRR: events supplied with --events."
            % (bin_size // 1000, ", GC-corrected" if gc_corrected else "", DIRECTION))


def chrom_caption(chrom, trio_name, bin_size, events):
    """A chromosome figure's caption, with the chromosome's calls summarised; also used by `triokaryo report`."""
    evs = [e for e in events if e.chrom == chrom]
    what = "; ".join("%s: %s %s f %.2f%s" % (e.role, e.type, e.span, e.f if np.isfinite(e.f) else float("nan"), (", " + (e.origin_phase or e.origin)) if (e.origin_phase or e.origin) else "")
                     for e in evs) or "no event called"
    return ("%s in trio %s. Row 1: LRR per %d-kb bin per member with its step fit and calls as lines (gain, loss; copy-neutral events at 0); supplied events as "
            "brackets. Row 2: BAF at each member's heterozygous sites; in the child, informative sites (parents opposite homozygotes) are coloured by the "
            "parent of the alt allele. Row 3: the phased fraction (maternal allele along the child at informative sites; transmitted allele along a parent "
            "where the child is homozygous) per site (faint), pooled by window, with its step fit: 1/2 with equal homologue copy, 1/2 +/- d along an event "
            "with the sign giving the parent of origin (%s). Thin lines: the auxiliary tracks (child: sites where one parent is homozygous; parent: sites "
            "where the child is heterozygous), which depart from the main track where the child carries two different homologues of one parent; dashed "
            "vertical lines mark the crossovers where that state changes along a whole-chromosome event. Row 4: "
            "per-homologue copy number, the LRR-derived copy number split by the phased fraction (maternal and paternal in the child; transmitted and "
            "untransmitted in a parent). Row 5: heterozygosity rate per bin, zero under a constitutional loss of heterozygosity. Dashed line: centromere. "
            "Calls: %s." % (chrom, trio_name, bin_size // 1000, DIRECTION, what))


def fig_genome(trio, bins, scan, events, external, genome, out_dir, name="genome", max_points=30000, min_dp=8, min_gq=20, tracks=None):
    """Genome-wide figure: per member the LRR (with step fit and calls) and the BAF (with the pooled phased fraction); the
    child's depth over the parents' mean; the child's maternal and paternal copy number; the supplied events."""
    off, total = _genome_axis(bins, genome)
    roles = ("child", "father", "mother")
    fig, axes = plt.subplots(8, 1, figsize=(180 * MM, 165 * MM), sharex=True, gridspec_kw=dict(height_ratios=[1, 0.8, 1, 0.8, 1, 0.8, 1, 1], hspace=0.12))
    xs = np.array([off[c] for c in bins.chrom]) + (bins.start + bins.end) / 2.0
    rng = np.random.default_rng(1)
    letters = "abcdefgh"
    for m, role in enumerate(roles):
        ax = axes[2 * m]
        for k, c in enumerate(bins.index):
            if k % 2:
                ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
        ax.axhline(0, color=PAL["ref"], lw=0.5, ls=":")
        ax.scatter(xs, bins.lrr_gc[m], s=1.2, c=PAL["depth"], lw=0, rasterized=True)
        ax.plot(xs, _lrr_step(bins, m), color=PAL["step"], lw=0.6, rasterized=True)
        _draw_events(ax, [e for e in events if e.role == role and e.type in ("gain", "loss")], off, lambda e: e.lrr)
        _draw_events(ax, [e for e in events if e.role == role and e.type == "LOH"], off, lambda e: 0.0)
        _draw_external(ax, external, off, 1.05, trio.members[m])
        ax.set_ylim(-1.5, 1.2)
        ax.set_ylabel("%s\nLRR" % role, fontsize=6)
        ax.text(-0.06, 1.0, letters[2 * m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        ax = axes[2 * m + 1]
        for k, c in enumerate(bins.index):
            if k % 2:
                ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
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
        gx, gy, _ = _phased_points(tracks, m, off=off)
        if len(gx):
            ax.scatter(gx, gy, s=0.8, c=PAL["phased"], lw=0, rasterized=True)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("BAF", fontsize=6)
        ax.text(-0.06, 1.0, letters[2 * m + 1], transform=ax.transAxes, fontweight="bold", fontsize=8)
    ax = axes[6]
    for k, c in enumerate(bins.index):
        if k % 2:
            ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
    ax.axhline(0, color=PAL["ref"], lw=0.5, ls=":")
    ax.scatter(xs, bins.child_vs_mid, s=1.2, c=PAL["child"], lw=0, rasterized=True)
    ax.set_ylim(-1.5, 1.2)
    ax.set_ylabel("child over\nparents' mean", fontsize=6)
    ax.text(-0.06, 1.0, "g", transform=ax.transAxes, fontweight="bold", fontsize=8)
    ax = axes[7]
    for k, c in enumerate(bins.index):
        if k % 2:
            ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
    for yv in (0, 1, 2):
        ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
    _copies_lines(ax, tracks, 0, off=off, lw=0.7)
    ax.set_ylim(-0.15, 2.6)
    ax.set_yticks([0, 1, 2])
    ax.set_ylabel("child's copies\nmaternal, paternal", fontsize=6)
    ax.text(-0.06, 1.0, "h", transform=ax.transAxes, fontweight="bold", fontsize=8)
    ax.set_xlim(0, total)
    ax.set_xticks([off[c] + genome.length[c] / 2 for c in bins.index])
    ax.set_xticklabels([c[3:] for c in bins.index], fontsize=5)
    ax.set_xlabel("chromosome")
    caption = genome_caption(bins.bin_size, bool(np.isfinite(bins.gc).any()))
    return _save(fig, out_dir, name, "Large chromosomal events in trio %s, genome-wide" % trio.name, caption, GENOME_KEYS)


def fig_chrom(trio, chrom, bins, scan, events, external, genome, out_dir, min_dp=8, min_gq=20, max_points=6000, tracks=None):
    """One chromosome: per member the LRR with calls and step fit, the BAF (the child's informative sites coloured by the
    parent of the alt allele), the phased fraction (sites, pooled windows, step fit, auxiliary tracks), per-homologue copy
    number, and the heterozygosity rate."""
    sl = bins.of(chrom)
    sites = scan.sites(chrom)
    L = genome.length[chrom]
    pe = genome.p_end.get(chrom, 0)
    roles = ("child", "father", "mother")
    fig, axes = plt.subplots(5, 3, figsize=(180 * MM, 175 * MM), sharex=True, gridspec_kw=dict(height_ratios=[1, 1, 1, 0.8, 0.5], hspace=0.15, wspace=0.25))
    x = (bins.start[sl] + bins.end[sl]) / 2e6
    rng = np.random.default_rng(2)
    letters = "abcdefghijklmno"
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
        ax.plot(x, _lrr_step(bins, m)[sl], color=PAL["step"], lw=0.7)
        evs = [e for e in events if e.role == role and e.chrom == chrom]
        for e in evs:
            col = PAL["gain" if e.type == "gain" else "loss" if e.type == "loss" else "loh"]
            yv = e.lrr if e.type in ("gain", "loss") and e.source != "phased" else 0.0
            ax.plot([e.start / 1e6, e.end / 1e6], [yv, yv], color=col, lw=2.2, solid_capstyle="butt")
        for xe in [v for v in external if v.sample == trio.members[m] and v.chrom == chrom]:
            ax.plot([xe.start / 1e6, xe.end / 1e6], [1.05, 1.05], color=PAL["ext"], lw=0.9)
            ax.plot([xe.start / 1e6] * 2, [0.97, 1.05], color=PAL["ext"], lw=0.9)
            ax.plot([xe.end / 1e6] * 2, [0.97, 1.05], color=PAL["ext"], lw=0.9)
        ax.set_ylim(-1.5, 1.2)
        ax.set_ylabel("LRR" if m == 0 else "")
        ax.set_title(role, fontsize=7, loc="left")
        ax.text(-0.18, 1.02, letters[m], transform=ax.transAxes, fontweight="bold", fontsize=8)
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
        ax.text(-0.18, 1.02, letters[3 + m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        # the phased fraction: every site faint, the pooled windows, the step fit
        ax = axes[2, m]
        t = tracks[m].get(chrom) if tracks and m < len(tracks) else None
        if t is not None and sites is not None:
            j = np.arange(len(t.idx))
            if len(j) > max_points:
                j = rng.choice(j, max_points, replace=False)
            ax.scatter(sites.pos[t.idx[j]] / 1e6, np.clip(t.frac[j], 0, 1), s=1.2, c=PAL["psite"], lw=0, rasterized=True)
            for name, col in (("father_hom", PAL["mat"]), ("mother_hom", PAL["pat"]), ("child_het", PAL["untrans"])):
                if name in t.aux:
                    af, an = t.aux[name]
                    ax.plot(t.w_mid / 1e6, np.where(an >= 5, af, np.nan), color=col, lw=0.7, alpha=0.9)
            sh = t.w_shared if t.w_shared is not None else np.zeros(len(t.w_mid), dtype=bool)
            ax.scatter(t.w_mid[~sh] / 1e6, t.w_frac[~sh], s=3, c=PAL["phased"], lw=0, rasterized=True)
            if sh.any():
                ax.scatter(t.w_mid[sh] / 1e6, t.w_frac[sh], s=4, facecolors="none", edgecolors=PAL["phased"], lw=0.4, rasterized=True)
            ax.step(t.w_mid / 1e6, t.w_step, where="mid", color=PAL["step"], lw=0.8)
        for e in [e for e in events if e.role == role and e.chrom == chrom and e.crossovers]:
            for xo in e.crossovers.split(";"):
                if xo:
                    ax.axvline(float(xo), color=PAL["step"], lw=0.7, ls="--")               # a crossover placed by the auxiliary track
        for yv in (1 / 3, 0.5, 2 / 3):
            ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
        ax.axvline(pe / 1e6, color=PAL["ref"], lw=0.5, ls="--")
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel("maternal allele\nfraction" if m == 0 else "transmitted allele\nfraction", fontsize=6)
        ax.text(-0.18, 1.02, letters[6 + m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        # the two homologues' copies
        ax = axes[3, m]
        for yv in (0, 1, 2):
            ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
        ax.axvline(pe / 1e6, color=PAL["ref"], lw=0.5, ls="--")
        _copies_lines(ax, tracks, m, chrom=chrom)
        ax.set_ylim(-0.15, 2.6)
        ax.set_yticks([0, 1, 2])
        ax.set_ylabel("copies\nmaternal, paternal" if m == 0 else "copies\npassed, not passed", fontsize=6)
        ax.text(-0.18, 1.02, letters[9 + m], transform=ax.transAxes, fontweight="bold", fontsize=8)
        ax = axes[4, m]
        hr = bins.het_rate[m][sl]
        ax.scatter(x, hr, s=4, c=PAL["depth"], lw=0)
        ax.axvline(pe / 1e6, color=PAL["ref"], lw=0.5, ls="--")
        ax.set_ylim(0, max(0.05, float(np.nanmax(hr)) * 1.2 if np.isfinite(hr).any() else 0.5))
        ax.set_ylabel("het rate" if m == 0 else "")
        ax.set_xlabel("%s (Mb)" % chrom)
        ax.text(-0.18, 1.02, letters[12 + m], transform=ax.transAxes, fontweight="bold", fontsize=8)
    axes[0, 0].set_xlim(0, L / 1e6)
    caption = chrom_caption(chrom, trio.name, bins.bin_size, events)
    return _save(fig, out_dir, "chrom_%s" % chrom, "Trio %s, %s" % (trio.name, chrom), caption, CHROM_KEYS)


def fig_landscape(rows, genome, out_dir, name="landscape", max_labels=60):
    """The cohort's events on one genome axis: (a) events per chromosome by type; (b) one row per trio, each event a bar
    coloured by type (the child's thick, a parent's thin), matched events marked. rows: [(trio, Event)]."""
    trios = sorted({t for t, _ in rows})
    chroms = [c for c in genome.chroms if c != "chrY"]
    off, x = {}, 0
    for c in chroms:
        off[c] = x
        x += genome.length[c]
    total = x
    n = max(1, len(trios))
    fig, axes = plt.subplots(2, 1, figsize=(180 * MM, (45 + min(6 * n, 180)) * MM), sharex=True,
                             gridspec_kw=dict(height_ratios=[1, max(1.2, min(n, 30) / 6.0)], hspace=0.08))
    ax = axes[0]
    for k, c in enumerate(chroms):
        if k % 2:
            for a in axes:
                a.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
    width = 0.8
    for i, c in enumerate(chroms):
        xm = off[c] + genome.length[c] / 2
        bottom = 0
        for kind, col in (("gain", PAL["gain"]), ("loss", PAL["loss"]), ("LOH", PAL["loh"]), ("UPD", PAL["loh"])):
            k = sum(1 for _, e in rows if e.chrom == c and e.type == kind)
            if k:
                ax.bar(xm, k, width=width * genome.length[c], bottom=bottom, color=col, lw=0)
                bottom += k
    ax.set_ylabel("events", fontsize=6)
    ax.text(-0.06, 1.0, "a", transform=ax.transAxes, fontweight="bold", fontsize=8)
    ax = axes[1]
    yof = {t: i for i, t in enumerate(trios)}
    for t, e in rows:
        if e.chrom not in off:
            continue
        col = PAL["gain" if e.type == "gain" else "loss" if e.type == "loss" else "loh"]
        y = yof[t] + {"child": 0.0, "father": 0.28, "mother": -0.28}.get(e.role, 0.0)
        ax.plot([off[e.chrom] + e.start, off[e.chrom] + e.end], [y, y], color=col, lw=2.6 if e.role == "child" else 1.2, solid_capstyle="butt", alpha=1.0 if e.role == "child" else 0.75)
        if e.external:
            ax.plot([off[e.chrom] + e.start, off[e.chrom] + e.end], [y + 0.14, y + 0.14], color=PAL["ext"], lw=0.5)
    ax.set_ylim(-0.8, n - 0.2)
    if n <= max_labels:
        ax.set_yticks(range(n))
        ax.set_yticklabels(trios, fontsize=5)
    else:
        ax.set_yticks([])
        ax.set_ylabel("%d trios" % n, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlim(0, total)
    ax.set_xticks([off[c] + genome.length[c] / 2 for c in chroms])
    ax.set_xticklabels([c[3:] for c in chroms], fontsize=5)
    ax.set_xlabel("chromosome")
    ax.text(-0.06, 1.0, "b", transform=ax.transAxes, fontweight="bold", fontsize=8)
    caption = ("(a) Events per chromosome, stacked by type (gain, loss, copy-neutral LOH or heterodisomy). (b) One row per trio: each event a bar over its span, "
               "coloured by type; the child's thick on the row's centre line, the father's thin above it, the mother's thin below; a thin black line over a bar "
               "marks an event matched by a supplied event. %d trios, %d events." % (n, len(rows)))
    return _save(fig, out_dir, name, "Large chromosomal events in the cohort", caption, LANDSCAPE_KEYS)

