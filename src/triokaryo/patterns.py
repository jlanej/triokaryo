"""Pattern cards for the guide: the expected appearance of each event type in the five rows of a chromosome figure,
drawn from idealised tracks with the noise of a 30x genome (no data: every point is simulated from the expected values)."""
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .model import f_from_d  # noqa: E402,F401  (the formulas the captions quote)
from .plots import MM, PAL, write_legend, write_sidecar  # noqa: E402

PATTERN_KEYS = ["depth", "step", "baf", "psite", "phased", "aux_mat", "aux_pat", "cmat", "cpat", "gain", "loss", "loh", "ref"]
L = 100.0   # Mb


def _const(v):
    return lambda x: np.full(len(x), float(v))


def _piece(inside, v_in, v_out=0.5):
    lo, hi = inside
    return lambda x: np.where((x >= lo) & (x < hi), float(v_in), float(v_out))


def _cards():
    """Each card: title, the event's span, the called type, and the rows as functions of position (Mb):
    lrr; bands (a list of band centres for the raw BAF, or a function returning None where there are no heterozygous
    calls); main (the maternal fraction's expectation, None where the main track has no sites); aux_f / aux_m (the tracks
    read at the father's / the mother's homozygous sites); copies (maternal, paternal); het (the rate relative to normal)."""
    ev = (30.0, 70.0)
    f1, f2, f3 = 0.40, 0.60, 0.40
    g = lambda f: (1 + f) / (2 + f)   # noqa: E731
    cards = [
        dict(title="gain, paternal homologue\n40% of cells", span=ev, kind="gain", lrr=_piece(ev, np.log2(1 + f1 / 2), 0.0),
             bands=lambda x: [1 - g(f1), g(f1)] if ev[0] <= x < ev[1] else [0.5], main=_piece(ev, 1 - g(f1)), aux_f=_piece(ev, 1 - g(f1)), aux_m=_piece(ev, 1 - g(f1)),
             copies=(_const(1.0), _piece(ev, 1 + f1, 1.0)), het=_const(1.0)),
        dict(title="loss, paternal copy\n60% of cells", span=ev, kind="loss", lrr=_piece(ev, np.log2(1 - f2 / 2), 0.0),
             bands=lambda x: [(1 - f2) / (2 - f2), 1 / (2 - f2)] if ev[0] <= x < ev[1] else [0.5], main=_piece(ev, 1 / (2 - f2)), aux_f=_piece(ev, 1 / (2 - f2)), aux_m=_piece(ev, 1 / (2 - f2)),
             copies=(_const(1.0), _piece(ev, 1 - f2, 1.0)), het=_const(1.0)),
        dict(title="copy-neutral LOH\nmaternal copy retained, 40%", span=ev, kind="LOH", lrr=_const(0.0),
             bands=lambda x: [(1 - f3) / 2, (1 + f3) / 2] if ev[0] <= x < ev[1] else [0.5], main=_piece(ev, (1 + f3) / 2), aux_f=_piece(ev, (1 + f3) / 2), aux_m=_piece(ev, (1 + f3) / 2),
             copies=(_piece(ev, 1 + f3, 1.0), _piece(ev, 1 - f3, 1.0)), het=_const(1.0)),
        dict(title="maternal isodisomy\nwhole chromosome", span=(0.0, L), kind="LOH", lrr=_const(0.0),
             bands=lambda x: None, main=_const(1.0), aux_f=None, aux_m=None, copies=(_const(2.0), _const(0.0)), het=_const(0.0)),
        dict(title="maternal heterodisomy\nwhole chromosome", span=(0.0, L), kind="UPD", lrr=_const(0.0),
             bands=lambda x: [0.5], main=_const(1.0), aux_f=_const(0.5), aux_m=None, copies=(_const(2.0), _const(0.0)), het=_const(1.0)),
        dict(title="run of homozygosity\n(identical by descent)", span=ev, kind="LOH", lrr=_const(0.0),
             bands=lambda x: None if ev[0] <= x < ev[1] else [0.5], main=lambda x: np.where((x >= ev[0]) & (x < ev[1]), np.nan, 0.5), aux_f=None, aux_m=None,
             copies=(lambda x: np.where((x >= ev[0]) & (x < ev[1]), np.nan, 1.0), lambda x: np.where((x >= ev[0]) & (x < ev[1]), np.nan, 1.0)), het=_piece(ev, 0.0, 1.0)),
        dict(title="maternal meiotic trisomy\na crossover at 60 Mb", span=(0.0, L), kind="gain", lrr=_const(np.log2(1.5)),
             bands=lambda x: [1 / 3, 2 / 3], main=_const(2 / 3), aux_f=lambda x: np.where(x < 60, 1 / 3, 2 / 3), aux_m=_const(2 / 3),
             copies=(_const(2.0), _const(1.0)), het=lambda x: np.where(x < 60, 1.3, 1.0)),
        dict(title="maternal mitotic\n(or meiosis II) trisomy", span=(0.0, L), kind="gain", lrr=_const(np.log2(1.5)),
             bands=lambda x: [1 / 3, 2 / 3], main=_const(2 / 3), aux_f=_const(2 / 3), aux_m=_const(2 / 3),
             copies=(_const(2.0), _const(1.0)), het=_const(1.0)),
        # the sex chromosomes: the X of a son (baseline one maternal X: LRR -1, no heterozygous sites) or of a daughter (baseline 0)
        dict(title="47,XXY, maternal (a son)\nmeiosis I", span=(0.0, L), kind="gain", lrr=_const(0.0),
             bands=lambda x: [0.5], main=_const(1.0), aux_f=_const(0.5), aux_m=None,
             copies=(_const(2.0), _const(0.0)), het=_const(1.0)),
        dict(title="47,XXY, paternal (a son)\nX and Y transmitted together", span=(0.0, L), kind="gain", lrr=_const(0.0),
             bands=lambda x: [0.5], main=_const(0.5), aux_f=_const(0.5), aux_m=None,
             copies=(_const(1.0), _const(1.0)), het=_const(1.0)),
        dict(title="mosaic 46,XY/47,XXY (a son)\nextra X paternal, 40% of cells", span=(0.0, L), kind="gain", lrr=_const(np.log2(1.4 / 2)),
             bands=lambda x: [0.4 / 1.4, 1 / 1.4], main=_const(1 / 1.4), aux_f=_const(1 / 1.4), aux_m=None,
             copies=(_const(1.0), _const(0.4)), het=_const(0.8)),
        dict(title="45,X (a daughter)\npaternal X lost", span=(0.0, L), kind="loss", lrr=_const(-1.0),
             bands=lambda x: None, main=_const(1.0), aux_f=None, aux_m=None,
             copies=(_const(1.0), _const(0.0)), het=_const(0.0)),
    ]
    return cards


def _draw_card(axes, card, rng, dp=30):
    xb = np.arange(0.5, L, 1.0)
    # LRR
    ax = axes[0]
    ax.axhline(0, color=PAL["ref"], lw=0.5, ls=":")
    lrr = card["lrr"](xb)
    ax.scatter(xb, lrr + rng.normal(0, 0.03, len(xb)), s=3, c=PAL["depth"], lw=0)
    ax.plot(xb, lrr, color=PAL["step"], lw=0.7)
    col = PAL["gain" if card["kind"] == "gain" else "loss" if card["kind"] == "loss" else "loh"]
    lo, hi = card["span"]
    yv = float(card["lrr"](np.array([(lo + hi) / 2]))[0]) if card["kind"] in ("gain", "loss") else 0.0
    if card["title"].startswith("run of"):                       # a run of homozygosity is not drawn as a call here
        yv = None
    if yv is not None:
        ax.plot([lo, hi], [yv, yv], color=col, lw=2.2, solid_capstyle="butt")
    ax.set_ylim(-1.2, 1.0)
    ax.set_yticks([-1, 0, 1])
    # raw BAF
    ax = axes[1]
    xs = np.sort(rng.uniform(0, L, 700))
    px, py = [], []
    for x in xs:
        b = card["bands"](x)
        if not b:
            continue
        p = b[rng.integers(len(b))]
        k = rng.binomial(dp, p)
        if 0.15 <= k / dp <= 0.85:
            px.append(x)
            py.append(k / dp)
    ax.scatter(px, py, s=1.2, c=PAL["baf"], lw=0, rasterized=True)
    for yv in (1 / 3, 0.5, 2 / 3):
        ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    # the phased fraction: sites, windows, step, auxiliary tracks
    ax = axes[2]
    for yv in (1 / 3, 0.5, 2 / 3):
        ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
    xm = np.sort(rng.uniform(0, L, 240))
    mu = card["main"](xm)
    ok = np.isfinite(mu)
    k = rng.binomial(dp, np.clip(mu[ok], 0, 1))
    ax.scatter(xm[ok], k / dp, s=1.2, c=PAL["psite"], lw=0, rasterized=True)
    wx, wy = [], []
    for i in range(0, len(xm) - 9, 10):
        sel = ok[i:i + 10]
        if sel.sum() >= 5:
            wx.append(xm[i:i + 10][sel].mean())
            wy.append(rng.binomial(dp * int(sel.sum()), np.clip(mu[i:i + 10][sel].mean(), 0, 1)) / (dp * int(sel.sum())))
    for name, col in (("aux_f", PAL["mat"]), ("aux_m", PAL["pat"])):
        fn = card[name]
        if fn is not None:
            ax.plot(xb, fn(xb), color=col, lw=1.0)
    ax.scatter(wx, wy, s=4, c=PAL["phased"], lw=0)
    ax.plot(xb, mu if len(mu) == len(xb) else card["main"](xb), color=PAL["step"], lw=0.8)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    # the copies
    ax = axes[3]
    for yv in (0, 1, 2):
        ax.axhline(yv, color=PAL["ref"], lw=0.4, ls=":")
    cm, cp = card["copies"]
    ax.plot(xb, cm(xb), color=PAL["mat"], lw=1.1)
    ax.plot(xb, cp(xb), color=PAL["pat"], lw=1.1)
    ax.set_ylim(-0.15, 2.5)
    ax.set_yticks([0, 1, 2])
    # the heterozygosity rate
    ax = axes[4]
    h = card["het"](xb) * 0.6
    ax.scatter(xb, np.clip(h + rng.normal(0, 0.04, len(xb)) * (h > 0), 0, 1), s=3, c=PAL["depth"], lw=0)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    for ax in axes:
        ax.set_xlim(0, L)
    axes[4].set_xticks([0, 50, 100])


def fig_patterns(out_dir):
    """Three figures of four cards each: patterns_copy (gain, loss, copy-neutral LOH, isodisomy), patterns_disomy (heterodisomy,
    a run of homozygosity, a meiotic and a mitotic trisomy) and patterns_sex (47,XXY maternal and paternal, a mosaic 46,XY/47,XXY,
    45,X). Returns the figure records."""
    os.makedirs(out_dir, exist_ok=True)
    cards = _cards()
    out = []
    rows = ("LRR", "BAF", "maternal\nfraction", "copies\nmaternal, paternal", "het rate")
    for name, sel, caption in (
            ("patterns_copy", cards[:4], "Idealised cards (simulated with the noise of a 30x genome): the child's five rows, position in Mb. Gain of the paternal homologue in "
                                        "40% of cells: LRR up by log2(1 + f/2), the raw bands split to 1/(2 + f) and (1 + f)/(2 + f), the maternal fraction on the lower band "
                                        "(paternal homologue in excess), paternal copy number 1.4. Loss of the paternal copy in 60%: LRR down by log2(1 - f/2), the maternal "
                                        "fraction on the upper band at 1/(2 - f), paternal copy number 0.4. Copy-neutral LOH with the maternal copy retained in 40%: LRR flat, "
                                        "bands at (1 +- f)/2, maternal fraction (1 + f)/2, copy numbers 1.4 and 0.6, heterozygosity retained. Maternal isodisomy: LRR flat, no "
                                        "heterozygous calls, maternal fraction 1 (the child homozygous for the maternal allele at every informative site), copy numbers 2 and "
                                        "0, Mendelian errors at every informative site."),
            ("patterns_disomy", cards[4:8], "Maternal heterodisomy: LRR flat and heterozygosity retained (the child heterozygous wherever the mother is), yet the maternal "
                                          "fraction is 1 and the copy numbers 2 and 0; the track read at the father's homozygous sites (pink line) stays at 1/2, since both of "
                                          "the child's alleles there are maternal. Run of homozygosity: no heterozygous calls, no Mendelian errors, no sites for the maternal "
                                          "fraction (the parents share the haplotype), copy number not drawn. Maternal meiotic trisomy with a crossover at 60 Mb: LRR at "
                                          "log2(1.5), bands at 1/3 and 2/3, maternal fraction 2/3, and the pink auxiliary track at 1/3 where the two maternal copies are "
                                          "different homologues (heterozygosity raised there), returning to 2/3 beyond the crossover where they are identical. Mitotic or "
                                          "meiosis II trisomy: the same, with the auxiliary track on the main track throughout."),
            ("patterns_sex", cards[8:12], "The X chromosome, position in Mb; a son's X normally sits at LRR -1 (one maternal copy) with no heterozygous sites and a maternal "
                                         "fraction of 1, a daughter's at 0 with a maternal fraction of 1/2. 47,XXY from a maternal meiosis I error: the X at LRR 0 (two "
                                         "copies), heterozygous wherever the mother is, the maternal fraction 1 (every allele maternal) with the auxiliary track read at the "
                                         "father's sites at 1/2 (two different maternal homologues), copy numbers maternal 2 and paternal 0, and Mendelian errors at the "
                                         "informative sites under the autosomal rules. 47,XXY with the extra X paternal (a paternal meiosis I error: X and Y transmitted "
                                         "together): the X at 0, bands at 1/2, the maternal fraction 1/2, copy numbers 1 and 1. A mosaic 46,XY/47,XXY with the extra X "
                                         "paternal in 40% of cells: the X at log2(1.4/2), bands at 0.29 and 0.71, the maternal fraction 0.71, the paternal copy number 0.4. "
                                         "45,X with the paternal X lost: the X at -1 in a daughter, no heterozygous calls, the maternal fraction 1, copy numbers 1 and 0.")):
        fig, axes = plt.subplots(5, 4, figsize=(180 * MM, 120 * MM), sharex=True, gridspec_kw=dict(height_ratios=[1, 1, 1, 0.8, 0.5], hspace=0.18, wspace=0.28))
        rng = np.random.default_rng(7)
        for j, card in enumerate(sel):
            _draw_card(axes[:, j], card, rng)
            axes[0, j].set_title(card["title"], fontsize=6, loc="left")
            axes[4, j].set_xlabel("Mb")
        for i, r in enumerate(rows):
            axes[i, 0].set_ylabel(r, fontsize=6)
        for ext in ("png", "svg", "pdf"):
            fig.savefig(os.path.join(out_dir, "%s.%s" % (name, ext)), bbox_inches="tight")
        plt.close(fig)
        sub = {"patterns_copy": " (copy number)", "patterns_disomy": " (disomies, trisomies)", "patterns_sex": " (sex chromosomes)"}[name]
        write_sidecar(out_dir, name, "Expected appearance of each event type" + sub, caption, PATTERN_KEYS)
        write_legend(out_dir, name, PATTERN_KEYS)
        out.append(dict(name=name, title="Expected appearance of each event type" + sub, caption=caption, keys=PATTERN_KEYS, png=os.path.join(out_dir, name + ".png")))
    return out
