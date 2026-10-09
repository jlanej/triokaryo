"""The digest: one self-contained page of the cohort's largest events, to hand off. Every trio's run is read (summary.json,
events.tsv, bins.tsv, its figures); the events are ranked by impact - a whole chromosome (an aneuploidy, a uniparental
disomy) first, then an arm, then a stretch by size, constitutional before mosaic - and the first `top` are shown, each with
its facts in plain words, the trio's genome-wide figure (every chromosome of the child, the father and the mother: LRR,
BAF, the phased fraction) and the chromosome's own figure, embedded. Within a tier a change of copy number comes before a
uniparental disomy, that before a run of homozygosity, each constitutional before mosaic. The long tail is counted and left
in events.all.tsv.

What is left out, and why it says so: a stretch under `min_mb` (the floor of what is worth a page), a mosaic segment
under `min_f` of the cells, a run of homozygosity (an LOH stretch) under `min_mb_loh` (consanguinity's runs, not events),
and a parent's X or Y lost in part of the cells (loss with age, common in blood), unless asked for. A whole chromosome is
never left out for its share of cells: a trisomy in a fifth of the cells is an event. Runs without a panel are the first
pass, whose calls are not to be read; they are refused unless allowed."""
import base64
import html
import io
import json
import os
import time

import numpy as np

from .guide import key_table
from .model import NA
from .plots import GENOME_KEYS, KEY, MM, PAL
from .report import CSS, MEMBERS, _finite, _read_summary_tsv, events_of, karyotypes_of, write_tsv

TIER = {"whole": 0, "p": 1, "q": 1, "stretch": 2}
TIER_NAME = {0: "a whole chromosome", 1: "an arm", 2: "a stretch"}
DIGEST_COLS = ("rank", "trio", "sample", "role", "chrom", "start", "end", "span", "bands", "type", "size_mb", "f", "source", "inheritance", "origin",
               "stage", "lrr", "n_bins", "mie_rate", "external", "note", "karyotype", "genome_figure", "chrom_figure")


def _f(x, nd=2):
    return ("%.*f" % (nd, x)) if _finite(x) else "NA"


def size_mb(e):
    return (e.end - e.start) / 1e6


TYPE_ORDER = {"gain": 0, "loss": 0, "UPD": 1, "LOH": 2}


def rank_key(e):
    """Impact order: whole chromosome, arm, stretch; within, constitutional (f >= 0.8) before mosaic; within that, a change of copy
    number before a uniparental disomy before a loss of heterozygosity; then size, then the share of cells."""
    return (TIER.get(e.span, 2), 0 if _finite(e.f) and e.f >= 0.8 else 1, TYPE_ORDER.get(e.type, 3), -size_mb(e), -(e.f if _finite(e.f) else 0.0))


def left_out_reason(e, min_mb, min_f, min_mb_loh, parent_sex_mosaics):
    """'' where the event is listed, else why not."""
    whole = e.span == "whole"
    if e.chrom in ("chrX", "chrY") and e.role != "child" and _finite(e.f) and e.f < 0.8 and not parent_sex_mosaics:
        return "a parent's %s in part of the cells: loss with age in blood (--parent-sex-mosaics lists them)" % ("X" if e.chrom == "chrX" else "Y")
    if e.type in ("LOH",) and not whole and size_mb(e) < min_mb_loh:
        return "a run of homozygosity under %g Mb (--min-mb-loh)" % min_mb_loh
    if not whole and e.span == "stretch" and size_mb(e) < min_mb:
        return "a stretch under %g Mb (--min-mb)" % min_mb
    if not whole and _finite(e.f) and e.f < min_f:
        return "in under %.0f%% of the cells (--min-f)" % (100 * min_f)
    return ""


def what(e):
    """The event in words: 'loss of 2q22.1-q31.1 (32 Mb) in 67% of the cells'."""
    t = {"gain": "gain", "loss": "loss", "LOH": "copy-neutral loss of heterozygosity", "UPD": "uniparental disomy"}.get(e.type, e.type)
    c = e.chrom[3:]
    if e.span == "whole":
        where = ("chromosome " + c) if e.type not in ("UPD", "LOH") else ("chromosome " + c)
        head = {"gain": "an extra copy of %s", "loss": "one copy of %s missing", "UPD": "%s, both copies from one parent", "LOH": "%s without heterozygosity"}.get(e.type, t + " of %s") % where
    elif e.span in ("p", "q"):
        head = "%s of the %s arm of %s (%s%s)" % (t, "short" if e.span == "p" else "long", c, c, e.span)
    else:
        bands = (e.bands or "").strip()
        head = "%s of %s%s (%.1f-%.1f Mb)" % (t, c, (" " + bands) if bands else "", e.start / 1e6, e.end / 1e6)
    cells = "" if not _finite(e.f) or e.f >= 0.8 else ", in %.0f%% of the cells" % (100 * e.f)
    return head + cells


def reading(e, s):
    """One paragraph: inheritance and parent of origin, the meiotic stage, what NGS-DOSE saw, the member's sex chromosomes."""
    bits = []
    if e.role == "child":
        inh = e.inheritance or ""
        if inh.startswith("new"):
            bits.append("New in the child: neither parent carries it.")
        elif inh.startswith("inherited"):
            bits.append(inh[0].upper() + inh[1:] + ".")
        elif inh:
            bits.append(inh[0].upper() + inh[1:] + ".")
    else:
        inh = e.inheritance or ""
        if inh:
            bits.append(("The %s's event, " % e.role) + inh + ".")
    if e.origin:
        bits.append("Parent of origin from the phased alleles: %s%s." % (e.origin, (" (%s)" % e.stage) if e.stage else ""))
    if e.external:
        bits.append("Other methods at this place: %s." % e.external)
    elif e.type in ("gain", "loss"):
        bits.append("No supplied event (NGS-DOSE) overlaps it.")
    if e.note:
        bits.append("Flag: %s." % e.note)
    sx = (s.get("sex_chromosomes") or {}).get(e.role) or {}
    chk = sx.get("sex_check") or ""
    if chk and chk != "agrees":
        bits.append("Sex chromosomes of this member: %s." % chk)
    return " ".join(bits)


def embed_png(path, max_width=1500):
    """A PNG on disk as a data URI, downscaled to at most max_width pixels wide (the figures are drawn at 300 dpi)."""
    try:
        from PIL import Image
        im = Image.open(path)
        if im.width > max_width:
            im = im.convert("RGB")
            im.thumbnail((max_width, int(im.height * max_width / im.width)))
            buf = io.BytesIO()
            im.save(buf, format="PNG", optimize=True)
            data = buf.getvalue()
        else:
            data = open(path, "rb").read()
    except Exception:                                          # no Pillow, or a file Pillow cannot read: as it is
        data = open(path, "rb").read()
    return "data:image/png;base64," + base64.b64encode(data).decode()


def _bins(run_dir):
    import csv
    p = os.path.join(run_dir, "bins.tsv")
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(open(p), delimiter="\t"))


def fallback_figure(rows, events, trio_name, genome, highlight=None, chrom=None):
    """Without the run's own figures (a run made with --no-figures): every chromosome - or one - of the three members from
    bins.tsv, the LRR per bin and the band deviation (median |BAF - 1/2| at the heterozygous sites, the BAF's departure from
    one half), the called events as lines and the event in question shaded. Returns PNG bytes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [r for r in rows if (chrom is None or r["chrom"] == chrom)]
    chroms = [c for c in genome.chroms if any(r["chrom"] == c for r in rows)]
    if not rows or not chroms:
        return None
    off, total = {}, 0
    for c in chroms:
        off[c] = total
        total += genome.length[c]
    xs = np.array([off[r["chrom"]] + (float(r["start"]) + float(r["end"])) / 2.0 for r in rows])
    fig, axes = plt.subplots(6, 1, figsize=(180 * MM, (118 if chrom is None else 95) * MM), sharex=True,
                             gridspec_kw=dict(height_ratios=[1, 0.7, 1, 0.7, 1, 0.7], hspace=0.14))
    for m, role in enumerate(MEMBERS):
        lrr = np.array([float(r.get("%s_lrr_gc" % role, "nan") or "nan") if r.get("%s_lrr_gc" % role, "NA") != "NA" else np.nan for r in rows])
        bd = np.array([float(r["%s_bdev" % role]) if r.get("%s_bdev" % role, "NA") != "NA" else np.nan for r in rows])
        for ax, y, lim, lab in ((axes[2 * m], lrr, (-1.5, 1.2), "%s\nLRR" % role), (axes[2 * m + 1], bd, (0, 0.5), "|BAF-1/2|")):
            for k, c in enumerate(chroms):
                if k % 2:
                    ax.axvspan(off[c], off[c] + genome.length[c], color=PAL["band"], lw=0)
            if highlight is not None and highlight.chrom in off:
                ax.axvspan(off[highlight.chrom] + highlight.start, off[highlight.chrom] + highlight.end, color=PAL["child"], alpha=0.10, lw=0)
            ax.axhline(0 if lim[0] < 0 else 0.0, color=PAL["ref"], lw=0.5, ls=":")
            ax.scatter(xs, y, s=1.4, c=PAL["depth"] if lim[0] < 0 else PAL["baf"], lw=0, rasterized=True)
            if lim[0] < 0:
                for e in events:
                    if e.role == role and e.chrom in off and e.type in ("gain", "loss", "LOH", "UPD"):
                        ax.plot([off[e.chrom] + e.start, off[e.chrom] + e.end], [e.lrr if (e.type in ("gain", "loss") and _finite(e.lrr)) else 0.0] * 2,
                                color=PAL["gain"] if e.type == "gain" else PAL["loss"] if e.type == "loss" else PAL["loh"], lw=2.0, solid_capstyle="butt")
            ax.set_ylim(*lim)
            ax.set_ylabel(lab, fontsize=6)
    axes[-1].set_xlim(0, total)
    axes[-1].set_xticks([off[c] + genome.length[c] / 2 for c in chroms])
    axes[-1].set_xticklabels([c[3:] for c in chroms], fontsize=5 if chrom is None else 7)
    axes[-1].set_xlabel("chromosome" if chrom is None else "%s (bins)" % chrom[3:])
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def write_digest(out, run_dirs, events_path=None, genome_name="grch38", top=25, min_mb=10.0, min_f=0.2, min_mb_loh=20.0, parent_sex_mosaics=False,
                 allow_first_pass=False, max_width=1500, title="", log=None):
    """digest.html, digest.tsv (the listed events) and digest_left_out.tsv (every other event, with why) under `out`."""
    from .external import match_external, read_events
    from .genome import genome as load_genome
    os.makedirs(out, exist_ok=True)
    G = load_genome(genome_name)
    runs, first_pass = [], []
    for d in run_dirs:
        sj = os.path.join(d, "summary.json")
        if not os.path.exists(sj):
            continue
        s = json.load(open(sj))
        s["run"] = d
        s["tsv"] = _read_summary_tsv(os.path.join(d, "summary.tsv"))
        s["events_obj"] = events_of(s)
        if str(s["tsv"].get("panel", "")).lower() not in ("1", "true"):
            first_pass.append(s["trio"])
            if not allow_first_pass:
                continue
        runs.append(s)
    if first_pass and not allow_first_pass and log:
        log("%d run(s) without a panel left out: the first pass, whose calls are not to be read (--allow-first-pass to keep them)" % len(first_pass))
    if not runs:
        raise SystemExit("no run with a panel under: " + " ".join(run_dirs) + ("" if not first_pass else " (%d first-pass runs; --allow-first-pass)" % len(first_pass)))
    events = [(s, e) for s in runs for e in s["events_obj"]]
    ext_n = 0
    if events_path:
        ext = read_events(events_path, G)
        samples = {m for s in runs for m in s["members"]}
        ext = [x for x in ext if x.sample in samples]
        match_external([e for _, e in events], ext)
        ext_n = len(ext)
    listed, left = [], []
    for s, e in events:
        why = left_out_reason(e, min_mb, min_f, min_mb_loh, parent_sex_mosaics)
        (left if why else listed).append((s, e, why))
    listed.sort(key=lambda t: rank_key(t[1]))
    shown = listed[:top]
    beyond = listed[top:]
    # the rows
    def row(i, s, e, why=""):
        kar = karyotypes_of(s["events_obj"], s["tsv"], G)
        return dict(rank=i, trio=s["trio"], sample=e.sample, role=e.role, chrom=e.chrom, start=e.start, end=e.end, span=e.span, bands=e.bands, type=e.type,
                    size_mb=round(size_mb(e), 1), f=e.f, source=e.source, inheritance=e.inheritance, origin=e.origin, stage=e.stage, lrr=e.lrr, n_bins=e.n_bins,
                    mie_rate=e.mie_rate, external=e.external, note=e.note, karyotype=kar.get(e.role, ""), left_out=why)
    write_tsv(os.path.join(out, "digest.tsv"), list(DIGEST_COLS[:-2]), [row(i + 1, s, e) for i, (s, e, _) in enumerate(shown)])
    write_tsv(os.path.join(out, "digest_left_out.tsv"), ["left_out"] + [c for c in DIGEST_COLS[:-2] if c != "rank"],
              [row("", s, e, why) for s, e, why in left] + [row("", s, e, "beyond the first %d (--top)" % top) for s, e, _ in beyond])
    # the page
    w = []
    esc = html.escape
    w.append('<!doctype html><html><head><meta charset="utf-8"><title>%s</title><style>%s%s</style></head><body>' % (
        esc(title or "triokaryo digest"), CSS,
        ".ev{border:1px solid #ddd;border-radius:6px;padding:0.6em 1em;margin:1.2em 0}.ev h2{margin:0.2em 0 0.4em;font-size:1.15em}"
        ".facts td:first-child{color:#666;white-space:nowrap}.box{background:#f6f6f6;border:1px solid #e3e3e3;padding:0.6em 1em;margin:1em 0}"
        ".tag{display:inline-block;background:#eee;border-radius:9px;padding:0 8px;font-size:11px;margin-left:6px}.muted{color:#777}"))
    n_trios = len(runs)
    mock = [s.get("mock_note") for s in runs if s.get("mock_note")]
    if mock:
        w.append('<div class="mock">MOCK DATA: %s</div>' % esc(mock[0]))
    w.append("<h1>%s</h1>" % esc(title or "The largest chromosomal events in %d trios" % n_trios))
    w.append('<p class="muted">triokaryo digest, made %s from %d trio runs%s. Every event is a reading of the trio VCF\'s depth and allele fractions, '
             'to be confirmed by a clinical karyotype or a chromosomal microarray before it is reported.</p>' % (
                 time.strftime("%Y-%m-%d"), n_trios, (", %d supplied events (NGS-DOSE) matched" % ext_n) if ext_n else ""))
    w.append('<div class="box"><p><b>What is on this page.</b> The %d largest events of %d called in these trios, ranked by impact: a whole chromosome '
             '(an aneuploidy, a uniparental disomy) first, then an arm, then a stretch by size; within each, an event in every cell before one '
             'in part of the cells, and a change of copy number before a uniparental disomy before a run of homozygosity. Left out, and counted in digest_left_out.tsv: stretches under %g Mb, mosaic segments in under %.0f%% of the cells, '
             'runs of homozygosity under %g Mb%s. A whole chromosome is never left out for its share of cells.</p>'
             '<p><b>How to read a figure.</b> Three members, child then father then mother. <b>LRR</b> is the depth per bin against the member\'s own '
             'autosomal median, on a log2 scale: 0 is two copies, about +0.58 three, about -1 one; a step up is a gain, a step down a loss, and a '
             'mosaic sits in between in proportion to the share of cells. <b>BAF</b> is the alternate-allele fraction at the member\'s heterozygous '
             'sites: a band at one half means two balanced copies; a split into two bands means an imbalance (a gain, a loss, or a copy-neutral '
             'loss of heterozygosity when the LRR stays at 0), and the <b>phased fraction</b> (green) says which parent\'s homologue is in excess. '
             'The event under discussion is named by its place; in a figure drawn from the bins alone (a run made without figures) it is shaded.</p>%s</div>'
             % (len(shown), len(events), min_mb, 100 * min_f, min_mb_loh,
                "" if parent_sex_mosaics else ", and a parent's X or Y lost in part of the cells (loss with age)", key_table(GENOME_KEYS)))
    # the summary table
    w.append("<h2>The events</h2>")
    w.append("<table><tr>" + "".join("<th>%s</th>" % h for h in ("#", "trio", "member", "event", "size (Mb)", "cells", "how it is seen", "inheritance, parent of origin", "other methods")) + "</tr>")
    for i, (s, e, _) in enumerate(shown):
        w.append("<tr><td><a href=\"#ev%d\">%d</a></td><td>%s</td><td>%s (%s)</td><td>%s</td><td>%.1f</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>" % (
            i + 1, i + 1, esc(s["trio"]), esc(e.sample), esc(e.role), esc(what(e)), size_mb(e), ("%.0f%%" % (100 * e.f)) if _finite(e.f) else "NA",
            esc(e.source or ""), esc("; ".join(x for x in (e.inheritance or "", e.origin or "") if x)), esc((e.external or "")[:120])))
    w.append("</table>")
    if beyond:
        w.append('<p class="muted">%d more events pass the floors but lie beyond the first %d; %d are left out by the floors. All are in digest_left_out.tsv and the cohort\'s events.all.tsv.</p>' % (len(beyond), top, len(left)))
    else:
        w.append('<p class="muted">%d events are left out by the floors (digest_left_out.tsv).</p>' % len(left))
    # each event
    genome_shown = {}
    for i, (s, e, _) in enumerate(shown):
        kar = karyotypes_of(s["events_obj"], s["tsv"], G)
        w.append('<div class="ev" id="ev%d"><h2>%d. Trio %s, %s %s: %s<span class="tag">%s</span></h2>' % (
            i + 1, i + 1, esc(s["trio"]), esc(e.role), esc(e.sample), esc(what(e)), esc(TIER_NAME[TIER.get(e.span, 2)] + (", in every cell" if _finite(e.f) and e.f >= 0.8 else ", mosaic"))))
        w.append("<p>%s</p>" % esc(reading(e, s)))
        facts = [("cells carrying it", "%s (by depth %s, by BAF %s, by phase %s)" % (_f(e.f), _f(getattr(e, "f_lrr", NA)), _f(getattr(e, "f_baf", NA)), _f(getattr(e, "f_phase", NA)))),
                 ("LRR over the event", "%s (SE %s) over %s bins" % (_f(e.lrr, 3), _f(getattr(e, "lrr_se", NA), 3), e.n_bins)),
                 ("phased fraction shift", "%s (SE %s), %s phased sites; homologues: %s" % (_f(getattr(e, "phase_shift", NA), 3), _f(getattr(e, "phase_se", NA), 3), getattr(e, "n_phased", ""), getattr(e, "homologues", "") or "-")),
                 ("parent of origin", "%s (log-likelihood ratio %s over %s sites)%s" % (e.origin or "-", _f(getattr(e, "origin_llr", NA), 1), getattr(e, "origin_n", ""), (", " + e.stage) if e.stage else "")),
                 ("heterozygosity, Mendelian errors", "het rate %s (%s of the panel's), MIE rate %s over %s called sites" % (_f(getattr(e, "het_rate", NA), 3), _f(getattr(e, "het_rate_rel", NA)), _f(e.mie_rate, 3), getattr(e, "n_called", ""))),
                 ("crossovers along it", "%s %s" % (getattr(e, "n_crossovers", ""), getattr(e, "crossovers", "") or "")),
                 ("the member's karyotype", kar.get(e.role, "")),
                 ("the trio", "child %s, father %s, mother %s; sexes %s" % (s["members"][0], s["members"][1], s["members"][2], ",".join(s.get("sexes") or [])))]
        w.append('<table class="facts">' + "".join("<tr><td>%s</td><td>%s</td></tr>" % (esc(k), esc(str(v))) for k, v in facts) + "</table>")
        fdir = os.path.join(s["run"], "figures")
        gpng, cpng = os.path.join(fdir, "genome.png"), os.path.join(fdir, "chrom_%s.png" % e.chrom)
        if s["trio"] in genome_shown:
            w.append('<p class="muted">Genome-wide figure of this trio: see event <a href="#ev%d">%d</a>.</p>' % (genome_shown[s["trio"]], genome_shown[s["trio"]]))
        else:
            genome_shown[s["trio"]] = i + 1
            if os.path.exists(gpng):
                w.append('<figure><img src="%s" alt="genome-wide, trio %s"><figcaption>Trio %s, every chromosome: the child (a, b), the father (c, d) and the mother (e, f), LRR and BAF with the '
                         'phased fraction; (g) the child\'s depth over the parents\' mean; (h) the child\'s maternal and paternal copies. The event: %s.</figcaption></figure>'
                         % (embed_png(gpng, max_width), esc(s["trio"]), esc(s["trio"]), esc(what(e))))
            else:
                png = fallback_figure(_bins(s["run"]), s["events_obj"], s["trio"], G, highlight=e)
                if png:
                    w.append('<figure><img src="data:image/png;base64,%s" alt="genome-wide, trio %s"><figcaption>Trio %s, every chromosome, drawn from the bins (the run was made without figures): '
                             'per member the LRR and the band deviation |BAF - 1/2|; the event shaded.</figcaption></figure>' % (base64.b64encode(png).decode(), esc(s["trio"]), esc(s["trio"])))
        if os.path.exists(cpng):
            w.append('<figure><img src="%s" alt="%s, trio %s"><figcaption>%s in the three members: LRR, the BAF with the child\'s informative sites coloured by the parent of the alternate '
                     'allele, the phased fraction with its step fit, the copies per homologue, the heterozygosity rate.</figcaption></figure>' % (embed_png(cpng, max_width), esc(e.chrom), esc(s["trio"]), esc(e.chrom)))
        else:
            png = fallback_figure(_bins(s["run"]), s["events_obj"], s["trio"], G, highlight=e, chrom=e.chrom)
            if png:
                w.append('<figure><img src="data:image/png;base64,%s" alt="%s, trio %s"><figcaption>%s from the bins: per member the LRR and the band deviation; the event shaded.</figcaption></figure>'
                         % (base64.b64encode(png).decode(), esc(e.chrom), esc(s["trio"]), esc(e.chrom)))
        w.append("</div>")
    w.append('<h2>Left out</h2><p class="muted">%s</p>' % esc("; ".join("%d %s" % (n, why) for why, n in sorted(
        {why: sum(1 for _, _, x in left if x == why) for _, _, why in left}.items(), key=lambda kv: -kv[1])) or "nothing"))
    w.append('<p class="muted">Tables beside this page: digest.tsv (these events), digest_left_out.tsv (every other, with why). The cohort\'s full report, every event and '
             'every trio, is the cohort page (index.html) and events.all.tsv.</p></body></html>')
    with open(os.path.join(out, "digest.html"), "w") as fh:
        fh.write("\n".join(w))
    if log:
        log("digest: %d of %d events listed (%d left out by the floors, %d beyond --top %d) -> %s/digest.html" % (len(shown), len(events), len(left), len(beyond), top, out))
    return shown, left
