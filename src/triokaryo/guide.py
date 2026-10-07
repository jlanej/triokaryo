"""The guide: how to read every figure, colour, call and column triokaryo writes - one self-contained page, with the
pattern cards (what each kind of event looks like) drawn from idealised tracks. Written beside every trio's page and
the cohort report (guide.html), and by `triokaryo guide --out`."""
import base64
import html
import os
import shutil
import tempfile

from . import __version__
from .plots import CHROM_KEYS, DIRECTION, GENOME_KEYS, KEY, MARK

SITE = "https://jlanej.github.io/triokaryo/"

COLUMN_DOCS = [
    ("sample, role", "the member: child, father or mother"),
    ("chrom, start, end", "the event's segment in bp: the start is the first bin's left edge, the end the last bin's right edge (a phased find: its windows' first and last sites)"),
    ("start_fine, end_fine, edge_sites", "the edges at site resolution, from the phased sites either side of each bin edge (NA where the event runs to the chromosome's end, or the sites are too few); the sites the narrower edge used"),
    ("span", "whole (90% of the chromosome's usable bins), p or q (90% of an arm and under half of the other), else stretch"),
    ("type", "gain, loss, LOH (a copy-neutral loss of heterozygosity), UPD (a uniparental heterodisomy: both copies from one parent, two different homologues, the heterozygosity kept)"),
    ("source", "depth (the LRR called it), bands (the folded bands or the heterozygosity rate called it), phased (the phased scan: what the depth did not call)"),
    ("f", "the share of cells: the depth's for a gain or loss called by the depth, else the bands' (2 d for an LOH; 1 where the heterozygosity rate says in every cell)"),
    ("f_lrr, f_baf, f_phase", "the share from the depth (gain 2(2^LRR - 1), loss 2(1 - 2^LRR)), from the folded bands' deviation d (gain 4d/(1 - 2d), loss 4d/(1 + 2d), LOH 2d) and from the phased shift (the same formulas on |shift|)"),
    ("lrr, lrr_se, n_bins", "the segment's mean LRR, its standard error from the chromosome's noise, its bins"),
    ("d_hat, llr_baf", "the folded bands' deviation from one half at the segment's heterozygous sites (a depth-aware maximum likelihood) and its log-likelihood ratio against zero"),
    ("phase_shift, phase_se, n_phased", "the phased fraction's shift from one half over the event - the maternal allele's in the child, the transmitted allele's in a parent - its binomial error, and the sites"),
    ("origin, origin_llr, origin_n", "the child's parent of origin from the opposite-homozygote sites by a likelihood ratio (maternal over paternal for a gain or LOH; paternal-lost over maternal-lost for a loss), and the sites"),
    ("origin_phase", "what the phased shift's sign says: for the child the extra, lost or retained copy's parent (or, for a UPD, whose both copies are); for a parent whether the duplicated, lost or retained homologue is the one passed to the child. Needs 40 sites and 3 standard errors"),
    ("homologues, hetero_share", "a child's gain, LOH or UPD with a parent named and a shift of 0.05 or more: whether that parent's two copies are one homologue (the auxiliary track follows the main one) or two (it parts: 1/3 against 2/3 along a trisomy), and the share of the event's windows where they differ"),
    ("het_rate, het_rate_rel, n_het, n_called", "the heterozygosity rate in the segment, its ratio to the member's own (or to the panel's), and the counts"),
    ("mie_rate", "Mendelian errors among the segment's confident sites: a deletion or an isodisomy breaks Mendel at every opposite-homozygote site, a heterodisomy too, a trisomy or a mosaic does not, a run of homozygosity has none"),
    ("inheritance", "the child's: inherited from the father / the mother (the same event in that parent; 'in a share of the parent's cells' when mosaic there), or new; a parent's: passed to the child or not"),
    ("external", "the labels of the given events (--events) the segment overlaps by half or more"),
    ("note", "the flags: 'a run of homozygosity ... not a uniparental disomy' (an LOH in every cell with no Mendelian errors); 'Mendelian errors at the informative sites: a uniparental isodisomy'; 'no heterozygous calls: in every cell'; 'joined across n bin(s) the panel left out'; 'the phased bands read a share of x% against the depth's y%: the depth's call may be an artefact'; 'from the phased bands: the depth leans the same way but under its own threshold' / 'the depth is flat'; 'a uniparental heterodisomy ...'; 'the phased bands and the opposite-homozygote sites name different parents'; 'a loss on a single X'"),
]

ROWS = [
    ("1 · LRR", "log2 of each bin's median depth over the member's autosomal median, GC- and panel-corrected, with its step fit (black) and the calls as lines: a gain at the segment's mean (vermillion), a loss (blue), a copy-neutral LOH or heterodisomy drawn at zero (reddish purple). Brackets above: events given from elsewhere.",
     "a plateau away from zero over several bins; log2(1 + f/2) for a gain, log2(1 - f/2) for a loss; a flat LRR under a parted band is copy-neutral"),
    ("2 · BAF", "the raw alt-allele fraction at the member's heterozygous sites (a sample). In the child, the sites where the parents are opposite homozygotes are coloured by the parent of the alt allele (pink the mother's, orange the father's).",
     "two bands parting from one half; under a gain the duplicated parent's alleles sit on the upper band; under a loss the retained parent's at the top; an empty stretch is a loss of heterozygosity in every cell"),
    ("3 · phased fraction", "the maternal allele's fraction along the child, the transmitted allele's along a parent: every phased site (faint grey), the pooled windows (green; hollow where parted in two or more members) and the step fit (black); the thin pink and orange lines are the auxiliary tracks.",
     "one band, not two: %s A thin line leaving the main track means the child carries two different homologues of that parent there (a meiotic error, a heterodisomy), returning at each crossover." % DIRECTION),
    ("4 · copies", "the LRR step fit's copies (2 x 2^LRR) split by the fraction's step fit: maternal (pink) and paternal (orange) along the child, passed to the child (vermillion) and not passed (sky blue) along a parent.",
     "a trisomy reads 2 and 1 with its parent named; a deletion 1 and 0; an isodisomy or heterodisomy 2 and 0; a mosaic event the fractional copies; a run of homozygosity is blank (no phased sites)"),
    ("5 · het rate", "heterozygous calls per confident call, per bin.", "zero under a loss of heterozygosity in every cell (an isodisomy, a run of homozygosity, a deletion); raised where the child carries two different homologues of one parent"),
]

KINDS = [
    ("gain (one homologue duplicated), share f", "up by log2(1 + f/2)", "1/(2 + f) and (1 + f)/(2 + f)", "on the duplicated parent's band: (1 + f)/(2 + f) maternal, 1/(2 + f) paternal", "the parent's copies 1 + f, the other's 1", "unchanged", "none", "the extra copy's parent"),
    ("loss of one copy, share f", "down by log2(1 - f/2)", "(1 - f)/(2 - f) and 1/(2 - f)", "on the retained parent's band: 1/(2 - f) maternal if the paternal copy is lost", "the lost parent's copies 1 - f", "unchanged (f < 1); zero at f = 1", "at f = 1, at every opposite-homozygote site", "the lost copy's parent"),
    ("copy-neutral LOH, share f", "flat", "(1 - f)/2 and (1 + f)/2", "(1 + f)/2 on the retained parent's side", "retained 1 + f, replaced 1 - f", "unchanged (f < 1)", "none (f < 1)", "the retained copy's parent"),
    ("uniparental isodisomy (one homologue twice)", "flat", "none: no heterozygous calls", "1 (the child homozygous for that parent's allele)", "2 and 0", "zero", "at every opposite-homozygote site", "whose two copies"),
    ("uniparental heterodisomy (two different homologues)", "flat", "one band at one half: the heterozygosity kept", "1 on the main track; the auxiliary track at one half", "2 and 0", "unchanged or raised", "at every opposite-homozygote site", "whose two copies; found by the phased scan alone"),
    ("run of homozygosity (identical by descent)", "flat", "none", "no sites (the parents share the haplotype)", "blank", "zero", "none", "no parent: not a disomy"),
    ("meiotic trisomy (two different homologues)", "up by log2(1.5)", "1/3 and 2/3", "2/3 on the main track; the auxiliary track at 1/3 where the two copies differ, 2/3 where a crossover made them one", "2 and 1", "raised where the copies differ", "none", "the parent of the two copies; 'two different homologues'"),
    ("mitotic or meiosis II trisomy (one homologue twice)", "up by log2(1.5)", "1/3 and 2/3", "2/3 on both tracks", "2 and 1", "unchanged", "none", "'one homologue throughout'"),
]

RULES = [
    ("The depth's calls", "a segment of at least five bins whose mean LRR is at least 0.07 in size (a gain or loss in 10% of cells) and three standard errors from zero; bins the panel cannot pin are left out; pieces parted by a masked gap or a small step are joined."),
    ("The bands' calls", "a copy-neutral LOH where the bands part (d at least 0.04 with a likelihood ratio of 10: a share of about 8%) or the heterozygosity rate falls to 0.35 of the member's own (in every cell), over bins the depth did not call; the heterozygosity rate takes precedence over a d from a few sites."),
    ("The phased scan", "on the main phased track's windows not already inside an event: a segment of eight windows (or a whole short chromosome) shifted by at least 0.015 (a gain or loss in 6% of cells; 0.025, 5%, where the depth is flat and the bands are the only evidence), five empirical standard errors, three quarters of its windows and its median the same way, spanning at least 2 Mb; typed by the depth's lean (three standard errors) or read as a copy-neutral LOH; a shift of 0.4 or more with the heterozygosity kept is a heterodisomy. Windows parted in two or more members are left out first. Segments set aside are in phased_rejected.tsv with the reason."),
    ("The direction", DIRECTION),
    ("The share of cells three ways", "from the depth, from the folded bands and from the phased shift (the formulas under the columns below). They should agree; a phased share under half the depth's, over 500 sites, flags the depth's call as a likely artefact (GC, the panel's edge)."),
    ("The parent of origin two ways", "the likelihood ratio over the opposite-homozygote sites (origin, origin_llr) and the phased sign (origin_phase, over three to four times the sites). A disagreement is noted."),
    ("One homologue or two", "for a child's gain, LOH or UPD with its parent named: the auxiliary track read at the other parent's homozygous sites follows the main one where the two copies are one homologue and parts from it where they are two; the share of windows where it parts is hetero_share."),
    ("A run of homozygosity against a disomy", "an LOH in every cell with no Mendelian errors is a run of homozygosity (both copies identical by descent: the parents share the haplotype); with errors at the opposite-homozygote sites it is a uniparental isodisomy. A heterodisomy keeps the heterozygosity and has the errors."),
    ("The X", "read against the member's own X median, so a 47,XXY's X is 'two copies' and the check against the pedigree's sex says so; a single X has no heterozygous sites to phase; the pseudoautosomal bins are left out."),
    ("Parents' events", "a parent's gain, loss or LOH is read the same way; 'passed to the child' when the child carries the same event; the phased sign says whether it lies on the homologue the child received."),
]


def _b64(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def _swatch(col, mark):
    if mark == "point":
        return '<span style="display:inline-block;width:11px;height:11px;border-radius:50%%;background:%s;vertical-align:-1px"></span>' % col
    if mark == "line":
        return '<span style="display:inline-block;width:22px;height:3px;background:%s;vertical-align:3px"></span>' % col
    return '<span style="display:inline-block;width:22px;height:8px;border:2px solid %s;border-bottom:none;vertical-align:-1px"></span>' % col


def key_table(keys):
    """The colour key as an HTML table: swatch, mark, meaning."""
    h = ["<table class=\"key\"><tr><th></th><th>mark</th><th>meaning</th></tr>"]
    for k in keys:
        col, mk, text = KEY[k]
        h.append("<tr><td>%s</td><td>%s <code>%s</code></td><td>%s</td></tr>" % (_swatch(col, mk), MARK[mk], html.escape(col), html.escape(text)))
    h.append("</table>")
    return "\n".join(h)


CSS = ("html{background:#fff;color-scheme:light}body{font-family:Helvetica,Arial,sans-serif;max-width:1100px;margin:1.5em auto;padding:0 1em;color:#222;background:#fff;line-height:1.45}table{border-collapse:collapse;font-size:13px;"
       "margin:0.6em 0}th,td{border:1px solid #ddd;padding:3px 7px;text-align:left;vertical-align:top}th{background:#f3f3f3}img{max-width:100%}figure{margin:1.2em 0}"
       "figcaption{font-size:12.5px;color:#444}h2{margin-top:1.6em;border-bottom:1px solid #ddd;padding-bottom:3px}code{background:#f3f3f3;padding:1px 4px;border-radius:3px}"
       "table.key td:first-child{width:28px;text-align:center}.toc a{margin-right:14px}p.lead{font-size:15px}dl dt{font-weight:bold;margin-top:0.6em}dl dd{margin-left:0}")


def write_guide(out_html, figures_dir=None):
    """The guide page. figures_dir: where to keep the pattern figures (PNG, SVG, PDF, sidecar, legend); a temporary
    directory, deleted afterwards, when not given (the page embeds them either way)."""
    from .patterns import fig_patterns
    tmp = None
    if figures_dir is None:
        tmp = figures_dir = tempfile.mkdtemp(prefix="triokaryo_guide_")
    try:
        figs = fig_patterns(figures_dir)
        w = ['<!doctype html><html><head><meta charset="utf-8"><title>triokaryo: how to read it</title><style>%s</style></head><body>' % CSS,
             "<h1>How to read triokaryo's output</h1>",
             "<p class=\"lead\">triokaryo reads a trio's VCF for large chromosomal events with three signals: the <b>depth</b> at each site (the LRR per bin, the copy "
             "state), the <b>B-allele bands</b> at heterozygous sites (the share of cells a second time, and the copy-neutral events the depth cannot see) and "
             "<b>transmission phasing</b> (which parent each allele came from: a signed track that names the parent of origin, tells one homologue from two, "
             "finds what the depth misses, and gives the edges at site resolution). This page says what every row, colour, call and column means. Worked examples "
             "on two public trios: <a href=\"%s\">%s</a>.</p>" % (SITE, SITE),
             "<p class=\"toc\"><a href=\"#rows\">The rows</a><a href=\"#colours\">The colours</a><a href=\"#kinds\">What each event looks like</a><a href=\"#rules\">The rules</a>"
             "<a href=\"#columns\">The columns</a><a href=\"#cohort\">The cohort report</a></p>",
             "<h2 id=\"rows\">The rows of a figure</h2>",
             "<p>A chromosome figure has three columns (child, father, mother) and five rows. The genome figure has, per member, the first two rows (the phased "
             "fraction's pooled windows drawn over the raw B-allele fraction), then the child's depth over the parents' mean (the within-family difference: zero "
             "where the child inherited what the parents carry) and the child's maternal and paternal copies along the genome.</p>",
             "<table><tr><th>row</th><th>what is drawn</th><th>what to look for</th></tr>"]
        for r, what, look in ROWS:
            w.append("<tr><td><b>%s</b></td><td>%s</td><td>%s</td></tr>" % (html.escape(r), html.escape(what), html.escape(look)))
        w.append("</table>")
        w.append("<h2 id=\"colours\">The colours</h2><p>Okabe-Ito throughout; the same colour means the same thing in every figure. The image carries no legend: "
                 "each figure's sidecar (<code>&lt;name&gt;.txt</code>) and legend image (<code>legends/</code>) hold its key.</p>")
        w.append(key_table(list(dict.fromkeys(GENOME_KEYS + CHROM_KEYS))))
        w.append("<h2 id=\"kinds\">What each kind of event looks like</h2>")
        for n, fg in enumerate(figs, 1):
            w.append('<figure><img src="%s" alt="%s"><figcaption><b>Pattern cards %d.</b> %s</figcaption></figure>' % (_b64(fg["png"]), html.escape(fg["name"]), n, html.escape(fg["caption"])))
        w.append("<table><tr><th>event</th><th>LRR</th><th>raw bands</th><th>phased fraction</th><th>copies</th><th>het rate</th><th>Mendelian errors</th><th>what is read</th></tr>")
        for row in KINDS:
            w.append("<tr>" + "".join("<td>%s</td>" % html.escape(c) for c in row) + "</tr>")
        w.append("</table>")
        w.append("<h2 id=\"rules\">The rules the calls follow</h2><dl>")
        for t, text in RULES:
            w.append("<dt>%s</dt><dd>%s</dd>" % (html.escape(t), html.escape(text)))
        w.append("</dl>")
        w.append("<h2 id=\"columns\">The columns of events.tsv</h2><table><tr><th>column</th><th>meaning</th></tr>")
        for c, text in COLUMN_DOCS:
            w.append("<tr><td><code>%s</code></td><td>%s</td></tr>" % (html.escape(c), html.escape(text)))
        w.append("</table>")
        w.append("<p><code>summary.tsv</code>: the trio, the records and sites used, the genome's Mendelian-error rate, per member the events by type, the X copies and "
                 "the check against the pedigree's sex, the median depth, the phased sites, the windows parted in everyone, the phased finds; the parameters. "
                 "<code>bins.tsv</code>: every bin's depth, LRR, calls, heterozygosity and band deviation per member, with the panel's values and masks. "
                 "<code>phased.tsv</code>: every window of the phased tracks (the fraction, its error, the step fit, the copies, the auxiliary tracks). "
                 "<code>phased_rejected.tsv</code>: the phased scan's segments set aside, with the reason. <code>external.tsv</code>: the given events and their matches.</p>")
        w.append("<h2 id=\"cohort\">The cohort report</h2><p>Written by <code>triokaryo cohort</code> over every trio's run: the counts (trios, events by type, by "
                 "member, by source, new and inherited, the phased finds, the depth calls the bands doubt, the X readings against the pedigree); the landscape "
                 "figure (events per chromosome, and one row per trio with each event a bar coloured by type - the child's thick, a parent's thin above or below, "
                 "a black line over a bar for an event matched from elsewhere); the table of every event (click a heading to sort, type in the box to filter), "
                 "each with links to its trio's page and chromosome figure; the trios with their quality readings; the concordance with the given events; the "
                 "segments the phased scan set aside, by region; and this guide.</p>")
        w.append("<p style=\"color:#666;font-size:12px\">triokaryo %s</p></body></html>" % html.escape(__version__))
        with open(out_html, "w") as fh:
            fh.write("\n".join(w))
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return out_html
