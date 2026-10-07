"""The guide: the meaning of every figure row, colour, call and column, as one self-contained page with pattern cards
(the expected appearance of each event type) drawn from idealised tracks. Written beside every trio page and the
cohort report (guide.html), and by `triokaryo guide --out`."""
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
    ("chrom, start, end", "the event in bp: the first bin's start (0-based) and the last bin's end; for a phased-scan event, its first and last window sites"),
    ("start_fine, end_fine, edge_sites", "boundaries at site resolution from the phased sites on either side of each bin edge (NA at a chromosome end or with too few sites), and the smaller of the two site counts used"),
    ("span", "whole (90% of the chromosome's usable bins), p or q (90% of the arm's bins and under half of the other arm's), else stretch"),
    ("bands", "the cytogenetic bands spanned (GRCh38, UCSC cytoBand) in ISCN form, from the site-resolution boundaries where available: q22.1q31.1, or one band"),
    ("type", "gain, loss, LOH (copy-neutral loss of heterozygosity), UPD (uniparental heterodisomy: both copies from one parent as two different homologues, heterozygosity retained)"),
    ("source", "depth (LRR segmentation), bands (band deviation or heterozygosity rate), phased (the phased scan: events the depth did not call)"),
    ("f", "cell fraction: f_lrr for a depth-called gain or loss, otherwise f_baf (2d for LOH; 1 where the heterozygosity rate indicates a constitutional event)"),
    ("f_lrr, f_baf, f_phase", "cell fraction from the depth (gain 2(2^LRR - 1), loss 2(1 - 2^LRR)), from the folded band deviation d (gain 4d/(1 - 2d), loss 4d/(1 + 2d), LOH 2d) and from the phased shift (the same formulas on |shift|)"),
    ("lrr, lrr_se, n_bins", "the segment's mean LRR, its standard error from the chromosome's noise scale, and its bin count"),
    ("d_hat, llr_baf", "band deviation at the segment's heterozygous sites (depth-aware maximum likelihood) and its log-likelihood ratio against d = 0"),
    ("phase_shift, phase_se, n_phased", "the phased fraction's shift from 1/2 over the event (maternal allele in the child; transmitted allele in a parent), its binomial error, and the site count"),
    ("origin, origin_llr, origin_n", "the child's parent of origin from the informative sites by likelihood ratio (maternal over paternal for a gain or LOH; paternal-lost over maternal-lost for a loss), and the site count"),
    ("origin_phase", "from the sign of the phased shift: for the child the parent of the extra, lost or retained copy (for UPD, of both copies); for a parent whether the duplicated, lost or retained homologue is the transmitted one (read folded: the transmitted homologue switches at the child's crossovers, and the statement covers the first run with the switch positions). Requires 40 sites and 3 standard errors"),
    ("homologues, hetero_share", "for a child's gain, LOH or UPD with a parent named and |shift| of 0.03 or more (a gain in about 13% of cells): whether that parent's two copies are one homologue (the auxiliary track follows the main one) or two (it departs: 1/3 against 2/3 along a trisomy), and the share of the event's windows in which they differ"),
    ("stage, centromere, n_crossovers, crossovers, crossover_states", "for a child's whole-chromosome gain or heterodisomy with a parent named: the meiotic stage from the state of the two copies nearest the centromere (heterodisomic: meiosis I; isodisomic with a heterodisomic segment elsewhere: meiosis II; isodisomic throughout: mitotic, or meiosis II without a crossover), the centromeric state itself, and the crossovers as changes of state along the chromosome (positions in Mb, and the state before and after each: hetero>iso or iso>hetero; one row per crossover in crossovers.tsv). For a parent's event: the child's crossovers within it, where the event's homologue changes between the transmitted and the untransmitted one (transmitted>untransmitted or the reverse)"),
    ("het_rate, het_rate_rel, n_het, n_called", "heterozygosity rate in the segment, its ratio to the member's own (or to the panel's), and the counts"),
    ("mie_rate", "Mendelian-error rate among the segment's confident sites: a constitutional deletion, an isodisomy or a heterodisomy produces errors at the informative sites; a trisomy or a mosaic event does not; a run of homozygosity has none"),
    ("inheritance", "child: inherited from the father / the mother (an event of the same type in that parent with reciprocal overlap of at least 50%; 'in a share of the parent's cells' when mosaic there), or new; parent: passed to the child or not, by the same rule"),
    ("external", "labels of the supplied events (--events) whose intersection with the segment covers at least half of the shorter of the two"),
    ("note", "annotations: 'a run of homozygosity ... not a uniparental disomy' (constitutional LOH without Mendelian errors); 'Mendelian errors at the informative sites: a uniparental isodisomy'; 'no heterozygous calls: in every cell'; 'joined across n bin(s) ...'; 'the phased bands read a share of x% against the depth's y%: the depth's call may be an artefact'; 'from the phased bands: ...' (how a phased-scan event was typed); 'a uniparental heterodisomy ...'; 'the phased bands and the opposite-homozygote sites name different parents'; 'a loss on a single X'; 'whole-chromosome X/Y: n copies against m expected for a reported male/female (47,XXY ...)' with, for the Y, the source of the estimate"),
]

ROWS = [
    ("1 · LRR", "log2 of each bin's trimmed-mean depth over the member's autosomal median, GC- and panel-corrected where applicable, with its step fit (black) and the calls as lines: gain at the segment's mean LRR (vermillion), loss (blue), copy-neutral LOH or heterodisomy at 0 (reddish purple). Brackets above: supplied events.",
     "a plateau away from 0 over several bins: log2(1 + f/2) for a gain, log2(1 - f/2) for a loss; a flat LRR under split bands indicates a copy-neutral event"),
    ("2 · BAF", "alt-allele read fraction at the member's heterozygous sites (subsampled). In the child, informative sites (parents opposite homozygotes) are coloured by the parent of the alt allele (pink maternal, orange paternal).",
     "two bands splitting from 1/2; under a gain the duplicated parent's alleles lie on the upper band, under a loss the retained parent's; an empty interval is a constitutional loss of heterozygosity"),
    ("3 · phased fraction", "the maternal-allele fraction along the child, the transmitted-allele fraction along a parent: each phased site (faint grey), the pooled windows (green; hollow where shared by two or more members) and the step fit (black); thin pink and orange lines are the auxiliary tracks.",
     "one band rather than two. %s A thin line departing from the main track means the child carries two different homologues of that parent there (meiotic error, heterodisomy); it rejoins at each crossover, marked by a dashed vertical line along a whole-chromosome event." % DIRECTION),
    ("4 · copies", "the LRR step fit's copy number (2 x 2^LRR) split by the fraction's step fit: maternal (pink) and paternal (orange) along the child; transmitted (vermillion) and untransmitted (sky blue) along a parent.",
     "a trisomy reads 2 and 1 with its parent named; a deletion 1 and 0; an isodisomy or heterodisomy 2 and 0; a mosaic event fractional values; a run of homozygosity is blank (no phased sites)"),
    ("5 · het rate", "heterozygous calls per confident call, per bin.", "zero under a constitutional loss of heterozygosity (isodisomy, run of homozygosity, deletion); raised where the child carries two different homologues of one parent"),
]

KINDS = [
    ("gain (one homologue duplicated), cell fraction f", "up by log2(1 + f/2)", "1/(2 + f) and (1 + f)/(2 + f)", "on the duplicated parent's band: (1 + f)/(2 + f) if maternal, 1/(2 + f) if paternal", "duplicated parent 1 + f, other parent 1", "unchanged", "none", "the parent of the extra copy"),
    ("loss of one copy, cell fraction f", "down by log2(1 - f/2)", "(1 - f)/(2 - f) and 1/(2 - f)", "on the retained parent's band: 1/(2 - f) if the paternal copy is lost", "lost parent 1 - f, other parent 1", "unchanged (f < 1); zero at f = 1", "at f = 1, at every informative site", "the parent of the lost copy"),
    ("copy-neutral LOH, cell fraction f", "flat", "(1 - f)/2 and (1 + f)/2", "(1 + f)/2 on the retained parent's side", "retained 1 + f, replaced 1 - f", "unchanged (f < 1)", "none (f < 1)", "the parent of the retained copy"),
    ("uniparental isodisomy (one homologue twice)", "flat", "none: no heterozygous calls", "1 (the child homozygous for that parent's allele)", "2 and 0", "zero", "at every informative site", "the parent of both copies"),
    ("uniparental heterodisomy (two different homologues)", "flat", "one band at 1/2: heterozygosity retained", "1 on the main track; the auxiliary track at 1/2", "2 and 0", "unchanged or raised", "at every informative site", "the parent of both copies; detected by the phased scan alone"),
    ("run of homozygosity (identical by descent)", "flat", "none", "no sites (the parents share the haplotype)", "blank", "zero", "none", "no parent of origin: not a disomy"),
    ("meiotic trisomy (two different homologues)", "up by log2(1.5)", "1/3 and 2/3", "2/3 on the main track; the auxiliary track at 1/3 where the two copies differ and 2/3 where a crossover has made them identical", "2 and 1", "raised where the copies differ", "none", "the parent of the two copies; 'two different homologues'; meiosis I when heterodisomic at the centromere, meiosis II when isodisomic there with heterodisomic segments distally"),
    ("mitotic or meiosis II trisomy (one homologue twice)", "up by log2(1.5)", "1/3 and 2/3", "2/3 on both tracks", "2 and 1", "unchanged", "none", "'one homologue throughout'; stage 'mitotic, or meiosis II without a crossover'"),
    ("47,XXY, maternal (a son)", "the X at 0, where a son's X sits at -1 (two copies)", "one band at 1/2: heterozygous wherever the mother is", "1 (every allele maternal); the auxiliary track at the father's sites at 1/2 where the two maternal X's differ", "maternal 2, paternal 0", "raised from 0 (a son's X has no heterozygous sites)", "at every informative site under the autosomal rules (no paternal allele)", "both X's maternal: meiosis I if heterodisomic at the centromere, else meiosis II or post-zygotic; karyotype 47,XXY(mat,MI)"),
    ("47,XXY, paternal (a son)", "the X at 0", "one band at 1/2", "1/2 (one maternal, one paternal X)", "maternal 1, paternal 1", "raised from 0", "none", "the extra X paternal: a paternal meiosis I error (X and Y transmitted together); 47,XXY(pat,MI)"),
    ("mosaic 46,XY/47,XXY, extra X paternal, cell fraction f", "the X at log2((1 + f)/2)", "f/(1 + f) and 1/(1 + f)", "1/(1 + f)", "maternal 1, paternal f", "raised in proportion", "under the hemizygous rules, in proportion to f", "the extra X paternal, f from the depth and the phased fraction; mos 47,XXY(pat)[f]/46,XY"),
    ("45,X, paternal X lost (a daughter)", "the X at -1, where a daughter's X sits at 0", "none: no heterozygous calls", "1 (every allele maternal)", "maternal 1, paternal 0", "zero", "none under the hemizygous rules", "the lost X paternal; 45,X(pat)"),
    ("47,XXX, maternal (a daughter)", "the X at log2(1.5)", "1/3 and 2/3", "2/3; the auxiliary track at 1/3 where the two maternal X's differ", "maternal 2, paternal 1", "raised where the two maternal copies differ", "none", "the extra X maternal, staged from the centromere like a trisomy; 47,XXX(mat,MI)"),
]

RULES = [
    ("Depth calls", "a segment of at least five bins with |mean LRR| of at least max(0.07, 3 SE), that is a gain or loss in about 10% of cells or more; bins the panel cannot characterise are excluded; adjacent pieces separated by at most three bins and differing in LRR by at most 0.08 are joined."),
    ("Band calls", "copy-neutral LOH where the bands split (d of at least 0.04 with a log-likelihood ratio of at least 10: about 8% of cells or more) or the heterozygosity rate falls to 0.35 of the member's own or below (constitutional), in bins the depth did not call; the heterozygosity-rate criterion takes precedence over a d estimated from few sites."),
    ("Phased scan", "on the main phased track's windows outside the called events: a segment of eight windows (or a whole short chromosome) with a mean shift of at least 0.015 (a gain or loss in about 6% of cells; 0.025, about 5%, for a copy-neutral event with flat depth) and five empirical standard errors, whose median shift and at least three quarters of whose windows have the same sign, spanning at least 2 Mb; typed as a gain or loss by the depth (three standard errors), otherwise as copy-neutral LOH; a shift of 0.4 or more with heterozygosity retained is a heterodisomy. Windows shared by two or more members are excluded first. Rejected segments are listed in phased_rejected.tsv with the reason."),
    ("Sign convention", DIRECTION),
    ("Three estimates of the cell fraction", "from the depth, from the folded bands and from the phased shift (formulas under the columns below). They should agree; a phased estimate under half the depth's, over 500 or more sites, flags the depth call as a possible artefact (GC, a panel edge)."),
    ("Two estimates of the parent of origin", "the likelihood ratio over the informative sites (origin, origin_llr) and the sign of the phased shift (origin_phase, over three to four times as many sites). A disagreement is noted."),
    ("One homologue or two", "for a child's gain, LOH or UPD with its parent named: the auxiliary track read at the other parent's homozygous sites follows the main track where the two copies are one homologue and departs from it where they are two; the share of windows in which it departs is hetero_share."),
    ("Meiotic stage", "for a child's whole-chromosome gain or heterodisomy with its parent named: the per-window state (two homologues or one) is smoothed by a running majority of five windows and must persist over five windows to count (a shorter run is noise, not a crossover); the state of the eight windows nearest the centromere, within 15 Mb of it, classifies the event: heterodisomic, a meiosis I nondisjunction; isodisomic with a heterodisomic segment elsewhere, meiosis II; isodisomic throughout, a mitotic duplication or a meiosis II error without a crossover. Each change of state along the chromosome is a crossover, placed between the two windows. For an event with both copies from one parent (a heterodisomy, a maternal XXY) the auxiliary track is read at any confident child call, where it sits at 0 or 1 if the two copies are one homologue and near 1/2 if they differ."),
    ("Run of homozygosity versus disomy", "a constitutional LOH without Mendelian errors is a run of homozygosity (both copies identical by descent; the parents share the haplotype); with errors at the informative sites it is a uniparental isodisomy. A heterodisomy retains heterozygosity and shows the errors."),
    ("Sex chromosomes", "each member's X and Y copy numbers come from the depth (raw and rounded) and form a complement (XY, XX, XXY, X, XYY, ...) checked against the pedigree sex. Against the complement the pedigree sex implies, a whole-X or whole-Y deviation of 10% or more is an event with its cell fraction (47,XXY, 46,XY/47,XXY, 45,X, 45,X/46,XX, 47,XXX, 47,XYY, mosaic loss of Y); the Y estimate comes from the panel's Y rows, else from the father/son Y depth ratio, else only beyond 25%. A male's X has a hemizygous baseline: a gain's extra X is maternal when the maternal fraction stays at 1 (staged from the centromere like a trisomy) and paternal when it falls to 1/(1+f) (a paternal meiosis I error); a daughter\'s paternal extra X (47,XXX) is the father\'s single X twice (meiosis II or post-zygotic). Mendelian errors on a single X are a heterozygous call or an allele the mother lacks. Segmental X events are read relative to the member\'s own X level; the X copy-neutral LOH search is skipped unless the member has exactly two X copies; pseudoautosomal bins are excluded; the Y is not segmented. Each member\'s karyotype string (summary and page) summarises the complement, the constitutional terms and each mosaic as its own line: 47,XXY(mat,MI), mos 47,XXY(pat)[0.40]/46,XY, 47,XY,+21mat(MI)."),
    ("Parents' events", "a parent's gain, loss or LOH is analysed the same way; 'passed to the child' when the child carries an event of the same type with reciprocal overlap of at least 50%; the sign of the phased shift shows whether it lies on the transmitted homologue."),
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
             "<p class=\"lead\">triokaryo detects large chromosomal events in a trio's VCF from three signals: <b>depth</b> (the LRR per bin; copy number), the "
             "<b>B-allele bands</b> at heterozygous sites (a second estimate of the cell fraction, and the copy-neutral events depth cannot show) and "
             "<b>transmission phasing</b> (the parental origin of each allele, giving a signed track that names the parent of origin, distinguishes one "
             "homologue from two, detects events below the depth limit and places boundaries at site resolution). This page defines every figure row, colour, "
             "call and column. Worked examples on two public trios: <a href=\"%s\">%s</a>.</p>" % (SITE, SITE),
             "<p class=\"toc\"><a href=\"#rows\">The rows</a><a href=\"#colours\">The colours</a><a href=\"#kinds\">What each event looks like</a><a href=\"#rules\">The rules</a>"
             "<a href=\"#columns\">The columns</a><a href=\"#cohort\">The cohort report</a></p>",
             "<h2 id=\"rows\">The rows of a figure</h2>",
             "<p>A chromosome figure has three columns (child, father, mother) and five rows. The genome figure has, per member, the first two rows (the pooled "
             "phased fraction drawn over the raw BAF), then the child's depth over the parents' mean (zero where the child's copy number equals the parental "
             "mean) and the child's maternal and paternal copy number along the genome.</p>",
             "<table><tr><th>row</th><th>what is drawn</th><th>what to look for</th></tr>"]
        for r, what, look in ROWS:
            w.append("<tr><td><b>%s</b></td><td>%s</td><td>%s</td></tr>" % (html.escape(r), html.escape(what), html.escape(look)))
        w.append("</table>")
        w.append("<h2 id=\"colours\">The colours</h2><p>Okabe-Ito throughout; a colour has the same meaning in every figure. The image carries no legend: "
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
        w.append("<p><code>summary.tsv</code>: the trio, records and sites used, the genome-wide Mendelian-error rate, per member the event counts by type, the X "
                 "copy number with the check against the pedigree sex, the median depth, the phased-site count, the shared windows and the phased-scan events; "
                 "the parameters. <code>bins.tsv</code>: every bin's depth, LRR, call counts, heterozygosity rate and band deviation per member, with the panel's "
                 "values and masks. <code>phased.tsv</code>: every window of the phased tracks (fraction, error, step fit, per-homologue copy number, auxiliary "
                 "tracks). <code>phased_rejected.tsv</code>: the phased scan's rejected segments, with the reason. <code>external.tsv</code>: the supplied events "
                 "and their matches. <code>events.bed</code> and <code>events.vcf</code>: the events as BED6 (score 1000 x cell fraction) and as a structural-variant VCF. "
                 "<code>crossovers.tsv</code>: one row per crossover of the staged events, with the state on each side.</p>")
        w.append("<h2 id=\"cohort\">The cohort report</h2><p>Written by <code>triokaryo cohort</code> over many runs: counts (trios; events by type, member and "
                 "source; de novo and inherited; phased-scan events; depth calls doubted by the phased track; X copy numbers disagreeing with the pedigree); the "
                 "landscape figure (events per chromosome, and one row per trio with each event a bar coloured by type: the child's thick, a parent's thin above "
                 "or below, a black line over a bar for an event matched by a supplied one); a table of every event (click a heading to sort, type in the box to "
                 "filter) linked to its trio's page and chromosome figure; per-trio quality metrics with each member's karyotype string; the sex-chromosome "
                 "aneuploidies (one row per whole-X or whole-Y event, with counts by complement, parent of origin and stage; sex_aneuploidies.tsv); the "
                 "concordance with the supplied events; the rejected phased segments by region; and this guide.</p>")
        w.append("<p style=\"color:#666;font-size:12px\">triokaryo %s</p></body></html>" % html.escape(__version__))
        with open(out_html, "w") as fh:
            fh.write("\n".join(w))
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    return out_html
