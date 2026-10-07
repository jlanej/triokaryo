"""Tables and the per-trio page. Every number of the figures is in a table; the page embeds the figures and their keys."""
import base64
import html
import os

import numpy as np

from .model import MEMBERS, NA

EVENT_COLS = ("sample", "role", "chrom", "start", "end", "span", "type", "source", "f", "f_lrr", "f_baf", "f_phase", "lrr", "lrr_se", "n_bins", "d_hat", "llr_baf",
              "phase_shift", "phase_se", "n_phased", "homologues", "hetero_share", "start_fine", "end_fine", "edge_sites", "het_rate", "het_rate_rel", "n_het", "n_called", "mie_rate", "origin",
              "origin_llr", "origin_n", "origin_phase", "inheritance", "external", "note")


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        return ("%." + str(nd) + "g") % x if np.isfinite(x) else "NA"
    return str(x)


def write_tsv(path, cols, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c)) for c in cols) + "\n")


PHASED_COLS = ("role", "chrom", "start", "end", "mid", "n_sites", "depth", "frac", "se", "shared", "step", "lrr", "copies_tagged", "copies_other",
               "aux_mother_hom", "aux_mother_hom_sites", "aux_father_hom", "aux_father_hom_sites", "aux_child_het", "aux_child_het_sites")


def write_tables(out, trio, bins, events, x_copies, scan, base_mie, params, tracks=None, phase_info=None, rejected=None):
    os.makedirs(out, exist_ok=True)
    if rejected is not None:
        write_tsv(os.path.join(out, "phased_rejected.tsv"), ["sample", "role", "chrom", "start", "end", "shift", "windows", "reason"], rejected)
    if tracks:
        rows = []
        for m, role in enumerate(MEMBERS):
            for c, t in tracks[m].items():
                for i in range(len(t.w_mid)):
                    r = dict(role=role, chrom=c, start=int(t.w_start[i]), end=int(t.w_end[i]), mid=int(t.w_mid[i]), n_sites=int(t.w_n[i]), depth=int(t.w_dp[i]),
                             frac=t.w_frac[i], se=t.w_se[i], shared=bool(t.w_shared[i]) if t.w_shared is not None else False, step=t.w_step[i], lrr=t.w_lrr[i],
                             copies_tagged=t.w_copies_tag[i], copies_other=t.w_copies_other[i])
                    for k, (af, an) in t.aux.items():
                        r["aux_%s" % k] = af[i]
                        r["aux_%s_sites" % k] = int(an[i])
                    rows.append(r)
        write_tsv(os.path.join(out, "phased.tsv"), list(PHASED_COLS), rows)
    rows = []
    for i in range(bins.n):
        r = dict(chrom=bins.chrom[i], start=int(bins.start[i]), end=int(bins.end[i]), gc=bins.gc[i], child_vs_mid=bins.child_vs_mid[i],
                 father_vs_mother=bins.father_vs_mother[i], panel_median=bins.panel_median[i] if bins.panel_median is not None else NA,
                 panel_rsd=bins.panel_rsd[i] if bins.panel_rsd is not None else NA, masked=bool(bins.masked[i]) if bins.masked is not None else False,
                 panel_bdev=bins.panel_bdev[i] if bins.panel_bdev is not None else NA, masked_bands=bool(bins.masked_bands[i]) if bins.masked_bands is not None else False)
        for m, role in enumerate(MEMBERS):
            r.update({"%s_n_sites" % role: int(bins.n_dp[m, i]), "%s_depth" % role: bins.depth[m, i], "%s_lrr" % role: bins.lrr[m, i],
                      "%s_lrr_gc" % role: bins.lrr_gc[m, i], "%s_n_called" % role: int(bins.n_called[m, i]), "%s_n_het" % role: int(bins.n_het[m, i]),
                      "%s_het_rate" % role: bins.het_rate[m, i], "%s_bdev" % role: bins.bdev[m, i],
                      "%s_het_rel" % role: bins.het_rel[m, i] if bins.het_rel is not None else NA})
        rows.append(r)
    cols = ["chrom", "start", "end", "gc", "panel_median", "panel_rsd", "masked", "panel_bdev", "masked_bands"] + \
           ["%s_%s" % (role, k) for role in MEMBERS for k in ("n_sites", "depth", "lrr", "lrr_gc", "n_called", "n_het", "het_rate", "bdev", "het_rel")] + ["child_vs_mid", "father_vs_mother"]
    write_tsv(os.path.join(out, "bins.tsv"), cols, rows)
    write_tsv(os.path.join(out, "events.tsv"), list(EVENT_COLS), [e.as_dict() for e in events])
    summ = dict(trio=trio.name, child=trio.kid, father=trio.dad, mother=trio.mom, child_sex=trio.kid_sex, father_sex=trio.dad_sex, mother_sex=trio.mom_sex,
                records=scan.n_records, sites_used=scan.n_used, skipped="; ".join("%s %d" % kv for kv in sorted(scan.skipped.items())),
                mie_rate_genome=base_mie, bin_size=bins.bin_size, gc_corrected=bool(np.isfinite(bins.gc).any()),
                panel=bool(bins.panel_median is not None and np.isfinite(bins.panel_median).any()), bins_masked=int(bins.masked.sum()) if bins.masked is not None else 0)
    for m, role in enumerate(MEMBERS):
        ev = [e for e in events if e.role == role]
        summ["%s_events" % role] = len(ev)
        summ["%s_gains" % role] = sum(1 for e in ev if e.type == "gain")
        summ["%s_losses" % role] = sum(1 for e in ev if e.type == "loss")
        summ["%s_loh" % role] = sum(1 for e in ev if e.type == "LOH")
        xc = x_copies.get(role, NA)
        summ["%s_x_copies" % role] = xc
        sex = trio.sexes[m]
        summ["%s_x_check" % role] = ("" if not sex or not np.isfinite(xc) else "agrees" if (sex == "M" and xc == 1) or (sex == "F" and xc == 2)
                                     else "X copies %d in a reported %s" % (int(xc), "male (47,XXY or an XX male?)" if sex == "M" else "female (45,X?)"))
        summ["%s_depth" % role] = float(np.nanmedian(bins.depth[m][bins.autosomal])) if np.isfinite(bins.depth[m][bins.autosomal]).any() else NA
        if phase_info and role in phase_info:
            summ["%s_phased_sites" % role] = phase_info[role]["n_phased"]
            summ["%s_window_sites" % role] = phase_info[role]["window_sites"]
            summ["%s_ref_bias" % role] = phase_info[role]["ref_bias"]
            summ["%s_phased_finds" % role] = sum(1 for e in ev if e.source == "phased")
            summ["%s_windows_shared" % role] = phase_info[role].get("windows_shared", 0)
    summ.update({"param_" + k: v for k, v in params.items()})
    write_tsv(os.path.join(out, "summary.tsv"), list(summ.keys()), [summ])
    return summ


def _img(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def _table(cols, rows):
    h = "<table><tr>" + "".join("<th>%s</th>" % html.escape(c) for c in cols) + "</tr>"
    for r in rows:
        h += "<tr>" + "".join("<td>%s</td>" % html.escape(fmt(r.get(c), 3)) for c in cols) + "</tr>"
    return h + "</table>"


CSS = ("body{font-family:Helvetica,Arial,sans-serif;max-width:1200px;margin:1.5em auto;padding:0 1em;color:#222}table{border-collapse:collapse;font-size:12px;"
       "margin:0.5em 0}th,td{border:1px solid #ddd;padding:2px 6px;text-align:left}th{background:#f3f3f3}img{max-width:100%}figure{margin:1em 0}"
       "figcaption{font-size:12px;color:#444}.key{font-size:11px;color:#555}h2{margin-top:1.5em}.mock{background:#fff3cd;padding:0.5em;border:1px solid #e0c060}")


def write_html(out, trio, figs, events, summ, external, mock_note=""):
    w = []
    w.append('<!doctype html><html><head><meta charset="utf-8"><title>triokaryo %s</title><style>%s</style></head><body>' % (html.escape(trio.name), CSS))
    if mock_note:
        w.append('<div class="mock">%s</div>' % html.escape(mock_note))
    w.append("<h1>Trio %s: large chromosomal events from the VCF</h1>" % html.escape(trio.name))
    w.append("<p>child %s (%s), father %s (%s), mother %s (%s). %s sites used of %s records (%s). Genome-wide Mendelian-error rate %s. X copies: child %s %s, "
             "father %s %s, mother %s %s.</p>" % tuple(html.escape(str(v)) for v in (
                 trio.kid, trio.kid_sex or "sex not given", trio.dad, trio.dad_sex or "?", trio.mom, trio.mom_sex or "?", summ["sites_used"], summ["records"],
                 summ["skipped"], fmt(summ["mie_rate_genome"], 3), fmt(summ["child_x_copies"]), summ["child_x_check"], fmt(summ["father_x_copies"]),
                 summ["father_x_check"], fmt(summ["mother_x_copies"]), summ["mother_x_check"])))
    w.append("<h2>Events</h2>")
    if events:
        w.append(_table(["role", "chrom", "span", "start", "end", "start_fine", "end_fine", "type", "source", "f", "f_lrr", "f_baf", "f_phase", "lrr", "d_hat", "phase_shift",
                         "n_phased", "homologues", "het_rate_rel", "mie_rate", "origin", "origin_llr", "origin_n", "origin_phase", "inheritance", "external", "note"], [e.as_dict() for e in events]))
    else:
        w.append("<p>No event called.</p>")
    if external:
        w.append("<h3>Events given from elsewhere</h3>")
        w.append(_table(["sample", "chrom", "start", "end", "note", "inheritance"], [dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, note=x.note, inheritance=x.inheritance) for x in external]))
    for n, fg in enumerate(figs, 1):
        w.append('<figure><img src="%s" alt="%s"><figcaption><b>Figure %d. %s.</b> %s</figcaption><div class="key">%s</div></figure>' % (
            _img(fg["png"]), html.escape(fg["name"]), n, html.escape(fg["title"]), html.escape(fg["caption"]),
            " &middot; ".join('<span style="color:%s">&#9632;</span> %s' % (c, html.escape(t)) for c, _, t in (__import__("triokaryo.plots", fromlist=["KEY"]).KEY[k] for k in fg["keys"]))))
    w.append("<h2>How to read it</h2><p>LRR: log2 of the bin's median depth over the member's autosomal median; a gain of one copy in a share f of the cells "
             "reads log2(1 + f/2), a loss log2(1 - f/2). BAF: the alt-allele fraction at heterozygous sites; a gain parts the bands to 1/(2+f) and (1+f)/(2+f), "
             "a loss to (1-f)/(2-f) and 1/(2-f), a copy-neutral loss of heterozygosity to (1-f)/2 and (1+f)/2 - f estimated from the bands independently "
             "of the depth (f_baf beside f_lrr). The heterozygosity rate falls to zero under a loss of heterozygosity in every cell. Parent of origin: at "
             "sites where the parents are opposite homozygotes the child's alleles have known parents, and the alt fraction says whose copy is extra, lost "
             "or doubled. Phased: at every heterozygous site where at least one parent is homozygous the child's alleles have known parents (and at a parent's, "
             "where the child or the other parent is, the allele passed to the child is known); the fraction of the maternal allele along the child, of the "
             "transmitted allele along a parent, sits at 1/2 + d or 1/2 - d along an event - the sign is the parent of origin (phase_shift, f_phase, origin_phase), "
             "the pooled sites give the edges at site resolution (start_fine, end_fine), and the LRR's copies split by the fraction are the maternal and "
             "paternal copies. Source: depth (the LRR), bands (the folded bands or the heterozygosity rate), phased (a shift of the phased track the depth did "
             "not call: a few per cent of cells). Inheritance: the same event in a parent. Every number is in events.tsv, bins.tsv and phased.tsv beside this page.</p>")
    w.append("</body></html>")
    with open(os.path.join(out, "index.html"), "w") as fh:
        fh.write("\n".join(w))
