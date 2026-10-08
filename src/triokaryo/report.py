"""Tables, the per-trio page and the cohort report. Every quantity drawn in a figure is in a table; the pages embed the
figures with their captions and keys."""
import base64
import html
import json
import os

import numpy as np

from .model import MEMBERS, NA

EVENT_COLS = ("sample", "role", "chrom", "start", "end", "span", "bands", "type", "source", "f", "f_lrr", "f_baf", "f_phase", "lrr", "lrr_se", "n_bins", "d_hat", "llr_baf",
              "phase_shift", "phase_se", "n_phased", "homologues", "hetero_share", "stage", "centromere", "n_crossovers", "crossovers", "crossover_states", "start_fine", "end_fine", "edge_sites",
              "het_rate", "het_rate_rel", "n_het", "n_called", "mie_rate", "origin", "origin_llr", "origin_n", "origin_phase", "inheritance", "external", "note")


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


def write_bed(path, events, label="triokaryo", trio_of=None):
    """BED6 of the events for genome browsers and interval tools: name sample|role|type|bands|f|origin (spaces as underscores),
    score 1000 x the cell fraction. trio_of: an event -> trio name, prefixed to the name in a cohort file."""
    with open(path, "w") as fh:
        fh.write('track name="%s" description="large chromosomal events (triokaryo)" useScore=1\n' % label)
        for e in events:
            f = e.f if np.isfinite(e.f) else 0.0
            parts = [trio_of(e)] if trio_of else []
            parts += [e.sample, e.role, e.type, e.bands or e.span, "f=%.2f" % f, e.origin_phase or e.origin]
            name = "|".join(x for x in parts if x).replace(" ", "_")
            fh.write("%s\t%d\t%d\t%s\t%d\t.\n" % (e.chrom, e.start, e.end, name, int(round(1000 * min(f, 1.0)))))


def write_vcf(path, events, trio, genome):
    """The events as a structural-variant VCF (one record per event; ALT <DUP>, <DEL>, <CNLOH> or <UPD>; INFO END, SVLEN, SVTYPE,
    CF the cell fraction, SOURCE, ORIGIN, STAGE, BANDS, INHERITANCE, DOUBT for a depth call the phased track doubts) with the
    three members as samples: the carrier has GT 0/1 (1/1 for a constitutional UPD or LOH), CN the copy number implied by the
    cell fraction and CF; the others ./.. For tools that take CNV or SV VCFs (annotation, filtering, browsers)."""
    members = list(trio.members)
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=triokaryo\n")
        for c in genome.chroms:
            fh.write("##contig=<ID=%s,length=%d>\n" % (c, genome.length[c]))
        for alt, desc in (("DUP", "gain of one copy (whole chromosome, arm or segment)"), ("DEL", "loss of one copy"), ("CNLOH", "copy-neutral loss of heterozygosity"),
                          ("UPD", "uniparental heterodisomy")):
            fh.write("##ALT=<ID=%s,Description=\"%s\">\n" % (alt, desc))
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
                 '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">\n'
                 '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="DUP, DEL, CNLOH or UPD">\n'
                 '##INFO=<ID=CF,Number=1,Type=Float,Description="Cell fraction carrying the event">\n'
                 '##INFO=<ID=SOURCE,Number=1,Type=String,Description="depth, bands or phased">\n'
                 '##INFO=<ID=ORIGIN,Number=1,Type=String,Description="Parent of the extra, lost or retained copy (mat or pat) for the child">\n'
                 '##INFO=<ID=STAGE,Number=1,Type=String,Description="Meiotic stage of a whole-chromosome event (MI, MII, MII/mit)">\n'
                 '##INFO=<ID=BANDS,Number=1,Type=String,Description="Cytogenetic bands spanned">\n'
                 '##INFO=<ID=INHERITANCE,Number=1,Type=String,Description="inherited, new, passed or not_passed">\n'
                 '##INFO=<ID=DOUBT,Number=0,Type=Flag,Description="A depth call whose phased cell fraction is under half the depth one">\n'
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0/1 carrier, 1/1 constitutional UPD or LOH, ./. not the carrier">\n'
                 '##FORMAT=<ID=CN,Number=1,Type=Float,Description="Copy number implied by the cell fraction (2 for copy-neutral events)">\n'
                 '##FORMAT=<ID=CF,Number=1,Type=Float,Description="Cell fraction">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(members))
        from .sexchrom import _short_origin, _short_stage
        for e in sorted(events, key=lambda e: (genome.chroms.index(e.chrom) if e.chrom in genome.chroms else 99, e.start)):
            svtype = {"gain": "DUP", "loss": "DEL", "LOH": "CNLOH", "UPD": "UPD"}[e.type]
            f = e.f if np.isfinite(e.f) else 1.0
            cn = 2 + f if e.type == "gain" else 2 - f if e.type == "loss" else 2.0
            if e.chrom == "chrX" and e.span == "whole" and e.role != "mother" and "expected for a reported male" in e.note:
                cn = 1 + f if e.type == "gain" else 1 - f                       # a male's X: one copy expected
            if e.chrom == "chrY":
                cn = 1 + f if e.type == "gain" else 1 - f
            inh = "inherited" if e.inheritance.startswith("inherited") else "new" if e.inheritance.startswith("new") else "passed" if e.inheritance == "passed to the child" else "not_passed" if e.inheritance else ""
            info = ["END=%d" % e.end, "SVLEN=%d" % (e.end - e.start), "SVTYPE=%s" % svtype, "CF=%.3f" % f, "SOURCE=%s" % e.source]
            o, stg = _short_origin(e), _short_stage(e.stage)
            if o:
                info.append("ORIGIN=%s" % o)
            if stg:
                info.append("STAGE=%s" % stg.replace("/", "_"))
            if e.bands:
                info.append("BANDS=%s" % e.bands)
            if inh:
                info.append("INHERITANCE=%s" % inh)
            if "may be an artefact" in e.note:
                info.append("DOUBT")
            gt = "1/1" if (e.type in ("LOH", "UPD") and f >= 0.9) else "0/1"
            cols = ["%s:%.2f:%.3f" % (gt, cn, f) if m == e.sample else "./.:.:." for m in members]
            fh.write("%s\t%d\t%s\tN\t<%s>\t.\tPASS\t%s\tGT:CN:CF\t%s\n" % (e.chrom, e.start + 1, "%s_%s_%d_%s" % (e.sample, e.chrom, e.start + 1, svtype), svtype, ";".join(info), "\t".join(cols)))


CROSSOVER_COLS = ["trio", "sample", "role", "chrom", "position", "position_mb", "from_state", "to_state", "event_type", "event_f", "stage", "parent"]


def crossover_rows(events, trio=""):
    """One row per crossover of the events that carry them: the child's whole-chromosome gains and heterodisomies (the state of the
    two copies changes), and a parent's events spanning one of the child's crossovers (the event's homologue changes between the
    transmitted and the untransmitted one)."""
    out = []
    for e in events:
        if not e.crossovers:
            continue
        pos = [float(x) for x in e.crossovers.split(";") if x]
        states = (e.crossover_states or "").split(";")
        parent = "maternal" if "maternal" in (e.origin_phase or e.origin) else "paternal" if "paternal" in (e.origin_phase or e.origin) else ""
        for i, x in enumerate(pos):
            fr, to = (states[i].split(">") + ["", ""])[:2] if i < len(states) and states[i] else ("", "")
            out.append(dict(trio=trio, sample=e.sample, role=e.role, chrom=e.chrom, position=int(round(x * 1e6)), position_mb=x,
                            from_state={"hetero": "heterodisomic", "iso": "isodisomic"}.get(fr, fr), to_state={"hetero": "heterodisomic", "iso": "isodisomic"}.get(to, to),
                            event_type=e.type, event_f=e.f, stage=e.stage, parent=parent))
    return out


PHASED_COLS = ("role", "chrom", "start", "end", "mid", "n_sites", "depth", "frac", "se", "shared", "step", "lrr", "copies_tagged", "copies_other",
               "aux_mother_hom", "aux_mother_hom_sites", "aux_father_hom", "aux_father_hom_sites", "aux_child_het", "aux_child_het_sites")


def karyotypes_of(events, summ, genome):
    """Each member's ISCN-like karyotype string from its events and the complement in the summary (X copies alone when the
    complement is missing, as in older runs)."""
    from .sexchrom import karyotype_string
    order = {c: i for i, c in enumerate(genome.chroms)}
    out = {}
    for role in MEMBERS:
        comp = summ.get("%s_sex_karyotype" % role) or ""
        if not comp and _finite(summ.get("%s_x_copies" % role)):
            comp = "X" * int(float(summ["%s_x_copies" % role]))
        out[role] = karyotype_string([e for e in events if e.role == role], comp, order, genome)
    return out


def write_tables(out, trio, bins, events, x_copies, scan, base_mie, params, tracks=None, phase_info=None, rejected=None, sex_info=None, genome=None):
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
    write_bed(os.path.join(out, "events.bed"), events, label="triokaryo %s" % trio.name)
    if genome is not None:
        write_vcf(os.path.join(out, "events.vcf"), events, trio, genome)
    write_tsv(os.path.join(out, "crossovers.tsv"), CROSSOVER_COLS, crossover_rows(events, trio.name))
    summ = dict(trio=trio.name, child=trio.kid, father=trio.dad, mother=trio.mom, child_sex=trio.kid_sex, father_sex=trio.dad_sex, mother_sex=trio.mom_sex,
                records=scan.n_records, sites_used=scan.n_used, skipped="; ".join("%s %d" % kv for kv in sorted(scan.skipped.items())),
                mie_rate_genome=base_mie, bin_size=bins.bin_size, gc_corrected=bool(np.isfinite(bins.gc).any()),
                panel=bool(bins.panel_median is not None and np.isfinite(bins.panel_median).any()), bins_masked=int(bins.masked.sum()) if bins.masked is not None else 0,
                genotypes_from_pl=bool(getattr(scan, "genotypes_from_pl", False)), depth_sites=getattr(bins, "depth_sites", "all"))
    for m, role in enumerate(MEMBERS):
        ev = [e for e in events if e.role == role]
        summ["%s_events" % role] = len(ev)
        summ["%s_gains" % role] = sum(1 for e in ev if e.type == "gain")
        summ["%s_losses" % role] = sum(1 for e in ev if e.type == "loss")
        summ["%s_loh" % role] = sum(1 for e in ev if e.type == "LOH")
        xc = x_copies.get(role, NA)
        summ["%s_x_copies" % role] = xc
        hr = getattr(bins, "homref_depth_ratio", None)
        summ["%s_homref_depth_ratio" % role] = float(hr[m]) if hr is not None and m < len(hr) else NA
        sex = trio.sexes[m]
        st = (sex_info or {}).get("states", {}).get(role)
        if st:
            from .sexchrom import karyotype
            summ["%s_x_check" % role] = st.get("x_check", "")
            summ["%s_x_copies_raw" % role] = st["x_copies_raw"]
            summ["%s_y_copies" % role] = st["y_copies"]
            summ["%s_y_copies_raw" % role] = st["y_copies_raw"]
            summ["%s_sex_karyotype" % role] = karyotype(st["x_copies"], st["y_copies"])
            summ["%s_sex_check" % role] = st.get("sex_check", "")
        else:
            summ["%s_x_check" % role] = ("" if not sex or not np.isfinite(xc) else "agrees" if (sex == "M" and xc == 1) or (sex == "F" and xc == 2)
                                         else "X copies %d in a reported %s" % (int(xc), "male (47,XXY or an XX male?)" if sex == "M" else "female (45,X?)"))
        summ["%s_depth" % role] = float(np.nanmedian(bins.depth[m][bins.autosomal])) if np.isfinite(bins.depth[m][bins.autosomal]).any() else NA
        if phase_info and role in phase_info:
            summ["%s_phased_sites" % role] = phase_info[role]["n_phased"]
            summ["%s_window_sites" % role] = phase_info[role]["window_sites"]
            summ["%s_ref_bias" % role] = phase_info[role]["ref_bias"]
            summ["%s_phased_finds" % role] = sum(1 for e in ev if e.source == "phased")
            summ["%s_windows_shared" % role] = phase_info[role].get("windows_shared", 0)
    if sex_info:
        summ["y_father_son_log2"] = sex_info.get("y_father_son_log2", NA)
        summ["y_father_son_sites"] = sex_info.get("y_father_son_sites", 0)
        summ["y_panel"] = bool(sex_info.get("y_panel", False))
        summ["x_offset_trio"] = sex_info.get("x_offset_trio", 0.0)
    if genome is not None:
        for role, k in karyotypes_of(events, summ, genome).items():
            summ["%s_karyotype" % role] = k
    summ.update({"param_" + k: v for k, v in params.items()})
    write_tsv(os.path.join(out, "summary.tsv"), list(summ.keys()), [summ])
    return summ


def _finite(v):
    """True for a number (or numeric string) that is finite: the page is also built from summary.tsv's strings."""
    try:
        return v is not None and v != "" and np.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _img(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def _table(cols, rows):
    h = "<table><tr>" + "".join("<th>%s</th>" % html.escape(c) for c in cols) + "</tr>"
    for r in rows:
        h += "<tr>" + "".join("<td>%s</td>" % html.escape(fmt(r.get(c), 3)) for c in cols) + "</tr>"
    return h + "</table>"


CSS = ("html{background:#fff;color-scheme:light}body{font-family:Helvetica,Arial,sans-serif;max-width:1200px;margin:1.5em auto;padding:0 1em;color:#222;background:#fff}table{border-collapse:collapse;font-size:12px;"
       "margin:0.5em 0}th,td{border:1px solid #ddd;padding:2px 6px;text-align:left}th{background:#f3f3f3}img{max-width:100%}figure{margin:1em 0}"
       "figcaption{font-size:12px;color:#444}.key{font-size:11px;color:#555}h2{margin-top:1.5em}.mock{background:#fff3cd;padding:0.5em;border:1px solid #e0c060}")


def write_html(out, trio, figs, events, summ, external, mock_note="", genome=None):
    w = []
    w.append('<!doctype html><html><head><meta charset="utf-8"><title>triokaryo %s</title><style>%s</style></head><body>' % (html.escape(trio.name), CSS))
    if mock_note:
        w.append('<div class="mock">%s</div>' % html.escape(mock_note))
    w.append("<h1>Trio %s: large chromosomal events from the VCF</h1>" % html.escape(trio.name))
    if genome is not None:
        kar = karyotypes_of(events, summ, genome)
        w.append("<p><b>Karyotype.</b> %s.</p>" % "; ".join("%s <code>%s</code>" % (role, html.escape(kar[role])) for role in MEMBERS))
    w.append("<p>child %s (%s), father %s (%s), mother %s (%s). %s sites used of %s records (%s). Genome-wide Mendelian-error rate %s.</p>" % tuple(
        html.escape(str(v)) for v in (trio.kid, trio.kid_sex or "sex not given", trio.dad, trio.dad_sex or "?", trio.mom, trio.mom_sex or "?", summ["sites_used"],
                                      summ["records"], summ["skipped"], fmt(summ["mie_rate_genome"], 3))))
    w.append("<p><b>Sex chromosomes.</b> %s.%s</p>" % (
        "; ".join(html.escape("%s %s (X %s, Y %s copies: %s)" % (
            role, summ.get("%s_sex_karyotype" % role) or ("X" * int(float(summ["%s_x_copies" % role])) if _finite(summ.get("%s_x_copies" % role)) else "?"),
            fmt(summ.get("%s_x_copies_raw" % role, NA), 3), fmt(summ.get("%s_y_copies_raw" % role, NA), 3),
            summ.get("%s_sex_check" % role) or summ.get("%s_x_check" % role) or "no pedigree sex")) for role in MEMBERS),
        (" Y depth, father over son (log2, each relative to its autosomes): %s over %s sites." % (fmt(summ["y_father_son_log2"], 3), summ["y_father_son_sites"]))
        if _finite(summ.get("y_father_son_log2")) else ""))
    w.append("<h2>Events</h2>")
    if events:
        w.append(_table(["role", "chrom", "span", "bands", "start", "end", "start_fine", "end_fine", "type", "source", "f", "f_lrr", "f_baf", "f_phase", "lrr", "d_hat", "phase_shift",
                         "n_phased", "homologues", "stage", "crossovers", "het_rate_rel", "mie_rate", "origin", "origin_llr", "origin_n", "origin_phase", "inheritance", "external", "note"],
                        [e.as_dict() for e in events]))
    else:
        w.append("<p>No event called.</p>")
    if external:
        w.append("<h3>Supplied events (--events)</h3>")
        w.append(_table(["sample", "chrom", "start", "end", "label", "match"], [dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in external]))
    for n, fg in enumerate(figs, 1):
        w.append('<figure><img src="%s" alt="%s"><figcaption><b>Figure %d. %s.</b> %s</figcaption><div class="key">%s</div></figure>' % (
            _img(fg["png"]), html.escape(fg["name"]), n, html.escape(fg["title"]), html.escape(fg["caption"]),
            " &middot; ".join('<span style="color:%s">&#9632;</span> %s' % (c, html.escape(t)) for c, _, t in (__import__("triokaryo.plots", fromlist=["KEY"]).KEY[k] for k in fg["keys"]))))
    flagged = [e for e in events if e.note]
    if flagged:
        w.append("<h3>Flags</h3><ul>" + "".join("<li>%s %s %s %s: %s</li>" % tuple(html.escape(str(v)) for v in (e.role, e.chrom, e.type, e.span, e.note)) for e in flagged) + "</ul>")
    w.append("<h2>How to read it</h2>")
    from .guide import key_table
    from .plots import DIRECTION
    w.append("<p><b>Rows.</b> LRR: log2 of the bin's trimmed-mean depth over the member's autosomal median, with its step fit and the calls (a gain in a cell fraction f "
             "reads log2(1 + f/2), a loss log2(1 - f/2); copy-neutral events are drawn at 0). BAF: the alt-allele read fraction at heterozygous sites, with the "
             "child's informative sites coloured by the parent of the alt allele. Phased fraction: the maternal-allele fraction along the child, the "
             "transmitted-allele fraction along a parent (sites, pooled windows, step fit); thin lines are the auxiliary tracks, which depart from the main track "
             "where the child carries two different homologues of one parent. Copies: the LRR step fit's copy number split by the fraction's step fit. Het rate: "
             "heterozygous calls per confident call.</p>")
    w.append("<p><b>Sign convention.</b> %s</p>" % html.escape(DIRECTION))
    w.append("<p><b>Calls.</b> source: depth (LRR segmentation), bands (band deviation or heterozygosity rate), phased (a shift of the phased track the depth did not "
             "call: a gain or loss in a few per cent of cells, or a uniparental heterodisomy). f: the cell fraction, from the depth for a depth-called gain or loss, "
             "otherwise from the bands; f_lrr, f_baf and f_phase are the three estimates side by side. origin: the parent of origin from the informative sites; "
             "origin_phase: from the sign of the phased shift; homologues: one or two. A constitutional LOH without Mendelian errors is a run of homozygosity, "
             "with them a uniparental isodisomy. Every quantity is in events.tsv, bins.tsv and phased.tsv beside this page; the full guide with pattern cards is "
             "<a href=\"guide.html\">guide.html</a>.</p>")
    used = []
    for fg in figs:
        for k in fg["keys"]:
            if k not in used:
                used.append(k)
    if used:
        w.append("<p><b>Colours.</b></p>" + key_table(used))
    w.append("</body></html>")
    with open(os.path.join(out, "index.html"), "w") as fh:
        fh.write("\n".join(w))


def _trio_of(summary):
    from .pedigree import Trio
    m, sx = summary["members"], summary.get("sexes") or ["", "", ""]
    return Trio(m[0], m[1], m[2], sx[0] or "", sx[1] or "", sx[2] or "")


def events_of(summary):
    """The events of a run's summary.json back as Event objects."""
    from .segment import Event
    out = []
    for e in summary["events"]:
        ev = Event(e["sample"], e["role"], e["chrom"], int(e["start"]), int(e["end"]), e["span"], e["type"])
        for k, v in e.items():
            if hasattr(ev, k) and k not in ("f",):
                if isinstance(v, str) or v is None:
                    setattr(ev, k, v if v is not None else (float("nan") if isinstance(getattr(ev, k), float) else ""))
                else:
                    setattr(ev, k, v)
        ev.note = e.get("note") or ""
        ev.external = e.get("external") or ""
        out.append(ev)
    return out


def _read_summary_tsv(path):
    if not os.path.exists(path):
        return {}
    lines = open(path).read().splitlines()
    if len(lines) < 2:
        return {}
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


def _externals_of(path):
    """The supplied events of a run, read back from external.tsv (sample, chrom, start, end, label, match)."""
    from .segment import Event
    out = []
    if not os.path.exists(path):
        return out
    head = None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if head is None:
            head = f
            continue
        r = dict(zip(head, f))
        x = Event(r["sample"], "", r["chrom"], int(float(r["start"])), int(float(r["end"])), "", "", note=r.get("label", ""))
        x.inheritance = r.get("match", "")
        out.append(x)
    return out


def rebuild_run(run_dir, log=None):
    """Rebuild a run's page, figure sidecars, legend images and guide from its tables and existing figures; no VCF needed."""
    from .genome import genome as load_genome
    from .guide import write_guide
    from .plots import chrom_caption, genome_caption, read_sidecar, write_legend, write_sidecar
    summary = json.load(open(os.path.join(run_dir, "summary.json")))
    trio = _trio_of(summary)
    events = events_of(summary)
    summ = _read_summary_tsv(os.path.join(run_dir, "summary.tsv"))
    for k in ("mie_rate_genome", "child_x_copies", "father_x_copies", "mother_x_copies", "y_father_son_log2"):
        try:
            summ[k] = float(summ.get(k, "nan"))
        except ValueError:
            summ[k] = float("nan")
    summ.setdefault("y_father_son_sites", "NA")
    try:
        bin_size = int(float(summ.get("bin_size", "0") or 0))          # the captions are rebuilt from the tables when the bin size is known
    except ValueError:
        bin_size = 0
    gc_corrected = str(summ.get("gc_corrected", "0")).lower() in ("1", "true")
    external = _externals_of(os.path.join(run_dir, "external.tsv"))
    G = load_genome(summary.get("genome", "grch38"))
    fdir = os.path.join(run_dir, "figures")
    figs = []
    if os.path.isdir(fdir):
        names = [f[:-4] for f in os.listdir(fdir) if f.endswith(".png") and (f == "genome.png" or f.startswith("chrom_"))]
        names.sort(key=lambda n: (n != "genome", G.chroms.index(n[6:]) if n.startswith("chrom_") and n[6:] in G.chroms else 99))
        for n in names:
            sc = os.path.join(fdir, n + ".txt")
            title, caption, keys = read_sidecar(sc) if os.path.exists(sc) else ("", "", None)
            if keys is None:
                from .plots import CHROM_KEYS, GENOME_KEYS
                keys = GENOME_KEYS if n == "genome" else CHROM_KEYS
            if bin_size and n == "genome":
                title, caption = title or "Large chromosomal events in trio %s, genome-wide" % trio.name, genome_caption(bin_size, gc_corrected)
            elif bin_size and n.startswith("chrom_"):
                title, caption = title or "Trio %s, %s" % (trio.name, n[6:]), chrom_caption(n[6:], trio.name, bin_size, events)
            write_sidecar(fdir, n, title, caption, keys)
            write_legend(fdir, n, keys)
            figs.append(dict(name=n, title=title, caption=caption, keys=keys, png=os.path.join(fdir, n + ".png")))
    write_html(run_dir, trio, figs, events, summ, external, summary.get("mock_note", ""), genome=G)
    write_guide(os.path.join(run_dir, "guide.html"))
    if log:
        log("%s: page, %d figure sidecars and legends, guide rebuilt" % (run_dir, len(figs)))
    return figs


JS = """
function sortTable(t, n) {
  const tb = t.tBodies[0], rows = Array.from(tb.rows), asc = !(t.dataset.sortCol == n && t.dataset.sortAsc == '1');
  const num = v => { const x = parseFloat(v.replace(/,/g, '')); return isNaN(x) ? null : x; };
  rows.sort((a, b) => { const va = a.cells[n].textContent.trim(), vb = b.cells[n].textContent.trim(); const na = num(va), nb = num(vb);
    const c = (na !== null && nb !== null) ? na - nb : va.localeCompare(vb); return asc ? c : -c; });
  rows.forEach(r => tb.appendChild(r)); t.dataset.sortCol = n; t.dataset.sortAsc = asc ? '1' : '0';
}
document.querySelectorAll('table.sortable').forEach(t => Array.from(t.tHead.rows[0].cells).forEach((th, i) => { th.style.cursor = 'pointer'; th.title = 'sort'; th.addEventListener('click', () => sortTable(t, i)); }));
const box = document.getElementById('filter');
if (box) box.addEventListener('input', () => { const q = box.value.toLowerCase(); let n = 0;
  document.querySelectorAll('#events tbody tr').forEach(r => { const on = r.textContent.toLowerCase().includes(q); r.style.display = on ? '' : 'none'; if (on) n++; });
  document.getElementById('nshown').textContent = n; });
"""

COHORT_COLS = ["trio", "sample", "role", "chrom", "start", "end", "start_fine", "end_fine", "span", "bands", "type", "source", "f", "f_lrr", "f_baf", "f_phase", "lrr", "d_hat",
               "phase_shift", "n_phased", "homologues", "hetero_share", "stage", "centromere", "n_crossovers", "crossovers", "crossover_states", "het_rate_rel", "mie_rate", "origin", "origin_llr",
               "origin_n", "origin_phase", "inheritance", "external", "note"]


def write_cohort(out, run_dirs, events_path=None, genome_name="grch38", log=None):
    """The cohort report over every trio's run: the tables (events.all.tsv, summary.all.tsv, concordance.tsv, rejected.all.tsv,
    flags.tsv), the landscape figure, the page (index.html) and the guide (guide.html)."""
    from .external import match_external, read_events
    from .genome import genome as load_genome
    from .guide import key_table, write_guide
    from .plots import DIRECTION, LANDSCAPE_KEYS, fig_landscape
    os.makedirs(out, exist_ok=True)
    G = load_genome(genome_name)
    summaries, events, rows, rejected = [], [], [], []
    for d in run_dirs:
        sj = os.path.join(d, "summary.json")
        if not os.path.exists(sj):
            continue
        s = json.load(open(sj))
        s["run"] = d
        s["tsv"] = _read_summary_tsv(os.path.join(d, "summary.tsv"))
        s["events_obj"] = events_of(s)
        summaries.append(s)
        for e in s["events_obj"]:
            events.append((s["trio"], e))
        for e in s["events"]:
            rows.append(dict(e, trio=s["trio"]))
        rp = os.path.join(d, "phased_rejected.tsv")
        if os.path.exists(rp):
            head = None
            for line in open(rp):
                f = line.rstrip("\n").split("\t")
                if head is None:
                    head = f
                    continue
                rejected.append(dict(zip(head, f), trio=s["trio"]))
    if not summaries:
        raise SystemExit("no run with a summary.json under: " + " ".join(run_dirs))
    write_tsv(os.path.join(out, "events.all.tsv"), COHORT_COLS, rows)
    trio_by_sample = {e.sample: t for t, e in events}
    write_bed(os.path.join(out, "events.all.bed"), [e for _, e in events], label="triokaryo cohort", trio_of=lambda e: trio_by_sample.get(e.sample, ""))
    write_tsv(os.path.join(out, "crossovers.all.tsv"), CROSSOVER_COLS, [r for t, e in events for r in crossover_rows([e], t)])
    # the concordance with the supplied events
    conc, ext = None, []
    if events_path:
        ext = read_events(events_path, G)
        samples = {m for s in summaries for m in s["members"]}
        ext = [x for x in ext if x.sample in samples]
        match_external([e for _, e in events], ext)
        conc = dict(external=len(ext), matched=sum(1 for x in ext if x.inheritance.startswith("matched")),
                    new=sum(1 for _, e in events if e.type in ("gain", "loss") and not e.external), cn=sum(1 for _, e in events if e.type in ("LOH", "UPD")))
        write_tsv(os.path.join(out, "concordance.tsv"), ["sample", "chrom", "start", "end", "label", "match"],
                  [dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in ext])
    # per trio
    scols = ["trio", "members", "sexes", "sites_used", "records", "mie_rate_genome", "x_copies", "sex_karyotypes", "child_sex_check", "father_sex_check", "mother_sex_check",
             "child_karyotype", "father_karyotype", "mother_karyotype", "y_father_son_log2", "child_depth", "father_depth", "mother_depth", "child_phased_sites",
             "father_phased_sites", "mother_phased_sites", "windows_shared", "rejected", "n_events", "n_flagged", "seconds", "run"]
    srows = []
    for s in summaries:
        t = s["tsv"]
        chk = {r: t.get("%s_sex_check" % r) or t.get("%s_x_check" % r, "") for r in MEMBERS}
        kar = karyotypes_of(s["events_obj"], t, G)
        srows.append(dict(trio=s["trio"], members=",".join(s["members"]), sexes=",".join(s.get("sexes") or []), sites_used=s["sites_used"], records=s["records"],
                          mie_rate_genome=s["mie_rate_genome"], x_copies=json.dumps(s["x_copies"]),
                          child_karyotype=kar["child"], father_karyotype=kar["father"], mother_karyotype=kar["mother"],
                          sex_karyotypes=",".join(t.get("%s_sex_karyotype" % r) or ("X" * int(float(t["%s_x_copies" % r])) if _finite(t.get("%s_x_copies" % r)) else "?") for r in MEMBERS),
                          child_sex_check=chk["child"], father_sex_check=chk["father"], mother_sex_check=chk["mother"], y_father_son_log2=t.get("y_father_son_log2", "NA"),
                          child_depth=t.get("child_depth", "NA"), father_depth=t.get("father_depth", "NA"), mother_depth=t.get("mother_depth", "NA"),
                          child_phased_sites=t.get("child_phased_sites", "NA"), father_phased_sites=t.get("father_phased_sites", "NA"), mother_phased_sites=t.get("mother_phased_sites", "NA"),
                          windows_shared=t.get("child_windows_shared", "NA"), rejected=sum(1 for r in rejected if r["trio"] == s["trio"]), n_events=len(s["events"]),
                          n_flagged=sum(1 for e in s["events_obj"] if e.note), seconds=s.get("seconds"), run=s["run"]))
    write_tsv(os.path.join(out, "summary.all.tsv"), scols, srows)
    write_tsv(os.path.join(out, "rejected.all.tsv"), ["trio", "sample", "role", "chrom", "start", "end", "shift", "windows", "reason"], rejected)
    flags = [dict(trio=t, sample=e.sample, role=e.role, chrom=e.chrom, start=e.start, end=e.end, type=e.type, source=e.source, f=e.f, note=e.note) for t, e in events if e.note]
    write_tsv(os.path.join(out, "flags.tsv"), ["trio", "sample", "role", "chrom", "start", "end", "type", "source", "f", "note"], flags)
    # the sex-chromosome aneuploidies: one row per whole-X or whole-Y event, with the member's complement and karyotype string
    import re
    kar_by_trio = {s["trio"]: karyotypes_of(s["events_obj"], s["tsv"], G) for s in summaries}
    comp_by_trio = {s["trio"]: {r: s["tsv"].get("%s_sex_karyotype" % r, "") for r in MEMBERS} for s in summaries}
    sex_ane = []
    for t, e in events:
        if e.chrom not in ("chrX", "chrY") or e.span != "whole":
            continue
        m = re.search(r"\(([^)]*)\)", e.note)
        sex_ane.append(dict(trio=t, sample=e.sample, role=e.role, complement=comp_by_trio[t].get(e.role, ""), karyotype=kar_by_trio[t].get(e.role, ""),
                            label=m.group(1) if m else "%s %s" % (e.chrom, e.type), chrom=e.chrom, type=e.type, f=e.f, origin=e.origin_phase or e.origin, stage=e.stage,
                            centromere=e.centromere, n_crossovers=e.n_crossovers, crossovers=e.crossovers, mie_rate=e.mie_rate, note=e.note))
    write_tsv(os.path.join(out, "sex_aneuploidies.tsv"), ["trio", "sample", "role", "complement", "karyotype", "label", "chrom", "type", "f", "origin", "stage", "centromere",
                                                           "n_crossovers", "crossovers", "mie_rate", "note"], sex_ane)
    # the figure and the guide
    fdir = os.path.join(out, "figures")
    os.makedirs(fdir, exist_ok=True)
    land = fig_landscape(events, G, fdir)
    write_guide(os.path.join(out, "guide.html"), figures_dir=fdir)
    # the counts
    kids = [e for _, e in events if e.role == "child"]
    n_by = lambda key, vals: {v: sum(1 for _, e in events if getattr(e, key) == v) for v in vals}  # noqa: E731
    by_type = n_by("type", ("gain", "loss", "LOH", "UPD"))
    by_role = n_by("role", ("child", "father", "mother"))
    by_src = n_by("source", ("depth", "bands", "phased"))
    new = sum(1 for e in kids if e.inheritance.startswith("new"))
    inh = sum(1 for e in kids if e.inheritance.startswith("inherited"))
    doubted = sum(1 for _, e in events if "may be an artefact" in e.note)
    roh = sum(1 for _, e in events if "run of homozygosity" in e.note)
    xbad = [(s["trio"], r) for s in summaries for r in MEMBERS
            if (s["tsv"].get("%s_sex_check" % r) or s["tsv"].get("%s_x_check" % r, "")) not in ("", "agrees", "agrees (Y not in the VCF)")]
    quiet = sum(1 for s in summaries if not s["events"])
    # the page
    w = ['<!doctype html><html><head><meta charset="utf-8"><title>triokaryo cohort</title><style>%s'
         'input{font-size:13px;padding:3px 6px;width:22em}.tiles{display:flex;flex-wrap:wrap;gap:10px;margin:0.8em 0}.tile{background:#f6f6f8;border:1px solid #ddd;border-radius:8px;'
         'padding:8px 12px;min-width:140px}.tile b{display:block;font-size:20px}.tile span{font-size:12px;color:#555}</style></head><body>' % CSS,
         "<h1>triokaryo: %d trios, %d events</h1>" % (len(summaries), len(events)),
         "<p>Large chromosomal events in each trio from its VCF, from depth, B-allele bands and transmission phasing. The meaning of every figure, colour and "
         "column: <a href=\"guide.html\">guide.html</a>. Tables beside this page: events.all.tsv, summary.all.tsv, concordance.tsv, flags.tsv, rejected.all.tsv.</p>",
         '<div class="tiles">']
    tiles = [(len(summaries), "trios (%d without an event)" % quiet), (len(events), "events"), (by_type["gain"], "gains"), (by_type["loss"], "losses"),
             (by_type["LOH"], "copy-neutral LOH"), (by_type["UPD"], "heterodisomies"), (by_role["child"], "in children (%d de novo, %d inherited)" % (new, inh)),
             (by_role["father"] + by_role["mother"], "in parents"), (by_src["phased"], "from the phased scan alone"), (doubted, "depth calls doubted by the phased track"),
             (roh, "runs of homozygosity"), (len(sex_ane), "sex-chromosome aneuploidies"), (len(xbad), "sex-chromosome complements disagreeing with the pedigree sex")]
    if conc:
        tiles.append((conc["matched"], "of %d supplied events matched" % conc["external"]))
    for n, label in tiles:
        w.append('<div class="tile"><b>%s</b><span>%s</span></div>' % (n, html.escape(label)))
    w.append("</div>")
    w.append('<figure><img src="%s" alt="landscape"><figcaption><b>Figure 1. %s.</b> %s</figcaption><div class="key">%s</div></figure>' % (
        _img(land["png"]), html.escape(land["title"]), html.escape(land["caption"]),
        " &middot; ".join('<span style="color:%s">&#9632;</span> %s' % (c, html.escape(t)) for c, _, t in (__import__("triokaryo.plots", fromlist=["KEY"]).KEY[k] for k in LANDSCAPE_KEYS))))
    w.append("<h2>Events</h2><p>Click a heading to sort; type to filter. <input id=\"filter\" placeholder=\"filter: a trio, a chromosome, a type, a word of a note\"> "
             "<span id=\"nshown\">%d</span> shown. %s</p>" % (len(events), html.escape(DIRECTION)))
    w.append('<table class="sortable" id="events"><thead><tr>' + "".join("<th>%s</th>" % h for h in (
        "trio", "member", "chrom", "start (Mb)", "end (Mb)", "span", "bands", "type", "source", "f", "f depth", "f bands", "f phased", "origin", "homologues", "stage", "crossovers (Mb)",
        "inheritance", "supplied", "notes", "page", "figure")) + "</tr></thead><tbody>")
    for s in summaries:
        page = os.path.relpath(os.path.join(s["run"], "index.html"), out)
        for e in s["events_obj"]:
            fine = lambda v, b: ("%.2f" % (v / 1e6)) if (isinstance(v, (int, float)) and np.isfinite(v)) else ("%.0f" % (b / 1e6))  # noqa: E731
            figp = os.path.relpath(os.path.join(s["run"], "figures", "chrom_%s.png" % e.chrom), out)
            cells = [s["trio"], "%s (%s)" % (e.sample, e.role), e.chrom, fine(e.start_fine, e.start), fine(e.end_fine, e.end), e.span, e.bands, e.type, e.source, fmt(e.f, 2), fmt(e.f_lrr, 2),
                     fmt(e.f_baf, 2), fmt(e.f_phase, 2), e.origin_phase or e.origin, e.homologues, e.stage, e.crossovers, e.inheritance, e.external, e.note]
            w.append("<tr>" + "".join("<td>%s</td>" % html.escape(str(c)) for c in cells) + '<td><a href="%s">page</a></td><td><a href="%s">%s</a></td></tr>' % (
                html.escape(page), html.escape(figp), html.escape(e.chrom)))
    w.append("</tbody></table>")
    w.append("<h2>Trios</h2><table class=\"sortable\"><thead><tr>" + "".join("<th>%s</th>" % h for h in (
        "trio", "members", "sexes", "events", "flagged", "karyotype (child)", "sex chromosomes (child, father, mother)", "against the pedigree sex", "Y father/son (log2)",
        "MIE rate", "sites", "depth (child, father, mother)", "phased sites (child)", "rejected", "page")) + "</tr></thead><tbody>")
    for s, r in zip(summaries, srows):
        xchk = "; ".join("%s: %s" % (m, r["%s_sex_check" % m]) for m in MEMBERS if r["%s_sex_check" % m] not in ("", "agrees", "agrees (Y not in the VCF)")) or "agrees"
        w.append("<tr>" + "".join("<td>%s</td>" % html.escape(str(c)) for c in (
            s["trio"], ", ".join(s["members"]), ",".join(s.get("sexes") or []), r["n_events"], r["n_flagged"], r["child_karyotype"], r["sex_karyotypes"].replace(",", ", "), xchk,
            r["y_father_son_log2"], fmt(s["mie_rate_genome"], 3), s["sites_used"], "%s, %s, %s" % (r["child_depth"], r["father_depth"], r["mother_depth"]), r["child_phased_sites"],
            r["rejected"])) + '<td><a href="%s">page</a></td></tr>' % html.escape(os.path.relpath(os.path.join(s["run"], "index.html"), out)))
    w.append("</tbody></table>")
    w.append("<h2>Sex-chromosome aneuploidies</h2>")
    if sex_ane:
        by_label = {}
        for r in sex_ane:
            by_label.setdefault(r["label"], []).append(r)
        parts = []
        for label, rs in sorted(by_label.items(), key=lambda kv: -len(kv[1])):
            o = {k: sum(1 for r in rs if k in (r["origin"] or "")) for k in ("maternal", "paternal")}
            stg = {k: sum(1 for r in rs if (r["stage"] or "").startswith(k)) for k in ("meiosis I ", "meiosis I", "meiosis II", "mitotic")}
            detail = []
            if o["maternal"] or o["paternal"]:
                detail.append("%d maternal, %d paternal" % (o["maternal"], o["paternal"]))
            mi = sum(1 for r in rs if (r["stage"] or "").startswith("meiosis I") and not (r["stage"] or "").startswith("meiosis II"))
            if mi or stg["meiosis II"] or stg["mitotic"]:
                detail.append("meiosis I %d, meiosis II %d, mitotic or meiosis II without a crossover %d" % (mi, stg["meiosis II"], stg["mitotic"]))
            parts.append("%s: %d%s" % (label, len(rs), " (%s)" % "; ".join(detail) if detail else ""))
        w.append("<p>%d whole-X or whole-Y events in %d members (sex_aneuploidies.tsv). %s.</p>" % (len(sex_ane), len({(r["trio"], r["role"]) for r in sex_ane}), html.escape("; ".join(parts))))
        w.append(_table(["trio", "sample", "role", "complement", "karyotype", "label", "type", "f", "origin", "stage", "crossovers", "mie_rate", "note"], sex_ane))
    else:
        w.append("<p>None: every whole X and Y matches the complement the pedigree sex implies (sex_aneuploidies.tsv is empty).</p>")
    if conc:
        w.append("<h2>Concordance with the supplied events</h2><p>%d supplied events, %d matched by a triokaryo event of the same sample and chromosome whose "
                 "intersection covers at least half of the shorter segment; %d triokaryo gains or losses without a supplied event; %d copy-neutral events, which a "
                 "depth-based method cannot detect.</p>" % (conc["external"], conc["matched"], conc["new"], conc["cn"]))
        w.append(_table(["sample", "chrom", "start", "end", "label", "match"], [dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in ext]))
    if rejected:
        by_reg = {}
        for r in rejected:
            k = (r["chrom"], int(float(r["start"]) // 5e6))
            by_reg.setdefault(k, []).append(r)
        w.append("<h2>Segments rejected by the phased scan</h2><p>%d segments in %d trios (rejected.all.tsv): span under 2 Mb, windows inconsistent with the segment "
                 "mean, or a bands-only shift under its floor. A region recurring across trios indicates paralogous sequence or a dense cluster of sites rather "
                 "than an event.</p>" % (len(rejected), len({r["trio"] for r in rejected})))
        w.append("<table><tr><th>region</th><th>segments</th><th>trios</th><th>reasons</th></tr>")
        for (c, b), lst in sorted(by_reg.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            reasons = sorted({r["reason"].split(" (")[0] for r in lst})
            w.append("<tr><td>%s:%d-%d Mb</td><td>%d</td><td>%d</td><td>%s</td></tr>" % (c, b * 5, b * 5 + 5, len(lst), len({r["trio"] for r in lst}), html.escape("; ".join(reasons))))
        w.append("</table>")
    w.append("<h2>Colours</h2>" + key_table(["depth", "step", "baf", "phased", "gain", "loss", "loh", "mat", "pat", "cmat", "cpat", "ctrans", "cuntrans", "ext"]))
    w.append("<script>%s</script></body></html>" % JS)
    with open(os.path.join(out, "index.html"), "w") as fh:
        fh.write("\n".join(w))
    if log:
        log("%d trios, %d events -> %s" % (len(summaries), len(events), out))
    return dict(trios=len(summaries), events=len(events), concordance=conc, by_type=by_type, by_source=by_src)
