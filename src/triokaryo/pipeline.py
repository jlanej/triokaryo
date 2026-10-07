"""One trio end to end: scan the VCF, bin, segment and call each member, phase, run the trio analyses, write tables,
figures and the page."""
import json
import os
import sys
import time

import numpy as np

from .external import match_external, read_events
from .genome import genome as load_genome
from .model import MEMBERS, load_gc_track, make_bins
from .report import write_html, write_tables
from .segment import PARAMS, call_member
from .trio import read_trio
from .vcfscan import scan_vcf


def run_trio(vcf, trio, out, gc_track=None, events_path=None, bin_size=1_000_000, min_dp=8, min_gq=20, thin=1, genome_name="grch38", figures=True,
             params=None, log=None, mock_note="", panel=None):
    log = log or (lambda s: print("[triokaryo] " + s, file=sys.stderr))
    t0 = time.time()
    G = load_genome(genome_name)
    P = dict(PARAMS, **(params or {}), min_dp=min_dp, min_gq=min_gq)
    log("trio %s: child %s, father %s, mother %s" % (trio.name, trio.kid, trio.dad, trio.mom))
    scan = scan_vcf(vcf, trio.members, G, thin=thin, log=log)
    log("%d records, %d PASS biallelic SNVs used (%s); %.0f s" % (scan.n_records, scan.n_used, ", ".join("%s %d" % kv for kv in sorted(scan.skipped.items())), time.time() - t0))
    gc = load_gc_track(gc_track) if gc_track else None
    pan = None
    if panel:
        from .panel import load_panel
        pan = load_panel(panel)
        log("panel: %d bins from %s" % (len(pan), panel))
    bins = make_bins(scan, G, bin_size, min_dp, min_gq, gc, panel=pan)
    if pan is not None:
        from .model import PANEL_MAX_RSD, PANEL_MIN_N
        log("panel: %d bins masked (fewer than %d genomes, or robust SD above %.2f)" % (int(bins.masked.sum()), PANEL_MIN_N, PANEL_MAX_RSD))
    events, x_copies = [], {}
    for m, role in enumerate(MEMBERS):
        ev, xc = call_member(bins, scan, m, trio.members[m], G, P, trio.sexes[m])
        events += ev
        x_copies[role] = xc
    # the sex chromosomes: each member's X and Y copy number, the complement against the pedigree sex, whole-X and whole-Y events
    from .sexchrom import sex_check, sex_chromosome_events, sex_state, x_check, y_father_son
    states = {role: sex_state(bins, scan, m, P["min_sites"], P["min_len"]) for m, role in enumerate(MEMBERS)}
    y_ref = bool(pan is not None and np.isfinite(bins.panel_median[bins.of("chrY")]).sum() >= 3)
    y_ratio, y_ratio_n = y_father_son(scan, bins) if states["child"]["y_copies"] == 1 and states["father"]["y_copies"] == 1 else (float("nan"), 0)
    for m, role in enumerate(MEMBERS):
        states[role]["x_check"] = x_check(trio.sexes[m], states[role])
        states[role]["sex_check"] = sex_check(trio.sexes[m], states[role])
        events += sex_chromosome_events(m, trio.members[m], states, trio.sexes, bins, G, P["min_f"], y_ref=y_ref, y_ratio=y_ratio)
    sex_info = dict(states=states, y_father_son_log2=y_ratio, y_father_son_sites=y_ratio_n, y_panel=y_ref)
    log("sex chromosomes: %s; Y father/son log2 ratio %s" % ("; ".join("%s %s (%s)" % (role, states[role]["sex_check"] and (("X" * int(states[role]["x_copies"]) if np.isfinite(states[role]["x_copies"]) else "?") +
                                                                      ("Y" * int(states[role]["y_copies"]) if np.isfinite(states[role]["y_copies"]) else "")) or "unknown", states[role]["sex_check"] or "no pedigree sex")
                                                                      for role in MEMBERS), "%.3f over %d sites" % (y_ratio, y_ratio_n) if np.isfinite(y_ratio) else "NA"))
    child_x_baseline = 1 if trio.kid_sex == "M" else 2 if trio.kid_sex == "F" else (int(x_copies["child"]) if np.isfinite(x_copies["child"]) else 2)
    # the phased tracks: the maternal-allele fraction along the child, the transmitted-allele fraction along each parent
    from .phase import annotate_events, mask_rejected, mask_shared, phased_scan, phased_tracks
    tracks, phase_info, rejected = [], {}, []
    for m, role in enumerate(MEMBERS):
        t, info = phased_tracks(scan, bins, G, m, min_dp, min_gq, trio.sexes[m], child_x_copies=x_copies["child"])
        tracks.append(t)
        phase_info[role] = info
        log("%s: %d phased sites (%d more on the auxiliary tracks), windows of %d, reference bias %+.3f" % (role, info["n_phased"], info["n_aux"], info["window_sites"], info["ref_bias"]))
    mask_shared(tracks, bins, x_copies)
    for m, role in enumerate(MEMBERS):
        phase_info[role]["windows_shared"] = int(sum(int(t.w_shared.sum()) for t in tracks[m].values() if t.w_shared is not None))
        auto = bins.autosomal & (bins.n_called[m] >= 10)
        base_het = float(np.nanmedian(bins.het_rate[m][auto])) if auto.any() else float("nan")
        found = phased_scan(tracks[m], bins, scan, m, trio.members[m], G, events, None, min_dp, min_gq, trio.sexes[m], base_het, rejected, x_copies=x_copies[role])
        if found:
            log("%s: %d event(s) from the phased scan not called by the depth" % (role, len(found)))
        events += found
    mask_rejected(tracks, bins, rejected)
    if rejected:
        log("%d phased segment(s) rejected: %s" % (len(rejected), "; ".join("%s %s %.1f-%.1f Mb (%s)" % (r["role"], r["chrom"], r["start"] / 1e6, r["end"] / 1e6, r["reason"].split(" (")[0]) for r in rejected)))
    events.sort(key=lambda e: (MEMBERS.index(e.role), G.chroms.index(e.chrom), e.start))
    base_mie = read_trio(events, scan, G, min_dp, min_gq, child_x_copies=x_copies["child"], child_x_baseline=child_x_baseline)
    annotate_events(events, tracks, scan, bins, G, child_x_baseline=child_x_baseline)
    external = read_events(events_path, G) if events_path else []
    external = [x for x in external if x.sample in trio.members]
    match_external(events, external)
    os.makedirs(out, exist_ok=True)
    summ = write_tables(out, trio, bins, events, x_copies, scan, base_mie, P, tracks=tracks, phase_info=phase_info, rejected=rejected, sex_info=sex_info)
    if external:
        from .report import write_tsv
        write_tsv(os.path.join(out, "external.tsv"), ["sample", "chrom", "start", "end", "label", "match"],
                  [dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in external])
    figs = []
    if figures:
        from .plots import fig_chrom, fig_genome
        fdir = os.path.join(out, "figures")
        os.makedirs(fdir, exist_ok=True)
        figs.append(fig_genome(trio, bins, scan, events, external, G, fdir, min_dp=min_dp, min_gq=min_gq, tracks=tracks))
        chroms = sorted({e.chrom for e in events} | {x.chrom for x in external}, key=lambda c: G.chroms.index(c))
        for c in chroms:
            figs.append(fig_chrom(trio, c, bins, scan, events, external, G, fdir, min_dp, min_gq, tracks=tracks))
    write_html(out, trio, figs, events, summ, external, mock_note)
    from .guide import write_guide
    write_guide(os.path.join(out, "guide.html"))
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(dict(trio=trio.name, members=list(trio.members), sexes=list(trio.sexes), x_copies=x_copies, mie_rate_genome=base_mie, phasing=phase_info,
                       sex_chromosomes={role: {k: v for k, v in states[role].items() if k != "y_in_vcf"} for role in MEMBERS},
                       y_father_son=dict(log2=y_ratio, sites=y_ratio_n), genome=genome_name, mock_note=mock_note,
                       events=[e.as_dict() for e in events], external=[dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in external],
                       sites_used=scan.n_used, records=scan.n_records, skipped=scan.skipped, seconds=round(time.time() - t0, 1)),
                  fh, indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else str(o))
    log("%d event(s): %s; %.0f s -> %s" % (len(events), "; ".join("%s %s %s %s f %.2f" % (e.role, e.type, e.span, e.chrom, e.f) for e in events) or "none", time.time() - t0, out))
    return dict(events=events, x_copies=x_copies, bins=bins, scan=scan, summary=summ, figures=figs, external=external, tracks=tracks, phase_info=phase_info, sex=sex_info)
