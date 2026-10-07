"""One trio, end to end: scan the VCF, bin, call each member, read the trio, write tables, figures and the page."""
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
        log("panel: %d bins left out of the calls (unpinned: fewer than %d genomes, or spread beyond %.2f)" % (int(bins.masked.sum()), PANEL_MIN_N, PANEL_MAX_RSD))
    events, x_copies = [], {}
    for m, role in enumerate(MEMBERS):
        ev, xc = call_member(bins, scan, m, trio.members[m], G, P, trio.sexes[m])
        events += ev
        x_copies[role] = xc
    base_mie = read_trio(events, scan, G, min_dp, min_gq)
    external = read_events(events_path, G) if events_path else []
    external = [x for x in external if x.sample in trio.members]
    match_external(events, external)
    os.makedirs(out, exist_ok=True)
    summ = write_tables(out, trio, bins, events, x_copies, scan, base_mie, P)
    if external:
        from .report import write_tsv
        write_tsv(os.path.join(out, "external.tsv"), ["sample", "chrom", "start", "end", "label", "match"],
                  [dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in external])
    figs = []
    if figures:
        from .plots import fig_chrom, fig_genome
        fdir = os.path.join(out, "figures")
        os.makedirs(fdir, exist_ok=True)
        figs.append(fig_genome(trio, bins, scan, events, external, G, fdir, min_dp=min_dp, min_gq=min_gq))
        chroms = sorted({e.chrom for e in events} | {x.chrom for x in external}, key=lambda c: G.chroms.index(c))
        for c in chroms:
            figs.append(fig_chrom(trio, c, bins, scan, events, external, G, fdir, min_dp, min_gq))
    write_html(out, trio, figs, events, summ, external, mock_note)
    with open(os.path.join(out, "summary.json"), "w") as fh:
        json.dump(dict(trio=trio.name, members=list(trio.members), sexes=list(trio.sexes), x_copies=x_copies, mie_rate_genome=base_mie,
                       events=[e.as_dict() for e in events], external=[dict(sample=x.sample, chrom=x.chrom, start=x.start, end=x.end, label=x.note, match=x.inheritance) for x in external],
                       sites_used=scan.n_used, records=scan.n_records, skipped=scan.skipped, seconds=round(time.time() - t0, 1)),
                  fh, indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else str(o))
    log("%d event(s): %s; %.0f s -> %s" % (len(events), "; ".join("%s %s %s %s f %.2f" % (e.role, e.type, e.span, e.chrom, e.f) for e in events) or "none", time.time() - t0, out))
    return dict(events=events, x_copies=x_copies, bins=bins, scan=scan, summary=summ, figures=figs, external=external)
