"""Events from elsewhere (another caller, a depth tool such as NGS-DOSE's karyotype, a clinical karyotype) to draw beside
the trio's and to check against. A TSV with sample, chrom, start, end and a label; or NGS-DOSE's karyotype/events.tsv
(sample, chrom, span, start_mb, end_mb, label, kind), read by its columns."""
import csv

from .segment import Event


def read_events(path, genome):
    out = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = [c.lower() for c in rd.fieldnames or []]
        for r in rd:
            r = {k.lower(): v for k, v in r.items()}
            sample = r.get("sample") or r.get("id") or ""
            chrom = genome.normalize(r.get("chrom") or r.get("chr") or "")
            if not sample or chrom is None:
                continue
            if "start_mb" in cols:                       # NGS-DOSE's karyotype events
                span = (r.get("span") or "whole").lower()
                if span == "whole":
                    start, end = 0, genome.length[chrom]
                elif span in ("p", "q"):
                    a, b = genome.arm_bounds(chrom, span)
                    start, end = a - 1, b
                else:
                    try:
                        start, end = int(float(r["start_mb"]) * 1e6), int(float(r["end_mb"]) * 1e6)
                    except (ValueError, KeyError):
                        continue
                label = r.get("label") or r.get("kind") or span
                kind = (r.get("kind") or "").lower()
                etype = "gain" if ("gain" in kind or (label.startswith("+"))) else "loss" if ("loss" in kind or label.startswith("-")) else ""
            else:
                try:
                    start, end = int(float(r.get("start", 0))), int(float(r.get("end", 0)))
                except ValueError:
                    continue
                label = r.get("label") or r.get("type") or ""
                etype = (r.get("type") or "").lower()
            out.append(Event(sample, "", chrom, start, end, r.get("span", ""), etype or "", note=label))
    return out


def match_external(events, external, min_overlap=0.5):
    """Each of the trio's events gets the labels of the external events it overlaps (same sample); each external event is
    marked matched / unmatched. Returns the external list with .inheritance used as the match flag."""
    for e in events:
        hits = [x for x in external if x.sample == e.sample and x.chrom == e.chrom and e.overlap(x) >= min_overlap]
        e.external = "; ".join(x.note for x in hits) if hits else ""
    for x in external:
        hits = [e for e in events if e.sample == x.sample and e.chrom == x.chrom and x.overlap(e) >= min_overlap]
        x.inheritance = ("matched: " + "; ".join("%s %s %s f %.2f" % (e.type, e.span, e.chrom, e.f) for e in hits)) if hits else "unmatched"
    return external
