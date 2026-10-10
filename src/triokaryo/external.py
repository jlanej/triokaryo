"""Events supplied from another method (a CNV caller, a depth-based karyotype such as NGS-DOSE's, a clinical karyotype),
drawn beside the trio's events and matched to them. A TSV with sample, chrom, start, end and a label (optionally a
type); or NGS-DOSE's karyotype/events.tsv (sample, chrom, span, start_mb, end_mb, label, kind), recognised by its
columns. NGS-DOSE lists a sex-chromosome complement other than XX or XY (47,XXY, 45,X, 47,XYY ...) with chrom chrX/chrY:
such a row is placed on the chromosome the complement changes (complement_chrom)."""
import csv

from .segment import Event


def complement_chrom(label):
    """For a sex-chromosome complement written as a karyotype (47,XXY; 45,X; 47,XYY; 48,XXYY ...): the chromosome it changes and
    the direction - the X when its count is not the one the Y implies (one with a Y, two without), else the Y; (None, '') for a
    label that is no complement, or a normal one."""
    comp = label.strip().split(",")[-1].strip().upper()
    if not comp or set(comp) - {"X", "Y"}:
        return None, ""
    nx, ny = comp.count("X"), comp.count("Y")
    ex_x, ex_y = (1, 1) if ny else (2, 0)
    if nx != ex_x:
        return "chrX", "gain" if nx > ex_x else "loss"
    if ny != ex_y:
        return "chrY", "gain" if ny > ex_y else "loss"
    return None, ""


def read_events(path, genome):
    out = []
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = [c.lower() for c in rd.fieldnames or []]
        for r in rd:
            r = {k.lower(): v for k, v in r.items()}
            sample = r.get("sample") or r.get("id") or ""
            raw_chrom = r.get("chrom") or r.get("chr") or ""
            chrom = genome.normalize(raw_chrom)
            comp_type = ""
            if chrom is None and "/" in raw_chrom:      # a sex-chromosome complement row (NGS-DOSE: chrom chrX/chrY, label 47,XXY ...)
                chrom, comp_type = complement_chrom(r.get("label") or "")
            if not sample or chrom is None:
                continue
            if "start_mb" in cols:                       # NGS-DOSE's karyotype events
                span = "whole" if comp_type else (r.get("span") or "whole").lower()
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
                etype = comp_type or ("gain" if ("gain" in kind or (label.startswith("+"))) else "loss" if ("loss" in kind or label.startswith("-")) else "")
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
    """Each of the trio's events receives the labels of the supplied events of the same sample and chromosome whose
    intersection covers at least min_overlap of the shorter segment; each supplied event is marked matched or unmatched.
    Returns the supplied list with .inheritance used as the match flag."""
    for e in events:
        hits = [x for x in external if x.sample == e.sample and x.chrom == e.chrom and e.overlap(x) >= min_overlap]
        e.external = "; ".join(x.note for x in hits) if hits else ""
    for x in external:
        hits = [e for e in events if e.sample == x.sample and e.chrom == x.chrom and x.overlap(e) >= min_overlap]
        x.inheritance = ("matched: " + "; ".join("%s %s %s f %.2f" % (e.type, e.span, e.chrom, e.f) for e in hits)) if hits else "unmatched"
    return external
