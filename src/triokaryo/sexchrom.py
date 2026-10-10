"""Sex chromosomes: each member's X and Y copy number from the depth, the sex-chromosome complement (checked against the
pedigree sex), and whole-chromosome X and Y events, constitutional or mosaic, relative to the normal complement the member's
own Y implies - XY with a Y, XX without - with the pedigree sex standing in only where the Y is unreadable.
Segmental X events are called by the ordinary segmentation relative to the member's own X level (segment.call_member);
the Y is analysed as a whole chromosome only.

Copy number: 2 x 2^(median corrected LRR) over the chromosome's usable bins (outside the pseudoautosomal regions, at
least 20 sites, unmasked), kept raw and rounded. A member with depth at almost no Y site while the VCF carries Y sites
has 0 Y copies. Read against the autosomes alone, a male's Y carries a mappability deficit of roughly 5-10%; a panel
whose Y rows hold the male level on the diploid scale removes it, as does comparing a father's Y with his son's, whose
Y is the same sequence (y_father_son).

Whole-chromosome events: the expected complement is the one the member's Y implies (anchor_complement) - XY with a Y, XX
without, a Y in part of the cells going with the nearer of the two - and the pedigree sex only where the Y is unreadable
(effective_sex). The pedigree sex cannot be the anchor: called from X heterozygosity, as peddy's is, it reads a 47,XXY as
female (two X's) and a 45,X as male (one X), and against "female" an XXY's two X copies agree and its Y goes unexamined.
An X whose raw copy number deviates from the expectation by at least min_f is a whole-X gain or loss in a cell fraction
f = |copies - expected|: 47,XXY (f = 1) or 46,XY/47,XXY with a Y; 45,X, 45,X/46,XX or 47,XXX without. A Y deviating from
one copy by at least min_f is a whole-Y gain (47,XYY) or loss (mosaic loss of Y, the common age-related mosaicism in men),
from the panel-corrected level when the panel carries Y rows, else from the father/son ratio when the child is male, else
only beyond a deviation of 0.25; a Y in part of the cells of an XX genome (Y_PRESENT or more) is a whole-Y gain
(46,XX/47,XXY mosaic, or an XX/XY mixture). These events then receive the parent of origin, the phased statistics and the
meiotic stage like an autosome's, with a hemizygous baseline for the X of a member with a Y (trio.parent_of_origin,
phase.annotate_events). The complement against the pedigree sex (sex_check) stays the tell for a sample swap, an XX male
or an XY female, and names the 47,XXY or 45,X a heterozygosity-called sex mislabels."""
import numpy as np

from .model import NA, MEMBERS, chrom_level
from .segment import Event

Y_MIN_BINS = 3                 # usable Y bins for a Y copy number
X_OFFSET_MAX = 0.3             # the largest within-trio X correction applied (log2); beyond it the members' X levels contradict the pedigree
Y_EVENT_FLOOR_NO_REF = 0.25    # a whole-Y deviation this large is reported without a panel or a father/son comparison
Y_PRESENT = 0.25               # a Y read at this many copies or more is a Y, in all or part of the cells; less is stray depth
Y_FATHER_SON_MIN_SITES = 200   # Y sites with depth in both for the father/son ratio
EXPECTED = {"M": (1, 1), "F": (2, 0)}
KARYOTYPE_NAME = {"XXY": "47,XXY", "XYY": "47,XYY", "XXX": "47,XXX", "X": "45,X", "XXYY": "48,XXYY", "XXXY": "48,XXXY", "XXXX": "48,XXXX"}


def karyotype(x_copies, y_copies):
    """The sex-chromosome complement as a string (XY, XX, XXY, X, XYY, ...); '' when either copy number is unknown."""
    if not (np.isfinite(x_copies) and np.isfinite(y_copies)):
        return ""
    x, y = max(int(x_copies), 0), max(int(y_copies), 0)
    return ("X" * x + "Y" * y) if x + y else "0"


def anchor_complement(x_raw, y_raw, y_copies):
    """From the depth alone, the normal complement the member's sex chromosomes are read against: XY ('M') with a Y - the Y
    rounding to one copy or more - and XX ('F') with none, under Y_PRESENT copies. A Y in part of the cells (Y_PRESENT to
    half) goes with the nearer of the two: XY for an X at one copy (mosaic loss of Y), XX for an X at two (a Y in part of the
    cells of an XX genome). '' when the Y is unreadable. The pedigree sex plays no part: called from X heterozygosity, as
    peddy's is, it reads a 47,XXY as female and a 45,X as male."""
    if not np.isfinite(y_raw):
        return ""
    if np.isfinite(y_copies) and int(y_copies) >= 1:
        return "M"
    if y_raw < Y_PRESENT:
        return "F"
    if not np.isfinite(x_raw):
        return "M"
    return "M" if abs(x_raw - 1.0) + (1.0 - y_raw) <= abs(x_raw - 2.0) + y_raw else "F"


def sex_state(bins, scan, m, min_sites=20, min_bins=5):
    """The member's X and Y copy numbers (raw and rounded) and usable bins, and the complement the Y implies (anchor: M, F,
    or '' when the Y is unreadable). y_copies is the Y count of the base complement: 0 without a Y (the VCF carries Y sites
    but the member has depth at fewer than 5% of them, or the level is under Y_PRESENT), at least 1 with one, so that a Y
    lost in part of the cells still reads XY."""
    x_raw, x_cp, nx = chrom_level(bins, m, "chrX", min_sites, min_bins)
    y_raw, y_cp, ny = chrom_level(bins, m, "chrY", min_sites, Y_MIN_BINS)
    ys = scan.sites("chrY")
    y_in_vcf = ys is not None and ys.n > 0
    if y_in_vcf and not np.isfinite(y_raw) and float((ys.dp[m] > 0).mean()) < 0.05:
        y_raw, y_cp = 0.0, 0
    anchor = anchor_complement(x_raw, y_raw, y_cp)
    if anchor == "M" and y_cp < 1:
        y_cp = 1
    elif anchor == "F":
        y_cp = 0
    sl = bins.of("chrY")
    y_total = int(sl.stop - sl.start)
    y_masked = int(bins.masked[sl].sum()) if bins.masked is not None and y_total else 0
    y_thin = int((bins.n_dp[m][sl] < min_sites).sum()) if y_total else 0
    return dict(x_copies=x_cp, x_copies_raw=x_raw, x_bins=nx, y_copies=y_cp, y_copies_raw=y_raw, y_bins=ny, y_in_vcf=y_in_vcf, y_bins_total=y_total,
                y_bins_masked=y_masked, y_bins_thin=y_thin, y_min_sites=min_sites, anchor=anchor)


def y_unread(st):
    """Why the Y copy number is not read: the Y is not in the VCF, or its bins are unusable - masked by the panel, or holding fewer
    than min_sites depth sites each - and fewer than Y_MIN_BINS remain."""
    if not st.get("y_in_vcf"):
        return "Y not in the VCF"
    return "Y unreadable: %d usable Y bin(s) of %d, %d needed (%d masked by the panel, %d with fewer than %d depth sites)" % (
        st.get("y_bins", 0), st.get("y_bins_total", 0), Y_MIN_BINS, st.get("y_bins_masked", 0), st.get("y_bins_thin", 0), st.get("y_min_sites", 20))


def x_check(sex, st):
    """The X copy number against the pedigree sex, as summary.tsv has always reported it (Y-aware where the Y is known)."""
    xc = st["x_copies"]
    if not sex or not np.isfinite(xc):
        return ""
    if (sex == "M" and xc == 1) or (sex == "F" and xc == 2):
        return "agrees"
    k = karyotype(xc, st["y_copies"])
    if k:
        why = "an XX male, or a sample swap" if (sex == "M" and k == "XX") else "a sample swap, or an XY female" if (sex == "F" and k == "XY") else KARYOTYPE_NAME.get(k, k)
    else:
        why = {("M", 2): "47,XXY or an XX male?", ("F", 1): "45,X?", ("F", 3): "47,XXX?", ("M", 3): "48,XXXY?"}.get((sex, int(xc)), y_unread(st))
    return "X copies %d in a reported %s (%s)" % (int(xc), "male" if sex == "M" else "female", why)


def implied_sex(sex, st):
    """The pedigree sex, or, when it is not given, the sex the Y copy number implies (M with a Y, F without); '' when neither."""
    if sex:
        return sex, False
    yc = st.get("y_copies", NA)
    if np.isfinite(yc):
        return ("M" if int(yc) >= 1 else "F"), True
    return "", False


def effective_sex(sex, st):
    """The sex the member's whole X and Y are read against, and its source: the complement the Y implies where the Y is read
    (anchor: 'Y'), else the pedigree sex ('pedigree'), else ''."""
    a = st.get("anchor", "")
    if a:
        return a, "Y"
    return (sex, "pedigree") if sex else ("", "")


def _who(sex, src, ped):
    """The member in the event's note: how its sex was fixed for the expected complement (no parentheses: the cohort table takes
    the note's first parenthesis as the event's name)."""
    name = "male" if sex == "M" else "female"
    if src != "Y":
        return "a reported %s, the Y unread" % name
    if not ped:
        return "a %s with no pedigree sex, implied by the Y" % name
    if ped == sex:
        return "a reported %s" % name
    return "a reported %s %s" % ("male" if ped == "M" else "female", "with a Y" if sex == "M" else "without a Y")


def sex_check(sex, st):
    """The sex-chromosome complement against the pedigree sex: 'agrees', or the complement found and its name."""
    if not sex:
        return ""
    k = karyotype(st["x_copies"], st["y_copies"])
    expected = "XY" if sex == "M" else "XX"
    if not k:
        if not np.isfinite(st["x_copies"]):
            return ""
        why = y_unread(st)
        return ("agrees (%s)" % why) if int(st["x_copies"]) == EXPECTED[sex][0] else "X copies %d in a reported %s (%s)" % (
            int(st["x_copies"]), "male" if sex == "M" else "female", why)
    if k == expected:
        return "agrees"
    if sex == "M" and k == "XX":
        name = "an XX male, a sample swap, or an X without heterozygosity (an isodisomy) that a sex called from X heterozygosity reads male"
    elif sex == "F" and k == "XY":
        name = "a sample swap, or an XY female"
    elif sex == "F" and k == "XXY":
        name = "47,XXY; a sex called from X heterozygosity reads two X's as female"
    elif sex == "M" and k == "X":
        name = "45,X; a sex called from X heterozygosity reads one X as male"
    else:
        name = KARYOTYPE_NAME.get(k, k)
    return "%s in a reported %s (%s)" % (k, "male" if sex == "M" else "female", name)


def x_offset_within_trio(bins, scan, sexes, min_sites=20, min_bins=5):
    """Without a panel, the X read against the autosomes carries a mappability deficit of about 5-10%, which would read as a
    mosaic X loss in every female. The correction is the median, over the members with a pedigree or Y-implied sex, of the
    deviation of the X level (log2 of raw copies over 2) from its expectation (log2(1/2) for a male, 0 for a female): robust to
    one aneuploid member among three, and not applied beyond X_OFFSET_MAX. Returns (offset, {role: deviation})."""
    devs = {}
    for m, role in enumerate(MEMBERS):
        st = sex_state(bins, scan, m, min_sites, min_bins)
        sex, _ = implied_sex(sexes[m], st)
        if not sex or not np.isfinite(st["x_copies_raw"]) or st["x_copies_raw"] <= 0:
            continue
        devs[role] = float(np.log2(st["x_copies_raw"] / 2.0) - np.log2(EXPECTED[sex][0] / 2.0))
    if not devs:
        return 0.0, devs
    off = float(np.median(list(devs.values())))
    return (off if abs(off) <= X_OFFSET_MAX else 0.0), devs


def apply_x_offset(bins, offset):
    """Subtract the within-trio X offset from every member's corrected X LRR, in place."""
    sl = bins.of("chrX")
    if sl.start < sl.stop and offset:
        for m in range(bins.lrr_gc.shape[0]):
            v = bins.lrr_gc[m, sl]
            bins.lrr_gc[m, sl] = np.where(np.isfinite(v), v - offset, v)


def y_father_son(scan, bins, min_sites=Y_FATHER_SON_MIN_SITES):
    """log2 of the father's Y depth over the son's, each relative to its autosomal median depth, as the median over the Y sites
    with depth in both (outside the pseudoautosomal regions), and the site count. The two Y chromosomes are the same sequence,
    so mapping cancels: a negative value is a lower Y copy number in the father (mosaic loss of Y), a positive one in the son."""
    s = scan.sites("chrY")
    if s is None:
        return NA, 0
    ok = (s.dp[0] > 0) & (s.dp[1] > 0) & ~s.par
    if ok.sum() < min_sites:
        return NA, int(ok.sum())
    med = [float(np.nanmedian(bins.depth[i][bins.autosomal])) for i in (0, 1)]
    if not all(np.isfinite(med)) or min(med) <= 0:
        return NA, int(ok.sum())
    r = np.log2((s.dp[1][ok] / med[1]) / (s.dp[0][ok] / med[0]))
    return float(np.median(r)), int(ok.sum())


def sex_chromosome_events(m, sample, states, sexes, bins, genome, min_f, y_ref=False, y_ratio=NA):
    """The member's whole-X and whole-Y events against the complement its own Y implies (effective_sex: XY with a Y, XX without;
    the pedigree sex where the Y is unreadable). states: {role: sex_state}; y_ref: the panel carries Y rows, so the Y level is
    corrected; y_ratio: y_father_son's log2 ratio (NaN without a son)."""
    role = MEMBERS[m]
    st = states[role]
    out = []
    sex, src = effective_sex(sexes[m], st)
    if not sex:
        return out
    ex_x, ex_y = EXPECTED[sex]
    who = _who(sex, src, sexes[m])
    k = karyotype(st["x_copies"], st["y_copies"])
    x_raw = st["x_copies_raw"]
    if np.isfinite(x_raw) and abs(x_raw - ex_x) >= min_f:
        kind = "gain" if x_raw > ex_x else "loss"
        f = float(min(abs(x_raw - ex_x), 2.0 if kind == "gain" else 1.0))
        L = genome.length["chrX"]
        name = ""
        if f >= 0.9:
            name = (" (%s)" % KARYOTYPE_NAME[k]) if k in KARYOTYPE_NAME else {("M", "gain"): " (47,XXY)", ("F", "loss"): " (45,X)", ("F", "gain"): " (47,XXX)"}.get((sex, kind), "")
        elif kind == "gain":
            name = " (46,XY/47,XXY mosaic)" if sex == "M" else " (46,XX/47,XXX mosaic)"
        elif kind == "loss":
            name = " (45,X/46,%s mosaic)" % ("XY" if sex == "M" else "XX")
        note = "whole-chromosome X: %.2f copies against %d expected for %s%s" % (x_raw, ex_x, who, name)
        if sex == "M" and kind == "loss":
            note += "; the single X is maternal"
        out.append(Event(sample, role, "chrX", 0, L, "whole", kind, lrr=float(np.log2(x_raw / 2.0)) if x_raw > 0 else NA, n_bins=st["x_bins"], f_lrr=f, note=note))
    y_raw = st["y_copies_raw"]
    if not (np.isfinite(y_raw) and st["y_in_vcf"]):
        return out
    if sex == "M":
        dev = y_raw - ex_y
        kind = "gain" if dev > 0 else "loss"
        f_level = float(min(abs(dev), 2.0 if kind == "gain" else 1.0))
        f_ratio = NA
        if np.isfinite(y_ratio):
            other = states["child" if role == "father" else "father"]
            lower = (role == "father" and y_ratio < 0) or (role == "child" and y_ratio > 0)
            if lower and np.isfinite(other["y_copies_raw"]) and y_raw <= other["y_copies_raw"]:
                f_ratio = float(1.0 - 2.0 ** (-abs(y_ratio)))
        if y_ref and f_level >= min_f:
            f, how = f_level, "from the Y depth against the panel's male level"
            if np.isfinite(f_ratio):
                how += "; the father/son Y ratio gives %.0f%%" % (100 * f_ratio)
        elif np.isfinite(f_ratio) and f_ratio >= min_f and kind == "loss":
            f, how = f_ratio, "from the father/son Y ratio (the %s's Y taken as one copy; the Y depth alone reads %.2f copies, uncorrected for mappability)" % (
                "son" if role == "father" else "father", y_raw)
        elif not y_ref and not np.isfinite(y_ratio) and f_level >= Y_EVENT_FLOOR_NO_REF:
            f, how = f_level, "from the Y depth alone, uncorrected for mappability (no panel with Y rows, no son to compare)"
        else:
            f = NA
        if np.isfinite(f):
            name = (" (%s)" % KARYOTYPE_NAME.get(k, "47,XYY")) if kind == "gain" and f >= 0.9 else " (mosaic loss of Y)" if kind == "loss" and f < 0.9 else " (loss of Y)" if kind == "loss" else ""
            out.append(Event(sample, role, "chrY", 0, genome.length["chrY"], "whole", kind, lrr=float(np.log2(y_raw / 2.0)) if y_raw > 0 else NA, n_bins=st["y_bins"],
                             f_lrr=f, note="whole-chromosome Y%s: %.2f copies against 1 expected for %s, %s" % (name, y_raw, who, how)))
    elif y_raw >= Y_PRESENT:                                   # a Y in part of the cells of an XX genome
        f = float(min(y_raw, 1.0))
        if y_ref:
            how = "from the Y depth against the panel's male level"
        elif f >= Y_EVENT_FLOOR_NO_REF:
            how = "from the Y depth alone, uncorrected for mappability (no panel with Y rows)"
        else:
            return out
        out.append(Event(sample, role, "chrY", 0, genome.length["chrY"], "whole", "gain", lrr=float(np.log2(y_raw / 2.0)), n_bins=st["y_bins"], f_lrr=f,
                         note="whole-chromosome Y (46,XX/47,XXY mosaic, or an XX/XY mixture): a Y in %.0f%% of the cells of %s, %.2f copies against none expected, %s" % (
                             100 * f, who, y_raw, how)))
    return out


STAGE_SHORT = (("meiosis II", "MII"), ("meiosis I", "MI"), ("mitotic", "MII/mit"))


def _short_stage(stage):
    for k, v in STAGE_SHORT:
        if stage.startswith(k):
            return v
    return ""


def _short_origin(e):
    """mat or pat: the parent of the extra, lost or retained copy, from the child's origin labels; '' for a parent's event."""
    if e.role != "child":
        return ""
    s = e.origin_phase or e.origin
    i, j = s.find("maternal"), s.find("paternal")
    if i < 0 and j < 0:
        return ""
    return "mat" if (j < 0 or (0 <= i < j)) else "pat"


def event_bands(e, genome):
    """The cytogenetic bands an event spans, from its site-resolution boundaries where available ('' without a band table)."""
    a = e.start_fine if np.isfinite(e.start_fine) else e.start
    b = e.end_fine if np.isfinite(e.end_fine) else e.end
    return genome.band_range(e.chrom, int(a), int(b)) if genome is not None else ""


def _coords(e, genome=None):
    """A segment's extent for the karyotype string: its bands (ISCN) when the genome carries a band table, else Mb."""
    bands = e.bands or event_bands(e, genome)
    if bands:
        return bands
    a = e.start_fine if np.isfinite(e.start_fine) else e.start
    b = e.end_fine if np.isfinite(e.end_fine) else e.end
    return "%s%.1f-%.1fMb" % ((e.span + ":") if e.span in ("p", "q") else "", a / 1e6, b / 1e6)


def karyotype_string(events, complement, chrom_order, genome=None):
    """An ISCN-like summary of one member. The main line is the modal number and sex-chromosome complement, annotated for a
    constitutional sex-chromosome aneuploidy with the parent of the extra or lost copy and the meiotic stage (47,XXY(mat,MI)),
    followed by the constitutional autosomal terms: +N / -N with mat/pat and the stage, upd(N)mat(iso|hetero), roh(N)(bands),
    loh/dup/del(N)(bands) with mat/pat; the bands are cytogenetic (q22.1q31.1) with a band table, else Mb. Each mosaic event
    forms its own line relative to the base complement, with its cell fraction in brackets: mos 47,XXY(pat)[0.40]/46,XY. A
    term ending in ? is a call the phased track doubts."""
    comp = complement or "?"
    events = sorted(events, key=lambda e: (chrom_order.get(e.chrom, 99), e.start))
    is_mos = lambda e: np.isfinite(e.f) and e.f < 0.9  # noqa: E731
    autosomal_whole = lambda e: e.chrom not in ("chrX", "chrY") and e.span == "whole" and e.type in ("gain", "loss")  # noqa: E731
    # the modal number: the complement, plus the constitutional whole-chromosome gains and losses (every line carries them)
    modal = (44 + len(comp) if complement else 46) + sum((1 if e.type == "gain" else -1) for e in events if autosomal_whole(e) and not is_mos(e))
    const, mosaic, ann = [], [], ""
    for e in events:
        mos = is_mos(e)
        doubt = "?" if "may be an artefact" in e.note else ""
        n = e.chrom[3:]
        o, stg = _short_origin(e), _short_stage(e.stage)
        if e.chrom in ("chrX", "chrY") and e.span == "whole" and e.type in ("gain", "loss"):   # an X LOH or heterodisomy is a upd/loh term below
            a = ",".join(x for x in (o, stg) if x)
            if not mos:
                ann = "(%s)" % a if a else ""
                continue
            comp2, m2 = comp, modal
            if complement:
                comp2 = (comp + n) if e.type == "gain" else comp.replace(n, "", 1)
                comp2 = "X" * comp2.count("X") + "Y" * comp2.count("Y")
                m2 = modal + len(comp2) - len(comp)
            mosaic.append(("%d,%s%s%s" % (m2, comp2, "(%s)" % a if a else "", doubt), e.f))
            continue
        if e.type == "UPD" or (e.type == "LOH" and e.span == "whole" and "isodisomy" in e.note):
            term = "upd(%s)%s(%s)" % (n, o, "hetero" if e.type == "UPD" else "iso")
        elif e.type == "LOH" and "run of homozygosity" in e.note:
            term = "roh(%s)(%s)" % (n, _coords(e, genome))
        elif e.type == "LOH":
            term = "loh(%s)(%s)%s" % (n, _coords(e, genome), o)
        elif e.span == "whole":
            term = "%s%s%s%s" % ("+" if e.type == "gain" else "-", n, o, "(%s)" % stg if stg else "")
        else:
            term = "%s(%s)(%s)%s" % ("dup" if e.type == "gain" else "del", n, _coords(e, genome), o)
        term += doubt
        if mos:
            m2 = modal + (1 if (autosomal_whole(e) and e.type == "gain") else -1 if autosomal_whole(e) else 0)
            mosaic.append(("%d,%s,%s" % (m2, comp, term), e.f))
        else:
            const.append(term)
    main = "%d,%s%s" % (modal, comp, ann) + ("," + ",".join(const) if const else "")
    if mosaic:
        return "mos " + "/".join("%s[%.2f]" % (line, f) for line, f in mosaic) + "/" + main
    return main
