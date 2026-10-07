"""Trios: from a trios file (`#kid dad mom [kid_sex dad_sex mom_sex]`, the header naming the columns), a 6-column PED
(family id father mother sex phenotype; a child is a row with both parents), or the three IDs given outright."""
from dataclasses import dataclass

SEX = {"1": "M", "2": "F", "m": "M", "f": "F", "male": "M", "female": "F", "xy": "M", "xx": "F", "0": "", "-9": "", "": "", "na": "", "unknown": ""}
KID = ("kid", "child", "proband", "sample", "sample_id", "iid", "id")
DAD = ("dad", "father", "paternal_id", "fid_father")
MOM = ("mom", "mother", "maternal_id")


@dataclass
class Trio:
    kid: str
    dad: str
    mom: str
    kid_sex: str = ""
    dad_sex: str = ""
    mom_sex: str = ""
    family: str = ""

    @property
    def members(self):
        return (self.kid, self.dad, self.mom)

    @property
    def sexes(self):
        return (self.kid_sex, self.dad_sex, self.mom_sex)

    @property
    def name(self):
        return self.kid


def sex_of(v):
    return SEX.get(str(v or "").strip().lower(), str(v or "").strip().upper()[:1] if str(v or "").strip().upper()[:1] in ("M", "F") else "")


def read_trios(path):
    """Every trio in the file. A header line (starting with '#' or made of words) names the columns; without one, a PED."""
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    if not rows:
        return []
    hdr = [h.lstrip("#").strip().lower() for h in rows[0]]
    out = []
    if any(h in KID for h in hdr) and any(h in DAD for h in hdr) and any(h in MOM for h in hdr):
        ik = next(i for i, h in enumerate(hdr) if h in KID)
        idd = next(i for i, h in enumerate(hdr) if h in DAD)
        im = next(i for i, h in enumerate(hdr) if h in MOM)
        col = lambda names: next((i for i, h in enumerate(hdr) if h in names), None)
        iks, ids, ims = col(("kid_sex", "child_sex", "proband_sex", "sex")), col(("dad_sex", "father_sex")), col(("mom_sex", "mother_sex"))
        ifam = col(("family", "family_id", "fid", "fam"))
        for r in rows[1:]:
            if len(r) <= max(ik, idd, im) or r[idd] in ("", "0") or r[im] in ("", "0"):
                continue
            out.append(Trio(r[ik], r[idd], r[im], sex_of(r[iks]) if iks is not None and len(r) > iks else "",
                            sex_of(r[ids]) if ids is not None and len(r) > ids else "", sex_of(r[ims]) if ims is not None and len(r) > ims else "",
                            r[ifam] if ifam is not None and len(r) > ifam else ""))
        return out
    # a PED: family, id, father, mother, sex, phenotype
    recs = {}
    for r in rows:
        if r[0].startswith("#") or len(r) < 4:
            continue
        recs[r[1]] = r
    for r in recs.values():
        fam, sid, fa, mo = r[0], r[1], r[2], r[3]
        if fa in ("", "0") or mo in ("", "0"):
            continue
        sx = lambda s: sex_of(recs[s][4]) if s in recs and len(recs[s]) > 4 else ""
        out.append(Trio(sid, fa, mo, sx(sid), sx(fa), sx(mo), fam))
    return out


def trio_from_args(child, father, mother, sexes=""):
    """--child/--father/--mother, with --sex 'M,M,F' (child, father, mother) optional."""
    sx = [sex_of(x) for x in (sexes.split(",") if sexes else [])] + ["", "", ""]
    return Trio(child, father, mother, sx[0], sx[1], sx[2])
