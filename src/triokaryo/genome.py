"""The assembly: chromosome lengths, the end of each short arm (the centromere), the pseudoautosomal regions and the
cytogenetic bands (UCSC cytoBand). GRCh38, chr-prefixed names; positions 1-based. Other names ('1', 'X') are normalised
to these."""
import os
from dataclasses import dataclass

GRCH38_LENGTH = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
                 "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
                 "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
                 "chr19": 58617616, "chr20": 64444167, "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415}
# the end of the short arm: the boundary between the p and q centromeric (acen) bands of UCSC's hg38 cytoBand, i.e. the start of q11
GRCH38_P_END = {"chr1": 123400000, "chr2": 93900000, "chr3": 90900000, "chr4": 50000000, "chr5": 48800000, "chr6": 59800000,
                "chr7": 60100000, "chr8": 45200000, "chr9": 43000000, "chr10": 39800000, "chr11": 53400000, "chr12": 35500000,
                "chr13": 17700000, "chr14": 17200000, "chr15": 19000000, "chr16": 36800000, "chr17": 25100000, "chr18": 18500000,
                "chr19": 26200000, "chr20": 28100000, "chr21": 12000000, "chr22": 15000000, "chrX": 61000000, "chrY": 10400000}
# the pseudoautosomal regions (GRCh38): X and Y share them, so a male is diploid there
GRCH38_PAR = {"chrX": ((10001, 2781479), (155701383, 156030895)), "chrY": ((10001, 2781479), (56887903, 57217415))}
ACROCENTRIC = ("chr13", "chr14", "chr15", "chr21", "chr22")


@dataclass
class Genome:
    name: str
    length: dict
    p_end: dict
    par: dict
    bands: dict = None          # chrom -> [(start 0-based, end exclusive, name)] from the cytoband table, in order; None without one

    @property
    def chroms(self):
        return list(self.length)

    def band(self, chrom, pos):
        """The cytogenetic band holding a 1-based position ('q22.1'); '' without a table or outside it."""
        if not self.bands or chrom not in self.bands:
            return ""
        for s, e, name in self.bands[chrom]:
            if s < pos <= e:
                return name
        return ""

    def band_range(self, chrom, start, end):
        """The bands spanned by [start, end] (start 0-based, end inclusive) in ISCN form: 'q22.1q31.1', or one band name."""
        a, b = self.band(chrom, start + 1), self.band(chrom, max(end, start + 1))
        if not a or not b:
            return ""
        return a if a == b else a + b

    @property
    def autosomes(self):
        return [c for c in self.length if c not in ("chrX", "chrY", "chrM")]

    def normalize(self, contig):
        """'chr1', '1', 'X', 'chrx' -> 'chr1' ...; None for a contig the genome does not carry (alts, decoys, chrM)."""
        c = str(contig)
        if not c.startswith("chr"):
            c = "chr" + c
        c = c[:3] + c[3:].upper() if c[3:].upper() in ("X", "Y", "M", "MT") else c
        if c == "chrMT":
            c = "chrM"
        return c if c in self.length else None

    def is_par(self, chrom, pos):
        return any(a <= pos <= b for a, b in self.par.get(chrom, ()))

    def arm(self, chrom, pos):
        return "p" if pos <= self.p_end.get(chrom, 0) else "q"

    def arm_bounds(self, chrom, arm):
        if arm == "p":
            return 1, self.p_end[chrom]
        return self.p_end[chrom] + 1, self.length[chrom]


def load_bands(path):
    """{chrom: [(start, end, band)]} from a cytoband table (chrom, start, end, band, stain; '#' comments and a header allowed)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom\t") or not line.strip():
                continue
            c, s, e, b = line.rstrip("\n").split("\t")[:4]
            out.setdefault(c, []).append((int(s), int(e), b))
    for c in out:
        out[c].sort()
    return out


def genome(name="grch38"):
    if name.lower() in ("grch38", "hg38"):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "cytoband.hg38.tsv")
        return Genome("GRCh38", dict(GRCH38_LENGTH), dict(GRCH38_P_END), dict(GRCH38_PAR), load_bands(path) if os.path.exists(path) else None)
    raise ValueError("unknown genome %r (grch38 is built in; pass --lengths/--arms for another)" % name)
