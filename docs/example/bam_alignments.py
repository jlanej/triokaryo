#!/usr/bin/env python3
"""SAM on stdin (samtools view of contigs aligned to a reference) to the alignment table assembly_sexchrom.py reads: hap, chrom,
contig, start (1-based), end, reference span, mapq, kind (primary / secondary / supplementary).

  samtools view URL_or_BAM chrY | bam_alignments.py hap1 >> alignments.tsv"""
import re
import sys

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
label = sys.argv[1] if len(sys.argv) > 1 else "hap"
for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.split("\t", 9)
    flag, pos = int(f[1]), int(f[3])
    span = sum(int(n) for n, op in CIG.findall(f[5]) if op in "MDN=X")
    kind = "secondary" if flag & 256 else "supplementary" if flag & 2048 else "primary"
    sys.stdout.write("%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n" % (label, f[2], f[0], pos, pos + span, span, f[4], kind))
