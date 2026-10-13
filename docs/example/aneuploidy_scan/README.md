# Whole-chromosome aneuploidies across the 1000 Genomes DRAGEN cohort

A reproducible scan for trios worth analysing: which of the 3,202 genomes of Illumina's DRAGEN 3.7.6 re-analysis of the
1000 Genomes Project carry a whole-chromosome deviation, which of those are children of a trio in the cohort, and which
have a long-read assembly (HPRC year 1, HPRC release 2, HGSVC3) against which a trio reading of a trisomy can be held.

```bash
docs/example/aneuploidy_scan/fetch.sh work/      # about 6 MB of public downloads, then scan.py -> this directory
```

**Resource.** DRAGEN writes a per-sample ploidy VCF (`<sample>.ploidy.vcf.gz`, about 2 kB) with one record per chromosome:
the depth of coverage `DC`, the depth normalised to the ploidy DRAGEN assumed `NDC`, and `<DUP>` or `<DEL>` where it called
the chromosome's copy number off. `scan.py` reads the copies of every chromosome as 2 DC over the autosomal median DC, so the
X and the Y are on one scale whatever sex was assumed, and flags an autosome at 2.5 copies or more (a trisomy in at least
half the cells), at 2.2 to 2.5 (a mosaic gain in at least a fifth), at 1.5 or under (a monosomy) or 1.5 to 1.8 (a mosaic
loss), and a sex-chromosome complement other than the pedigree sex's, or one with the X or the Y a quarter copy or more off
an integer (a mosaic). The pedigree is the 1000 Genomes 3,202-sample table with populations; the assembly lists are the HPRC
year-1 and release-2 assembly indexes and the HGSVC3 assembly accessions table. Four samples are absent from the DRAGEN
bucket (HG03781, NA07357, NA18954, NA19377).

**Result** (`summary.md`, `aneuploidies.tsv`, `candidates.tsv`). Of 3,198 genomes, 109 deviate: 11 autosomal trisomies, 43
mosaic autosomal gains and 58 sex-chromosome complements other than the pedigree's. The trisomies are the ones
lymphoblastoid lines acquire (six trisomy 12, two trisomy 9, one each of 11 and 13 together, 14 and 15), in two children
(NA12739, +12; NA19208, +9) and five parents. Most mosaic gains sit on the GC-rich chromosomes 19 and 22 at 2.2 to 2.3
copies and DRAGEN itself calls many of them `DUP`; a library's GC bias lifts those chromosomes by a few per cent, so they are
marked as such and left to the trio analysis, whose GC correction decides. The sex-chromosome deviations are mostly mosaic
loss of an X in women's lines and of the Y in men's, with a constitutional 47,XXY (HG01683), 47,XYY (HG03456, HG02966 mosaic)
and two 47,XXX among them, and one pedigree male reading XX (HG02300: the HPRC release-2 sample table records her as female, so the 1000 Genomes pedigree's sex is the error, which the trio's sex check reports).

**No autosomal trisomy has an assembly**: the assembled genomes were karyotypically normal lines, or were selected to be. The
four aneuploid genomes with an assembly are sex-chromosomal: HG03456 (47,XYY; HGSVC3 verkko), HG02300 (XX in a pedigree male;
HPRC release 2 hifiasm trio), and two mothers with a mosaic X loss (HG00732, the mother of HG00733, HGSVC3; HG00864, HGSVC3).
Two of them are children of a trio in the cohort: HG03456 and HG02300, the candidates.

**Trios taken to the example pages** from this scan: NA19208 (trisomy 9 in the child), HG01683 (47,XXY), HG03456 (47,XYY,
assembled), HG02300 (XX in a pedigree male, assembled), HG01358 (the mother HG01357 carries a trisomy 12; the child has an
HPRC year-1 assembly) and NA10865 (the father NA11891 carries a trisomy 15). None of their six members is in the shipped
panel.

## The assemblies

`assembly.sh` holds the commands; `assembly_sexchrom.tsv` the result. For each haplotype of a diploid assembly and each
chromosome, the reference bases its contigs align to (primary and supplementary alignments at mapping quality 20 or
more, the pseudoautosomal regions excluded) and the bases covered once: an extra copy collapsed into one haplotype
gives 1× aligned and 1× covered in that haplotype and nothing in the other; two copies assembled separately give either
2× aligned over 1× covered in one haplotype, or 1× in each.

**HG03456, 47,XYY, HGSVC3 verkko (Strand-seq phased).** HGSVC3 ships the contigs aligned to GRCh38, which samtools reads
remotely one chromosome at a time. The X is one contig of 151 Mb in haplotype 2 (150 Mb of the 153 Mb outside the
PARs covered, 1.01× aligned over covered) and nothing in haplotype 1. The Y is one contig of 36 Mb in haplotype 1
(23.4 Mb covered, 1.59× aligned over covered: the Y's palindromes and ampliconic arrays map many-to-one) and 0.5 Mb
of PAR-adjacent fragments in haplotype 2. The control, HG00731 (46,XY) from the same batch, reads the same: one X
contig in one haplotype (150 Mb covered, 1.02×) and one Y contig of 34 Mb in the other (23.5 Mb covered, 1.46×).
The two Y chromosomes of HG03456 are identical copies, from a paternal meiosis II or a post-zygotic error, so the
assembly graph has nothing to separate them on and verkko collapsed them into one path: the assembly of the 47,XYY
genome reads as 46,XY. The second Y is visible only by dosage, which the trio analysis reads from the Y depth and
from the father/son Y depth ratio.

**HG02300, XX in a pedigree male, HPRC release 2 (hifiasm, trio-binned).** The two haplotype FASTAs aligned to the
GRCh38 X and Y with minimap2 (asm5). Each haplotype covers the X once: 149.7 Mb of the 152.9 Mb outside the PARs in the
paternal haplotype, 150.7 Mb in the maternal, with contigs of 67 and 69 Mb; neither holds a Y, the 9.4 Mb of Y "coverage"
per haplotype being satellite and PAR-boundary repeats at Yq11, Yq12 and the PAR2 edge (five times aligned over
covered) rather than the Y's euchromatin, which a real Y gives as one contig covering 23 Mb (HG00731 above). A paternal
X and a maternal X, no Y: a daughter, as HPRC's sample table records, and the 1000 Genomes pedigree has the sex wrong.
The aligned-over-covered of 1.9 on the X is an artefact of the GenBank release holding each chromosome both as a
scaffold (CM accessions) and as its contigs (JBHD accessions); the covered bases are the reading.
